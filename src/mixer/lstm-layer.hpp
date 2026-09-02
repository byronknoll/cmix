#include "lstm-layer.h"

#include "sigmoid.h"

#include <math.h>
#include <algorithm>
#include <numeric>

#ifdef P2_FASTMATH
#include "fastmath.h"
#define FAST_TANH fast_tanh
#define FAST_TANH_VEC fast_tanh_vec
#define LSTM_LOGISTIC fast_logistic
#define LSTM_EXPF fast_expf
#else
#define FAST_TANH tanh //fast_tanh
#define FAST_TANH_VEC tanh //fast_tanh_vec
#define LSTM_LOGISTIC Sigmoid::Logistic
#define LSTM_EXPF exp
#endif
namespace {
// inline float fast_tanh(const float x)
// {
//     const float ax = fabs(x);
//     const float x2 = x * x;

//     return(x * (2.45550750702956f + 2.45550750702956f * ax +
//         (0.893229853513558f + 0.821226666969744f * ax) * x2) /
//         (2.44506634652299f + (2.44506634652299f + x2) *
//             fabs(x + 0.814642734961073f * x * ax)));
// }

// float fast_tanh(float x){
//   float x2 = x * x;
//   float a = x * (135135.0f + x2 * (17325.0f + x2 * (378.0f + x2)));
//   float b = 135135.0f + x2 * (62370.0f + x2 * (3150.0f + x2 * 28.0f));
//   return a / b;
// }
// 
// template <class _Tp>
// inline std::valarray<_Tp> fast_tanh_vec(const std::valarray<_Tp>& __x) {
//   std::valarray<_Tp> __tmp(__x.size());
//   for (size_t __i = 0; __i < __x.size(); ++__i)
//     __tmp[__i] = fast_tanh(__x[__i]);
//   return __tmp;
// }

// Optimizer-constant knobs (code-dive 2026-07-11), all flag-off identical and all
// passed via CFLAGS_DEFINES (reaches all three flag groups; OLR_RAMP in lstm.h adds
// members, so the CFLAGS_DEFINES route is mandatory there — same ODR rule as GM_*).
// LSTM_BETA1/LSTM_BETA2: hoist the inherited Adam moments (0.025/0.9999, byte-identical
//   since >=2016, zero arms ever run). beta1=0.9 averages int8 bwd-quant spikes over
//   ~10 cycles (15.7-16.1x quiet-weight shock reduction) = decoupler candidate for the
//   confirmed OLR x quant overshoot. Two-sided: EG1 says full-cadence noise HELPS.
// LSTM_ALPHA_CAP_ONLY: decouple warmup length from plateau level — alpha keeps the
//   UPDATE_LIMIT cap but the pow bias-corrections track true t (plateau 2.80e-3 =
//   1.96x stock, 1.40x above what the coupled knob can express; float t exact to
//   2^24 so 1GB's t=7.8e6 is safe; pow(beta,large t) underflows to 0 -> correction 1).
// LSTM_EPS_OUTSIDE: standard sqrt(v)+eps placement (the inside placement floors the
//   denominator at 1e-3 = accidental plain SGD on small-gradient late-stream
//   coordinates). De-confound CONTROL for schedule arms, not an expected keep.
#ifndef LSTM_BETA1
#define LSTM_BETA1 0.025
#endif
#ifndef LSTM_BETA2
#define LSTM_BETA2 0.9999
#endif
#ifndef UPDATE_LIMIT
#define UPDATE_LIMIT 3000
#endif
#ifdef LSTM_EPS_OUTSIDE
#define LSTM_ADAM_DENOM(vhat) (sqrt(vhat) + 1e-8f)
#else
#define LSTM_ADAM_DENOM(vhat) (sqrt((vhat) + eps))
#endif
inline void Adam(std::valarray<float>* g, std::valarray<float>* m,
    std::valarray<float>* v, std::valarray<float>* w, float learning_rate,
    float t) {
  const float beta1 = LSTM_BETA1, beta2 = LSTM_BETA2, eps = 1e-6f;
  (void)eps;
  float alpha;
  if (t < UPDATE_LIMIT) {
    alpha = learning_rate * 0.1f / sqrt(5e-5f * t + 1.0f);
  } else {
    alpha = learning_rate * 0.1f / sqrt(5e-5f * UPDATE_LIMIT + 1.0f);
  }
  (*m) *= beta1;
  (*m) += (1.0f - beta1) * (*g);
  (*v) *= beta2;
  (*v) += (1.0f - beta2) * (*g) * (*g);
#ifdef LSTM_ALPHA_CAP_ONLY
  // bias corrections track true t in both regimes; only alpha stays capped
  (*w) -= alpha * (((*m) / (float)(1.0f - pow(beta1, t))) /
      LSTM_ADAM_DENOM((*v) / (float)(1.0f - pow(beta2, t))));
#else
  if (t < UPDATE_LIMIT) {
    (*w) -= alpha * (((*m) / (float)(1.0f - pow(beta1, t))) /
        LSTM_ADAM_DENOM((*v) / (float)(1.0f - pow(beta2, t))));
  } else {
    (*w) -= alpha * (((*m) / (float)(1.0f - pow(beta1, UPDATE_LIMIT))) /
        LSTM_ADAM_DENOM((*v) / (float)(1.0f - pow(beta2, UPDATE_LIMIT))));
  }
#endif
}

}

inline LstmLayer::LstmLayer(unsigned int input_size, unsigned int auxiliary_input_size,
    unsigned int output_size, unsigned int num_cells, int horizon,
    float gradient_clip, float learning_rate) :
    state_(num_cells), state_error_(num_cells), stored_error_(num_cells),
    tanh_state_(std::valarray<float>(num_cells), horizon),
    input_gate_state_(std::valarray<float>(num_cells), horizon),
    last_state_(std::valarray<float>(num_cells), horizon),
    gradient_clip_(gradient_clip), learning_rate_(learning_rate),
    num_cells_(num_cells), epoch_(0), horizon_(horizon),
    input_size_(auxiliary_input_size), output_size_(output_size),
    forget_gate_(input_size, num_cells, horizon, output_size_ + input_size_),
    input_node_(input_size, num_cells, horizon, output_size_ + input_size_),
    output_gate_(input_size, num_cells, horizon, output_size_ + input_size_) {
#ifdef P3_BATCH
  input_hist_.resize(horizon, nullptr);
  sym_hist_.resize(horizon, 0);
#endif
#ifdef LSTM_QUANT
  // input_size (ctor arg) = output_size embed columns + matvec columns; the
  // last matvec column multiplies the constant-1 bias input and stays float.
  qcols_ = input_size - output_size - 1;
  qstride_ = (qcols_ + (LSTM_QUANT_PAD - 1)) & ~(LSTM_QUANT_PAD - 1);
  if (qcols_ == 0 || qstride_ > 65536u) {
    fprintf(stderr, "LSTM_QUANT: unsupported width qcols=%u\n", qcols_);
    abort();  // int32 accumulator headroom audited for n <= 65536 only
  }
  qinput_.assign(qstride_, 0);
#ifdef LSTM_QUANT_KERNEL_VNNI
  qinput_biased_.assign(qstride_, 128);
#endif
  NeuronLayer* gates[3] = {&forget_gate_, &input_node_, &output_gate_};
  for (int g = 0; g < 3; ++g) {
    gates[g]->qweights_.assign((size_t)num_cells_ * qstride_, 0);
    gates[g]->qscale_.assign(num_cells_, 1.0f);
    gates[g]->qrowsum_.assign(num_cells_, 0);
  }
#endif
  float val = sqrt(6.0f / float(input_size_ + output_size_));
  float low = -val;
  float range = 2 * val;
  for (unsigned int i = 0; i < num_cells_; ++i) {
    for (unsigned int j = 0; j < forget_gate_.weights_[i].size(); ++j) {
      forget_gate_.weights_[i][j] = low + Rand() * range;
      input_node_.weights_[i][j] = low + Rand() * range;
      output_gate_.weights_[i][j] = low + Rand() * range;
    }
    forget_gate_.weights_[i][forget_gate_.weights_[i].size() - 1] = 1;
  }
}

inline void LstmLayer::ForwardPass(const std::valarray<float>& input, int input_symbol,
    std::valarray<float>* hidden, int hidden_start) {
  last_state_[epoch_] = state_;
#ifdef LSTM_QUANT
  // Weights change only in Adam (epoch 0 of the BPTT cycle): requantize
  // lazily, once per horizon steps. Activations are quantized once per step
  // and shared by all three gates.
  if (qdirty_) {
    QuantRequantize();
    qdirty_ = false;
  }
  QuantizeInput(input);
#endif
  ForwardPass(forget_gate_, input, input_symbol);
  ForwardPass(input_node_, input, input_symbol);
  ForwardPass(output_gate_, input, input_symbol);
  for (unsigned int i = 0; i < num_cells_; ++i) {
    forget_gate_.state_[epoch_][i] = LSTM_LOGISTIC(
        forget_gate_.state_[epoch_][i]);
    input_node_.state_[epoch_][i] = FAST_TANH(input_node_.state_[epoch_][i]);
    output_gate_.state_[epoch_][i] = LSTM_LOGISTIC(
        output_gate_.state_[epoch_][i]);
  }
  input_gate_state_[epoch_] = 1.0f - forget_gate_.state_[epoch_];
  state_ *= forget_gate_.state_[epoch_];
  state_ += input_node_.state_[epoch_] * input_gate_state_[epoch_];
  tanh_state_[epoch_] = FAST_TANH_VEC(state_);
  std::slice slice = std::slice(hidden_start, num_cells_, 1);
  (*hidden)[slice] = output_gate_.state_[epoch_] * tanh_state_[epoch_];
  ++epoch_;
  if (epoch_ == horizon_) epoch_ = 0;
}

inline void LstmLayer::ForwardPass(NeuronLayer& neurons,
    const std::valarray<float>& input, int input_symbol) {
#ifdef LSTM_QUANT
  // f = embed + bias + acc * (row_scale * act_scale); acc is the exact
  // integer dot of the quantized row and activations (deterministic:
  // integer accumulation, no FMA/reassociation ambiguity).
  const unsigned int bias_col = (unsigned int)neurons.weights_[0].size() - 1;
  for (unsigned int i = 0; i < num_cells_; ++i) {
    const LstmQuantAcc acc = QuantDotRow(neurons, i);
    neurons.norm_[epoch_][i] = neurons.weights_[i][input_symbol] +
        neurons.weights_[i][bias_col] +
        (float)acc * (neurons.qscale_[i] * qact_scale_);
  }
#else
  for (unsigned int i = 0; i < num_cells_; ++i) {
    float f = neurons.weights_[i][input_symbol];
    for (unsigned int j = 0; j < input.size(); ++j) {
      f += input[j] * neurons.weights_[i][output_size_ + j];
    }
    neurons.norm_[epoch_][i] = f;
  }
#endif
  neurons.ivar_[epoch_] = 1.0f / sqrt(((neurons.norm_[epoch_] *
      neurons.norm_[epoch_]).sum() / num_cells_) + 1e-5f);
  neurons.norm_[epoch_] *= neurons.ivar_[epoch_];
  neurons.state_[epoch_] = neurons.norm_[epoch_] * neurons.gamma_ +
      neurons.beta_;
}

inline void LstmLayer::ClipGradients(std::valarray<float>* arr) {
  for (unsigned int i = 0; i < arr->size(); ++i) {
    if ((*arr)[i] < -gradient_clip_) (*arr)[i] = -gradient_clip_;
    else if ((*arr)[i] > gradient_clip_) (*arr)[i] = gradient_clip_;
  }
}

inline void LstmLayer::BackwardPass(const std::valarray<float>&input, int epoch,
    int layer, int input_symbol, std::valarray<float>* hidden_error) {
  if (epoch == (int)horizon_ - 1) {
    stored_error_ = *hidden_error;
    state_error_ = 0;
  } else {
    stored_error_ += *hidden_error;
  }

  output_gate_.error_ = tanh_state_[epoch] * stored_error_ *
      output_gate_.state_[epoch] * (1.0f - output_gate_.state_[epoch]);
  state_error_ += stored_error_ * output_gate_.state_[epoch] * (1.0f -
      (tanh_state_[epoch] * tanh_state_[epoch]));
  input_node_.error_ = state_error_ * input_gate_state_[epoch] * (1.0f -
      (input_node_.state_[epoch] * input_node_.state_[epoch]));
  forget_gate_.error_ = (last_state_[epoch] - input_node_.state_[epoch]) *
      state_error_ * forget_gate_.state_[epoch] * input_gate_state_[epoch];

  *hidden_error = 0;
  if (epoch > 0) {
    state_error_ *= forget_gate_.state_[epoch];
    stored_error_ = 0;
  } else {
#ifdef LSTM_ALPHA_CAP_ONLY
    // true t must keep counting; alpha's cap lives inside Adam() under this knob
    ++update_steps_;
#else
    if (update_steps_ < UPDATE_LIMIT) {
      ++update_steps_;
    }
#endif
  }

#ifdef P3_BATCH
  input_hist_[epoch] = &input;
  sym_hist_[epoch] = input_symbol;
#endif
  BackwardPass(forget_gate_, input, epoch, layer, input_symbol, hidden_error);
  BackwardPass(input_node_, input, epoch, layer, input_symbol, hidden_error);
  BackwardPass(output_gate_, input, epoch, layer, input_symbol, hidden_error);

#ifdef LSTM_QUANT
  // Adam just ran inside the three per-gate passes: the float master
  // weights_ changed, so the quantized shadow is stale.
  if (epoch == 0) qdirty_ = true;
#endif

  ClipGradients(&state_error_);
  ClipGradients(&stored_error_);
  ClipGradients(hidden_error);
}

inline void LstmLayer::BackwardPass(NeuronLayer& neurons,
    const std::valarray<float>&input, int epoch, int layer, int input_symbol,
    std::valarray<float>* hidden_error) {
  if (epoch == (int)horizon_ - 1) {
    neurons.gamma_u_ = 0;
    neurons.beta_u_ = 0;
    for (unsigned int i = 0; i < num_cells_; ++i) {
      neurons.update_[i] = 0;
      int offset = output_size_ + input_size_;
      for (unsigned int j = 0; j < neurons.transpose_.size(); ++j) {
        neurons.transpose_[j][i] = neurons.weights_[i][j + offset];
      }
    }
  }
  neurons.beta_u_ += neurons.error_;
  neurons.gamma_u_ += neurons.error_ * neurons.norm_[epoch];
  neurons.error_ *= neurons.gamma_ * neurons.ivar_[epoch];
  neurons.error_ -= ((neurons.error_ * neurons.norm_[epoch]).sum() /
      num_cells_) * neurons.norm_[epoch];
  if (layer > 0) {
    for (unsigned int i = 0; i < num_cells_; ++i) {
      float f = 0;
      for (unsigned int j = 0; j < num_cells_; ++j) {
        f += neurons.error_[j] * neurons.transpose_[num_cells_ + i][j];
      }
      (*hidden_error)[i] += f;
    }
  }
  if (epoch > 0) {
    for (unsigned int i = 0; i < num_cells_; ++i) {
      float f = 0;
      for (unsigned int j = 0; j < num_cells_; ++j) {
        f += neurons.error_[j] * neurons.transpose_[i][j];
      }
      stored_error_[i] += f;
    }
  }
#ifndef P3_BATCH
  std::slice slice = std::slice(output_size_, input.size(), 1);
  for (unsigned int i = 0; i < num_cells_; ++i) {
    neurons.update_[i][slice] += neurons.error_[i] * input;
    neurons.update_[i][input_symbol] += neurons.error_[i];
  }
#else
  // Defer the rank-1 update accumulation: store this epoch's error vector,
  // flush the whole cycle in one cache-blocked pass at epoch 0. Per-element
  // accumulation order (epoch descending) matches the eager version exactly,
  // so the result is bit-identical; only memory traffic changes (update_
  // streaming drops from ~1.5MB/byte to ~16KB/byte).
  neurons.error_hist_[epoch] = neurons.error_;
  if (epoch == 0) {
#if   defined(P3_REGBLOCK)
    // Register-blocked cycle flush. Byte-identical to the batched path (MD5-
    // gated): weight-block OUTER, epoch INNER, so each block's update
    // accumulators live in registers across ALL horizon epochs instead of the
    // update_ row being loaded+stored once per epoch (~horizon-fold less
    // update_ load/store traffic). Per-weight accumulation stays epoch-
    // descending; the embed-column scatter writes the DISJOINT region
    // [0,output_size_), so it is hoisted out with its epoch order preserved =>
    // identical result bytes. FMA contraction is the same as the batched loop
    // (both `acc += fe*ip` under -ffp-model=fast).
    for (unsigned int i = 0; i < num_cells_; ++i) {
      float* up = &neurons.update_[i][0];
      float* upd = up + output_size_;
      for (int e = horizon_ - 1; e >= 0; --e)
        up[sym_hist_[e]] += neurons.error_hist_[e][i];
      const unsigned int n = (unsigned int)input_hist_[0]->size();
      constexpr unsigned int BLK = 32u;
      for (unsigned int jb = 0; jb < n; jb += BLK) {
        const unsigned int bn = (jb + BLK <= n) ? BLK : (n - jb);
        float acc[BLK];
        for (unsigned int k = 0; k < bn; ++k) acc[k] = upd[jb + k];
        for (int e = horizon_ - 1; e >= 0; --e) {
          const float fe = neurons.error_hist_[e][i];
          const float* ip = std::begin(*input_hist_[e]) + jb;
          for (unsigned int k = 0; k < bn; ++k) acc[k] += fe * ip[k];
        }
        for (unsigned int k = 0; k < bn; ++k) upd[jb + k] = acc[k];
      }
    }
#elif defined(P3_MICROKERNEL) && defined(__AVX512F__)
    // #1: true register-tiled GEMM microkernel for update_ += E^T @ X (C=E^T@X,
    // C=update_ CxN, E=error_hist_ HxC, X=input_hist HxN). Tiles MR cells; per
    // 16-col block keeps MR zmm accumulators in registers across ALL epochs
    // (cross-epoch update_ reuse, like REGBLOCK) AND loads X[e][block] ONCE per
    // epoch, reusing it across the MR cells (cross-cell reuse — the axis REGBLOCK
    // missed). Per-weight accumulation stays epoch-DESCENDING and contracts to the
    // same fma(fe,x,acc) => bit-identical to the batched path (MD5-gated). The
    // embed-column scatter [0,output_size_) is disjoint from the input columns, so
    // it is hoisted out with its per-cell epoch-descending order preserved.
    {
      const unsigned int H = (unsigned int)horizon_;
      const unsigned int n = (unsigned int)input_hist_[0]->size();
      const float* Xp[512]; const float* Ep[512];
      for (unsigned int e = 0; e < H; ++e) {
        Xp[e] = std::begin(*input_hist_[e]);
        Ep[e] = std::begin(neurons.error_hist_[e]);
      }
      for (unsigned int i = 0; i < num_cells_; ++i) {
        float* up = &neurons.update_[i][0];
        for (int e = (int)H - 1; e >= 0; --e) up[sym_hist_[e]] += Ep[e][i];
      }
      constexpr unsigned int MR = 8;
      unsigned int ib = 0;
      for (; ib + MR <= num_cells_; ib += MR) {
        float* upd[MR];
        for (unsigned int c = 0; c < MR; ++c)
          upd[c] = &neurons.update_[ib + c][output_size_];
        unsigned int j = 0;
        for (; j + 16 <= n; j += 16) {
          __m512 acc[MR];
          for (unsigned int c = 0; c < MR; ++c) acc[c] = _mm512_loadu_ps(upd[c] + j);
          for (int e = (int)H - 1; e >= 0; --e) {
            const __m512 x = _mm512_loadu_ps(Xp[e] + j);
            for (unsigned int c = 0; c < MR; ++c)
              acc[c] = _mm512_fmadd_ps(_mm512_set1_ps(Ep[e][ib + c]), x, acc[c]);
          }
          for (unsigned int c = 0; c < MR; ++c) _mm512_storeu_ps(upd[c] + j, acc[c]);
        }
        for (; j < n; ++j)
          for (int e = (int)H - 1; e >= 0; --e) {
            const float xj = Xp[e][j];
            for (unsigned int c = 0; c < MR; ++c) upd[c][j] += Ep[e][ib + c] * xj;
          }
      }
      for (; ib < num_cells_; ++ib) {
        float* upd = &neurons.update_[ib][output_size_];
        for (int e = (int)H - 1; e >= 0; --e) {
          const float fe = Ep[e][ib];
          const float* ip = Xp[e];
          for (unsigned int j = 0; j < n; ++j) upd[j] += fe * ip[j];
        }
      }
    }
#else
    for (unsigned int i = 0; i < num_cells_; ++i) {
      float* up = &neurons.update_[i][0];
      float* upd = up + output_size_;
      for (int e = horizon_ - 1; e >= 0; --e) {
        const float fe = neurons.error_hist_[e][i];
        const std::valarray<float>& in_e = *input_hist_[e];
        const float* ip = std::begin(in_e);
        const unsigned int n = in_e.size();
        for (unsigned int j = 0; j < n; ++j) upd[j] += fe * ip[j];
        up[sym_hist_[e]] += fe;
      }
    }
#endif
  }
#endif
  if (epoch == 0) {
    for (unsigned int i = 0; i < num_cells_; ++i) {
      Adam(&neurons.update_[i], &neurons.m_[i], &neurons.v_[i],
          &neurons.weights_[i], learning_rate_, update_steps_);
    }
    Adam(&neurons.gamma_u_, &neurons.gamma_m_, &neurons.gamma_v_,
        &neurons.gamma_, learning_rate_, update_steps_);
    Adam(&neurons.beta_u_, &neurons.beta_m_, &neurons.beta_v_,
        &neurons.beta_, learning_rate_, update_steps_);
  }
}


#ifdef LSTM_QUANT
inline void LstmLayer::QuantRequantizeGate(NeuronLayer& neurons) {
  for (unsigned int i = 0; i < num_cells_; ++i) {
    // valarray storage is contiguous; matvec block = cols [output_size_,
    // row_size-1) of the float master row.
    neurons.qscale_[i] = LstmQuantizeRow(&neurons.weights_[i][output_size_],
        qcols_, &neurons.qweights_[(size_t)i * qstride_], qstride_,
        &neurons.qrowsum_[i]);
  }
}

inline void LstmLayer::QuantRequantize() {
  QuantRequantizeGate(forget_gate_);
  QuantRequantizeGate(input_node_);
  QuantRequantizeGate(output_gate_);
}

inline void LstmLayer::QuantizeInput(const std::valarray<float>& input) {
  if ((unsigned int)input.size() != qcols_ + 1) {
    fprintf(stderr, "LSTM_QUANT: input width %u != qcols+1 %u\n",
        (unsigned int)input.size(), qcols_ + 1);
    abort();  // NDEBUG kills assert in this codebase; fail loudly instead
  }
  // Quantize all activations except the trailing constant-1 bias input
  // (handled float-exact via its own weight column).
  qact_scale_ = LstmQuantizeActs(&input[0], qcols_, qinput_.data(), qstride_);
#ifdef LSTM_QUANT_KERNEL_VNNI
  // vpdpbusd needs an unsigned operand: bias the activations by +128; the
  // dot kernel subtracts 128*rowsum(qweights). Padding maps 0 -> 128 and is
  // cancelled exactly because padded weights are 0.
  for (unsigned int j = 0; j < qstride_; ++j) {
    qinput_biased_[j] = (uint8_t)((int32_t)qinput_[j] + 128);
  }
#endif
}

inline LstmQuantAcc LstmLayer::QuantDotRow(const NeuronLayer& neurons,
    unsigned int i) const {
#ifdef LSTM_QUANT_KERNEL_VNNI
  return LstmQuantDot(&neurons.qweights_[(size_t)i * qstride_],
      qinput_.data(), qinput_biased_.data(), qstride_, neurons.qrowsum_[i]);
#else
  return LstmQuantDot(&neurons.qweights_[(size_t)i * qstride_],
      qinput_.data(), nullptr, qstride_, neurons.qrowsum_[i]);
#endif
}
#endif  // LSTM_QUANT

inline std::vector<std::valarray<std::valarray<float>>*> LstmLayer::Weights() {
  std::vector<std::valarray<std::valarray<float>>*> weights;
  weights.push_back(&forget_gate_.weights_);
  weights.push_back(&input_node_.weights_);
  weights.push_back(&output_gate_.weights_);
  return weights;
}
