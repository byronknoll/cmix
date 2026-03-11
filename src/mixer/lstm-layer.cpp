#include "lstm-layer.h"

#include "sigmoid.h"

#include <math.h>
#include <algorithm>
#include <numeric>

#ifdef __AVX2__
#include <immintrin.h>

namespace {

inline float HorizontalSum(__m256 v) {
  __m128 hi = _mm256_extractf128_ps(v, 1);
  __m128 lo = _mm256_castps256_ps128(v);
  lo = _mm_add_ps(lo, hi);
  hi = _mm_movehl_ps(hi, lo);
  lo = _mm_add_ps(lo, hi);
  hi = _mm_shuffle_ps(lo, lo, 1);
  lo = _mm_add_ss(lo, hi);
  return _mm_cvtss_f32(lo);
}

inline float AvxDotProduct(const float* a, const float* b, int n) {
  __m256 sum = _mm256_setzero_ps();
  int i = 0;
  for (; i + 7 < n; i += 8) {
    __m256 va = _mm256_loadu_ps(a + i);
    __m256 vb = _mm256_loadu_ps(b + i);
    sum = _mm256_fmadd_ps(va, vb, sum);
  }
  float result = HorizontalSum(sum);
  for (; i < n; ++i) {
    result += a[i] * b[i];
  }
  return result;
}

void Adam(std::valarray<float>* g, std::valarray<float>* m,
    std::valarray<float>* v, std::valarray<float>* w, float learning_rate,
    float t, unsigned long long update_limit) {
  const float beta1 = 0.025f, beta2 = 0.9999f, eps = 1e-6f;
  float alpha;
  if (t < update_limit) {
    alpha = learning_rate * 0.1f / sqrt(5e-5f * t + 1.0f);
  } else {
    alpha = learning_rate * 0.1f / sqrt(5e-5f * update_limit + 1.0f);
  }

  float bias_m, bias_v_corr;
  if (t < update_limit) {
    bias_m = 1.0f / (1.0f - pow(beta1, t));
    bias_v_corr = 1.0f / (1.0f - pow(beta2, t));
  } else {
    bias_m = 1.0f / (1.0f - pow(beta1, update_limit));
    bias_v_corr = 1.0f / (1.0f - pow(beta2, update_limit));
  }

  int n = g->size();
  float* gp = &(*g)[0];
  float* mp = &(*m)[0];
  float* vp = &(*v)[0];
  float* wp = &(*w)[0];

  __m256 vbeta1 = _mm256_set1_ps(beta1);
  __m256 v1mbeta1 = _mm256_set1_ps(1.0f - beta1);
  __m256 vbeta2 = _mm256_set1_ps(beta2);
  __m256 v1mbeta2 = _mm256_set1_ps(1.0f - beta2);
  __m256 valpha = _mm256_set1_ps(alpha);
  __m256 vbias_m = _mm256_set1_ps(bias_m);
  __m256 vbias_v = _mm256_set1_ps(bias_v_corr);
  __m256 veps = _mm256_set1_ps(eps);

  int i = 0;
  for (; i + 7 < n; i += 8) {
    __m256 vg = _mm256_loadu_ps(gp + i);
    __m256 vm = _mm256_loadu_ps(mp + i);
    __m256 vv = _mm256_loadu_ps(vp + i);
    __m256 vw = _mm256_loadu_ps(wp + i);

    // m = m * beta1 + (1 - beta1) * g
    vm = _mm256_fmadd_ps(vm, vbeta1, _mm256_mul_ps(v1mbeta1, vg));
    // v = v * beta2 + (1 - beta2) * g * g
    vv = _mm256_fmadd_ps(vv, vbeta2, _mm256_mul_ps(v1mbeta2, _mm256_mul_ps(vg, vg)));
    // w -= alpha * (m * bias_m) / sqrt(v * bias_v + eps)
    __m256 m_hat = _mm256_mul_ps(vm, vbias_m);
    __m256 v_hat = _mm256_fmadd_ps(vv, vbias_v, veps);
    __m256 update = _mm256_mul_ps(valpha, _mm256_div_ps(m_hat, _mm256_sqrt_ps(v_hat)));
    vw = _mm256_sub_ps(vw, update);

    _mm256_storeu_ps(mp + i, vm);
    _mm256_storeu_ps(vp + i, vv);
    _mm256_storeu_ps(wp + i, vw);
  }
  // Scalar tail
  for (; i < n; ++i) {
    mp[i] = mp[i] * beta1 + (1.0f - beta1) * gp[i];
    vp[i] = vp[i] * beta2 + (1.0f - beta2) * gp[i] * gp[i];
    wp[i] -= alpha * (mp[i] * bias_m) / sqrt(vp[i] * bias_v_corr + eps);
  }
}

}  // namespace

#else

namespace {

void Adam(std::valarray<float>* g, std::valarray<float>* m,
    std::valarray<float>* v, std::valarray<float>* w, float learning_rate,
    float t, unsigned long long update_limit) {
  const float beta1 = 0.025, beta2 = 0.9999, eps = 1e-6f;
  float alpha;
  if (t < update_limit) {
    alpha = learning_rate * 0.1f / sqrt(5e-5f * t + 1.0f);
  } else {
    alpha = learning_rate * 0.1f / sqrt(5e-5f * update_limit + 1.0f);
  }
  (*m) *= beta1;
  (*m) += (1.0f - beta1) * (*g);
  (*v) *= beta2;
  (*v) += (1.0f - beta2) * (*g) * (*g);
  if (t < update_limit) {
    (*w) -= alpha * (((*m) / (float)(1.0f - pow(beta1, t))) /
        (sqrt((*v) / (float)(1.0f - pow(beta2, t)) + eps)));
  } else {
    (*w) -= alpha * (((*m) / (float)(1.0f - pow(beta1, update_limit))) /
        (sqrt((*v) / (float)(1.0f - pow(beta2, update_limit)) + eps)));
  }
}

}  // namespace
#endif

LstmLayer::LstmLayer(unsigned int input_size, unsigned int auxiliary_input_size,
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

void LstmLayer::ForwardPass(const std::valarray<float>& input, int input_symbol,
    std::valarray<float>* hidden, int hidden_start) {
  last_state_[epoch_] = state_;
  ForwardPass(forget_gate_, input, input_symbol);
  ForwardPass(input_node_, input, input_symbol);
  ForwardPass(output_gate_, input, input_symbol);
  for (unsigned int i = 0; i < num_cells_; ++i) {
    forget_gate_.state_[epoch_][i] = Sigmoid::Logistic(
        forget_gate_.state_[epoch_][i]);
    input_node_.state_[epoch_][i] = tanh(input_node_.state_[epoch_][i]);
    output_gate_.state_[epoch_][i] = Sigmoid::Logistic(
        output_gate_.state_[epoch_][i]);
  }
  input_gate_state_[epoch_] = 1.0f - forget_gate_.state_[epoch_];
  state_ *= forget_gate_.state_[epoch_];
  state_ += input_node_.state_[epoch_] * input_gate_state_[epoch_];
  tanh_state_[epoch_] = tanh(state_);
  std::slice slice = std::slice(hidden_start, num_cells_, 1);
  (*hidden)[slice] = output_gate_.state_[epoch_] * tanh_state_[epoch_];
  ++epoch_;
  if (epoch_ == horizon_) epoch_ = 0;
}

void LstmLayer::ForwardPass(NeuronLayer& neurons,
    const std::valarray<float>& input, int input_symbol) {
#ifdef __AVX2__
  int input_n = input.size();
  const float* input_ptr = &input[0];
  for (unsigned int i = 0; i < num_cells_; ++i) {
    float f = neurons.weights_[i][input_symbol];
    const float* w_ptr = &neurons.weights_[i][output_size_];
    f += AvxDotProduct(input_ptr, w_ptr, input_n);
    neurons.norm_[epoch_][i] = f;
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

void LstmLayer::ClipGradients(std::valarray<float>* arr) {
  for (unsigned int i = 0; i < arr->size(); ++i) {
    if ((*arr)[i] < -gradient_clip_) (*arr)[i] = -gradient_clip_;
    else if ((*arr)[i] > gradient_clip_) (*arr)[i] = gradient_clip_;
  }
}

void LstmLayer::BackwardPass(const std::valarray<float>&input, int epoch,
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
    if (update_steps_ < update_limit_) {
      ++update_steps_;
    }
  }

  BackwardPass(forget_gate_, input, epoch, layer, input_symbol, hidden_error);
  BackwardPass(input_node_, input, epoch, layer, input_symbol, hidden_error);
  BackwardPass(output_gate_, input, epoch, layer, input_symbol, hidden_error);

  ClipGradients(&state_error_);
  ClipGradients(&stored_error_);
  ClipGradients(hidden_error);
}

void LstmLayer::BackwardPass(NeuronLayer& neurons,
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

#ifdef __AVX2__
  const float* err_ptr = &neurons.error_[0];
  int nc = num_cells_;

  if (layer > 0) {
    for (unsigned int i = 0; i < num_cells_; ++i) {
      const float* t_ptr = &neurons.transpose_[num_cells_ + i][0];
      (*hidden_error)[i] += AvxDotProduct(err_ptr, t_ptr, nc);
    }
  }
  if (epoch > 0) {
    for (unsigned int i = 0; i < num_cells_; ++i) {
      const float* t_ptr = &neurons.transpose_[i][0];
      stored_error_[i] += AvxDotProduct(err_ptr, t_ptr, nc);
    }
  }
#else
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
#endif

  std::slice slice = std::slice(output_size_, input.size(), 1);
  for (unsigned int i = 0; i < num_cells_; ++i) {
    neurons.update_[i][slice] += neurons.error_[i] * input;
    neurons.update_[i][input_symbol] += neurons.error_[i];
  }
  if (epoch == 0) {
    for (unsigned int i = 0; i < num_cells_; ++i) {
      Adam(&neurons.update_[i], &neurons.m_[i], &neurons.v_[i],
          &neurons.weights_[i], learning_rate_, update_steps_, update_limit_);
    }
    Adam(&neurons.gamma_u_, &neurons.gamma_m_, &neurons.gamma_v_,
        &neurons.gamma_, learning_rate_, update_steps_, update_limit_);
    Adam(&neurons.beta_u_, &neurons.beta_m_, &neurons.beta_v_,
        &neurons.beta_, learning_rate_, update_steps_, update_limit_);
  }
}

std::vector<std::valarray<std::valarray<float>>*> LstmLayer::Weights() {
  std::vector<std::valarray<std::valarray<float>>*> weights;
  weights.push_back(&forget_gate_.weights_);
  weights.push_back(&input_node_.weights_);
  weights.push_back(&output_gate_.weights_);
  return weights;
}
