#include "lstm.h"

#include <numeric>
#include <stdlib.h>
#include <fstream>
#include <iostream>

// F-A7: the output layer runs plain SGD with the raw constructor learning
// rate (0.03 from predictor.cpp), decoupled from the Adam schedule inside
// LstmLayer. -DOUTPUT_LR=<float> gives it its own knob; when the macro is
// not defined the code below expands to the stock member read, so the
// flag-off build is token-identical.
#ifdef OUTPUT_LR
#define LSTM_OUTPUT_LR_BASE ((float)(OUTPUT_LR))
#else
#define LSTM_OUTPUT_LR_BASE learning_rate_
#endif

// LSTM_LR_SCALE_OUT=<float> (track 3, per-layer LR, 2026-07-12): a pure
// multiplier on the EFFECTIVE output-layer LR, composing with both the
// OUTPUT_LR override and the OLR_RAMP/OLR_LATE schedule — it scales the
// early value (via LSTM_OUTPUT_LR below, which seeds olr_cur_) AND the late
// value (via LSTM_OLR_LATE), so the schedule's shape is preserved and the
// P1_MODE==2 replay stays exact (replay reads the same effective LR through
// LSTM_OUTPUT_LR_EFF either way). NOTE: with a flat OUTPUT_LR this knob is
// redundant (scale*0.06 == another OUTPUT_LR value); it exists so one
// multiplier can rescale a two-phase OLR schedule without retuning both
// phases. Flag-off, both macros collapse to the stock expansions —
// token-identical (no members added; passed via CFLAGS_DEFINES like the
// other LSTM_* knobs so all three flag groups agree).
#ifdef LSTM_LR_SCALE_OUT
#define LSTM_OUTPUT_LR (LSTM_OUTPUT_LR_BASE * (float)(LSTM_LR_SCALE_OUT))
#define LSTM_OLR_LATE ((float)(OLR_LATE) * (float)(LSTM_LR_SCALE_OUT))
#else
#define LSTM_OUTPUT_LR LSTM_OUTPUT_LR_BASE
#define LSTM_OLR_LATE ((float)(OLR_LATE))
#endif

// OLR_RAMP=<bytes> + OLR_LATE=<float> (code-dive 2026-07-11): position schedule on
// the output-layer LR — LSTM_OUTPUT_LR until olr_bytes_ >= OLR_RAMP, then OLR_LATE.
// The verified replay trap: P1_MODE==2 reconstructs past weights by replaying SGD
// updates at replay-time LR, so the switch may only happen at a cycle boundary AND
// after the replay block (the replay then uses the outgoing cycle's LR, exact).
// Flag-off, LSTM_OUTPUT_LR_EFF expands to the stock macro — token-identical.
#ifdef OLR_RAMP
#define LSTM_OUTPUT_LR_EFF olr_cur_
#else
#define LSTM_OUTPUT_LR_EFF LSTM_OUTPUT_LR
#endif

inline Lstm::Lstm(unsigned int input_size, unsigned int output_size, unsigned int
    num_cells, unsigned int num_layers, int horizon, float learning_rate,
    float gradient_clip) : input_history_(horizon),
    hidden_(num_cells * num_layers + 1), hidden_error_(num_cells),
    layer_input_(std::valarray<std::valarray<float>>(std::valarray<float>
    (input_size + 1 + num_cells * 2), num_layers), horizon),
#ifndef P1_MODE
    output_layer_(std::valarray<std::valarray<float>>(std::valarray<float>
   (num_cells * num_layers + 1), output_size), horizon),
#endif
    output_(std::valarray<float>(1.0 / output_size, output_size), horizon),
    learning_rate_(learning_rate), num_cells_(num_cells), epoch_(0),
    horizon_(horizon), input_size_(input_size), output_size_(output_size) {
  hidden_[hidden_.size() - 1] = 1;
#ifdef OLR_RAMP
  olr_cur_ = LSTM_OUTPUT_LR;
#endif
#ifdef P1_MODE
  output_w_.resize(output_size);
  for (unsigned int i = 0; i < output_size; ++i) {
    output_w_[i].resize(num_cells * num_layers + 1);
  }
#if P1_MODE == 2
  work_.resize(output_size);
  for (unsigned int i = 0; i < output_size; ++i) {
    work_[i].resize(num_cells * num_layers + 1);
  }
#endif
#endif
  for (int epoch = 0; epoch < horizon; ++epoch) {
    layer_input_[epoch][0].resize(1 + num_cells + input_size);
    for (unsigned int i = 0; i < num_layers; ++i) {
      layer_input_[epoch][i][layer_input_[epoch][i].size() - 1] = 1;
    }
  }
  for (unsigned int i = 0; i < num_layers; ++i) {
#ifdef LSTM_LR_SCALE_L2
    // Per-layer learning rate (track 3, 2026-07-12): the known depth
    // co-tuning gap — Gotty's deep configs failed without per-layer LR;
    // stock shares one Adam base rate across the stack. Layer 0 keeps the
    // base rate; every layer above it (the depth stack: index >= 1) gets
    // base * LSTM_LR_SCALE_L2. The scaled rate is baked in at construction
    // and flows into LstmLayer::learning_rate_ -> Adam() untouched anywhere
    // else (no replay interaction: P1_MODE replay reconstructs only the
    // output layer, whose LR is the separate LSTM_OUTPUT_LR stack above).
    const float layer_lr = (i == 0)
        ? learning_rate : learning_rate * (float)(LSTM_LR_SCALE_L2);
    layers_.emplace_back(layer_input_[0][i].size() + output_size, input_size_, output_size_,
        num_cells, horizon, gradient_clip, layer_lr);
#else
    layers_.emplace_back(layer_input_[0][i].size() + output_size, input_size_, output_size_,
        num_cells, horizon, gradient_clip, learning_rate);
#endif
  }
}

inline Lstm::~Lstm() {
  //SaveToDisk("lstm.dat");
}
/*
inline void Lstm::SaveToDisk(const std::string& path) {
  int last_epoch = epoch_ - 1;
  if (last_epoch == -1) last_epoch = horizon_ - 1;
  std::ofstream os(path, std::ios::binary | std::ios::out);
  if (!os.is_open()) return;
  for (int i = 0; i < output_size_; ++i) {
    os.write(reinterpret_cast<const char*>(&output_layer_[last_epoch][i][0]),
        std::streamsize(output_layer_[0][i].size() * sizeof(float)));
  }
  for (int i = 0; i < layers_.size(); ++i) {
    auto weights = layers_[i].Weights();
    for (int j = 0; j < weights.size(); ++j) {
      for (int k = 0; k < weights[j]->size(); ++k) {
        os.write(reinterpret_cast<const char*>(&(*weights[j])[k][0]),
          std::streamsize((*weights[j])[k].size() * sizeof(float)));
      }
    }
  }
  os.close();
}

inline void Lstm::LoadFromDisk(const std::string& path) {
  int last_epoch = epoch_ - 1;
  if (last_epoch == -1) last_epoch = horizon_ - 1;
  std::ifstream is(path, std::ios::binary | std::ios::in);
  if (!is.is_open()) return;
  for (int i = 0; i < output_size_; ++i) {
    is.read(reinterpret_cast<char*>(&output_layer_[last_epoch][i][0]),
        std::streamsize(output_layer_[0][i].size() * sizeof(float)));
  }
  for (int i = 0; i < layers_.size(); ++i) {
    auto weights = layers_[i].Weights();
    for (int j = 0; j < weights.size(); ++j) {
      for (int k = 0; k < weights[j]->size(); ++k) {
        is.read(reinterpret_cast<char*>(&(*weights[j])[k][0]),
          std::streamsize((*weights[j])[k].size() * sizeof(float)));
      }
    }
  }
  is.close();
}
*/
inline void Lstm::SetInput(const std::valarray<float>& input) {
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    std::copy(begin(input), begin(input) + input_size_,
        begin(layer_input_[epoch_][i]));
  }
}

inline std::valarray<float>& Lstm::Perceive(unsigned int input) {
  int last_epoch = epoch_ - 1;
  if (last_epoch == -1) last_epoch = horizon_ - 1;
  int old_input = input_history_[last_epoch];
  input_history_[last_epoch] = input;
  if (epoch_ == 0) {
#if defined(P1_MODE) && P1_MODE == 2
    for (unsigned int i = 0; i < output_size_; ++i) work_[i] = output_w_[i];
#endif
    for (int epoch = horizon_ - 1; epoch >= 0; --epoch) {
      for (int layer = layers_.size() - 1; layer >= 0; --layer) {
        int offset = layer * num_cells_;
#if defined(OERR_REGBLOCK)
#if !(defined(P1_MODE) && P1_MODE == 2)
#error "OERR_REGBLOCK implemented for the P1_MODE==2 (work_) seeding path only"
#endif
        // Register-blocked output-error seeding (byte-identical, MD5-gated):
        // block the hidden dim j OUTSIDE, loop the output symbols i INSIDE, so
        // hidden_error_[block] stays in registers across all output rows instead
        // of being loaded+stored once per symbol (~output_size-fold less
        // hidden_error_ traffic). Per-element accumulation order (i ascending)
        // and the (work_*error) product/FMA are unchanged => identical bytes.
        {
          const unsigned int H = (unsigned int)hidden_error_.size();
          const unsigned int OS = output_size_;
          float err[256];  // byte-level vocab: output_size_ <= 256
          for (unsigned int i = 0; i < OS; ++i)
            err[i] = (i == (unsigned int)input_history_[epoch])
                ? (output_[epoch][i] - 1.0f) : output_[epoch][i];
          constexpr unsigned int BLK = 32u;
          for (unsigned int jb = 0; jb < H; jb += BLK) {
            const unsigned int bn = (jb + BLK <= H) ? BLK : (H - jb);
            float acc[BLK];
            for (unsigned int k = 0; k < bn; ++k) acc[k] = hidden_error_[jb + k];
            for (unsigned int i = 0; i < OS; ++i) {
              const float e = err[i];
              const float* wr = &work_[i][jb + offset];
              for (unsigned int k = 0; k < bn; ++k) acc[k] += wr[k] * e;
            }
            for (unsigned int k = 0; k < bn; ++k) hidden_error_[jb + k] = acc[k];
          }
        }
#else
        for (unsigned int i = 0; i < output_size_; ++i) {
//          float error = 0;
//          if (i == input_history_[epoch]) error = output_[epoch][i] - 1;
//          else error = output_[epoch][i];
          float error = (i == input_history_[epoch]) ? (output_[epoch][i] - 1) : output_[epoch][i];
          for (unsigned int j = 0; j < hidden_error_.size(); ++j) {
#ifndef P1_MODE
            hidden_error_[j] += output_layer_[epoch][i][j + offset] * error;
#elif P1_MODE == 1
            hidden_error_[j] += output_w_[i][j + offset] * error;
#else
            hidden_error_[j] += work_[i][j + offset] * error;
#endif
          }
        }
#endif
        int prev_epoch = epoch - 1;
        if (prev_epoch == -1) prev_epoch = horizon_ - 1;
        int input_symbol = input_history_[prev_epoch];
        if (epoch == 0) input_symbol = old_input;
        layers_[layer].BackwardPass(layer_input_[epoch][layer], epoch, layer,
            input_symbol, &hidden_error_);
      }
#if defined(P1_MODE) && P1_MODE == 2
      if (epoch > 0) {
        // step work_ back from W[epoch] to W[epoch-1]: undo the rank-1 SGD
        // update lr*err(epoch-1) (x) h(epoch-1); h(epoch-1) is the hidden state
        // Predict copied into layer_input_[epoch] (bias input is constant 1)
        int pe = epoch - 1;
        for (unsigned int i = 0; i < output_size_; ++i) {
          float err = output_[pe][i] -
              ((i == input_history_[pe]) ? 1.0f : 0.0f);
          float s = LSTM_OUTPUT_LR_EFF * err;
          for (unsigned int l = 0; l < layers_.size(); ++l) {
            for (unsigned int c = 0; c < num_cells_; ++c) {
              work_[i][l * num_cells_ + c] +=
                  s * layer_input_[epoch][l][input_size_ + c];
            }
          }
          // bias column intentionally not reconstructed: the hidden_error_
          // seeding only reads columns [0, num_cells_*num_layers)
        }
      }
#endif
    }
#ifdef OLR_RAMP
    // cycle boundary, replay done: safe (and exact) point for the LR switch
    olr_bytes_ += horizon_;
    if (olr_bytes_ >= (unsigned long long)(OLR_RAMP)) olr_cur_ = LSTM_OLR_LATE;
#endif
  }

  for (unsigned int i = 0; i < output_size_; ++i) {
//    float error = 0;
//    if (i == input) error = output_[last_epoch][i] - 1;
//    else error = output_[last_epoch][i];
    float error = (i == input) ? (output_[last_epoch][i] - 1) : output_[last_epoch][i];
#ifndef P1_MODE
    output_layer_[epoch_][i] = output_layer_[last_epoch][i];
    output_layer_[epoch_][i] -= LSTM_OUTPUT_LR_EFF * error * hidden_;
#else
    output_w_[i] -= LSTM_OUTPUT_LR_EFF * error * hidden_;
#endif
  }
  return Predict(input);
}

inline std::valarray<float>& Lstm::Predict(unsigned int input) {
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    auto start = begin(hidden_) + i * num_cells_;
    std::copy(start, start + num_cells_, begin(layer_input_[epoch_][i]) +
        input_size_);
    layers_[i].ForwardPass(layer_input_[epoch_][i], input, &hidden_, i *
        num_cells_);
    if (i < layers_.size() - 1) {
      auto start2 = begin(layer_input_[epoch_][i + 1]) + num_cells_ +
          input_size_;
      std::copy(start, start + num_cells_, start2);
    }
  }
  for (unsigned int i = 0; i < output_size_; ++i) {
    float sum = 0;
    for (unsigned int j = 0; j < hidden_.size(); ++j) {
#ifndef P1_MODE
      sum += hidden_[j] * output_layer_[epoch_][i][j];
#else
      sum += hidden_[j] * output_w_[i][j];
#endif
    }
    output_[epoch_][i] = LSTM_EXPF(sum);
  }
  output_[epoch_] /= output_[epoch_].sum();
  int epoch = epoch_;
  ++epoch_;
  if (epoch_ == horizon_) epoch_ = 0;
  last_input_ = input;
  return output_[epoch];
}
