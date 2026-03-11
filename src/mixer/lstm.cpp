#include "lstm.h"

#include <numeric>
#include <stdlib.h>
#include <fstream>
#include <iostream>

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

// dst[i] += src[i] * scalar
inline void AvxScaledAdd(float* dst, const float* src, float scalar, int n) {
  __m256 vs = _mm256_set1_ps(scalar);
  int i = 0;
  for (; i + 7 < n; i += 8) {
    __m256 vd = _mm256_loadu_ps(dst + i);
    __m256 vsrc = _mm256_loadu_ps(src + i);
    vd = _mm256_fmadd_ps(vsrc, vs, vd);
    _mm256_storeu_ps(dst + i, vd);
  }
  for (; i < n; ++i) {
    dst[i] += src[i] * scalar;
  }
}

// dst = src; dst -= scalar * mul
inline void AvxCopyAndScaledSub(float* dst, const float* src, const float* mul,
    float scalar, int n) {
  __m256 vs = _mm256_set1_ps(scalar);
  int i = 0;
  for (; i + 7 < n; i += 8) {
    __m256 vsrc = _mm256_loadu_ps(src + i);
    __m256 vmul = _mm256_loadu_ps(mul + i);
    vsrc = _mm256_fnmadd_ps(vs, vmul, vsrc);
    _mm256_storeu_ps(dst + i, vsrc);
  }
  for (; i < n; ++i) {
    dst[i] = src[i] - scalar * mul[i];
  }
}

}  // namespace
#endif

Lstm::Lstm(unsigned int input_size, unsigned int output_size, unsigned int
    num_cells, unsigned int num_layers, int horizon, float learning_rate,
    float gradient_clip) : input_history_(horizon),
    hidden_(num_cells * num_layers + 1), hidden_error_(num_cells),
    layer_input_(std::valarray<std::valarray<float>>(std::valarray<float>
    (input_size + 1 + num_cells * 2), num_layers), horizon),
    output_layer_(std::valarray<std::valarray<float>>(std::valarray<float>
    (num_cells * num_layers + 1), output_size), horizon),
    output_(std::valarray<float>(1.0 / output_size, output_size), horizon),
    learning_rate_(learning_rate), num_cells_(num_cells), epoch_(0),
    horizon_(horizon), input_size_(input_size), output_size_(output_size) {
  hidden_[hidden_.size() - 1] = 1;
  for (int epoch = 0; epoch < horizon; ++epoch) {
    layer_input_[epoch][0].resize(1 + num_cells + input_size);
    for (unsigned int i = 0; i < num_layers; ++i) {
      layer_input_[epoch][i][layer_input_[epoch][i].size() - 1] = 1;
    }
  }
  for (unsigned int i = 0; i < num_layers; ++i) {
    layers_.push_back(std::unique_ptr<LstmLayer>(new LstmLayer(
        layer_input_[0][i].size() + output_size, input_size_, output_size_,
        num_cells, horizon, gradient_clip, learning_rate)));
  }
  //LoadFromDisk("lstm.dat");
}

Lstm::~Lstm() {
  //SaveToDisk("lstm.dat");
}

void Lstm::SaveToDisk(const std::string& path) {
  int last_epoch = epoch_ - 1;
  if (last_epoch == -1) last_epoch = horizon_ - 1;
  std::ofstream os(path, std::ios::binary | std::ios::out);
  if (!os.is_open()) return;
  for (int i = 0; i < output_size_; ++i) {
    os.write(reinterpret_cast<const char*>(&output_layer_[last_epoch][i][0]),
        std::streamsize(output_layer_[0][i].size() * sizeof(float)));
  }
  for (int i = 0; i < layers_.size(); ++i) {
    auto weights = layers_[i]->Weights();
    for (int j = 0; j < weights.size(); ++j) {
      for (int k = 0; k < weights[j]->size(); ++k) {
        os.write(reinterpret_cast<const char*>(&(*weights[j])[k][0]),
          std::streamsize((*weights[j])[k].size() * sizeof(float)));
      }
    }
  }
  os.close();
}

void Lstm::LoadFromDisk(const std::string& path) {
  int last_epoch = epoch_ - 1;
  if (last_epoch == -1) last_epoch = horizon_ - 1;
  std::ifstream is(path, std::ios::binary | std::ios::in);
  if (!is.is_open()) return;
  for (int i = 0; i < output_size_; ++i) {
    is.read(reinterpret_cast<char*>(&output_layer_[last_epoch][i][0]),
        std::streamsize(output_layer_[0][i].size() * sizeof(float)));
  }
  for (int i = 0; i < layers_.size(); ++i) {
    auto weights = layers_[i]->Weights();
    for (int j = 0; j < weights.size(); ++j) {
      for (int k = 0; k < weights[j]->size(); ++k) {
        is.read(reinterpret_cast<char*>(&(*weights[j])[k][0]),
          std::streamsize((*weights[j])[k].size() * sizeof(float)));
      }
    }
  }
  is.close();
}

void Lstm::SetInput(const std::valarray<float>& input) {
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    std::copy(begin(input), begin(input) + input_size_,
        begin(layer_input_[epoch_][i]));
  }
}

std::valarray<float>& Lstm::Perceive(unsigned int input) {
  int last_epoch = epoch_ - 1;
  if (last_epoch == -1) last_epoch = horizon_ - 1;
  int old_input = input_history_[last_epoch];
  input_history_[last_epoch] = input;
  if (epoch_ == 0) {
    for (int epoch = horizon_ - 1; epoch >= 0; --epoch) {
      for (int layer = layers_.size() - 1; layer >= 0; --layer) {
        int offset = layer * num_cells_;
#ifdef __AVX2__
        int hn = hidden_error_.size();
        for (unsigned int i = 0; i < output_size_; ++i) {
          float error = (i == input_history_[epoch]) ?
              (output_[epoch][i] - 1) : output_[epoch][i];
          AvxScaledAdd(&hidden_error_[0],
              &output_layer_[epoch][i][offset], error, hn);
        }
#else
        for (unsigned int i = 0; i < output_size_; ++i) {
          float error = (i == input_history_[epoch]) ? (output_[epoch][i] - 1) : output_[epoch][i];
          for (unsigned int j = 0; j < hidden_error_.size(); ++j) {
            hidden_error_[j] += output_layer_[epoch][i][j + offset] * error;
          }
        }
#endif
        int prev_epoch = epoch - 1;
        if (prev_epoch == -1) prev_epoch = horizon_ - 1;
        int input_symbol = input_history_[prev_epoch];
        if (epoch == 0) input_symbol = old_input;
        layers_[layer]->BackwardPass(layer_input_[epoch][layer], epoch, layer,
            input_symbol, &hidden_error_);
      }
    }
  }

#ifdef __AVX2__
  int hn = hidden_.size();
  for (unsigned int i = 0; i < output_size_; ++i) {
    float error = (i == input) ? (output_[last_epoch][i] - 1) : output_[last_epoch][i];
    float lr_error = learning_rate_ * error;
    AvxCopyAndScaledSub(&output_layer_[epoch_][i][0],
        &output_layer_[last_epoch][i][0], &hidden_[0], lr_error, hn);
  }
#else
  for (unsigned int i = 0; i < output_size_; ++i) {
    float error = (i == input) ? (output_[last_epoch][i] - 1) : output_[last_epoch][i];
    output_layer_[epoch_][i] = output_layer_[last_epoch][i];
    output_layer_[epoch_][i] -= learning_rate_ * error * hidden_;
  }
#endif
  return Predict(input);
}

std::valarray<float>& Lstm::Predict(unsigned int input) {
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    auto start = begin(hidden_) + i * num_cells_;
    std::copy(start, start + num_cells_, begin(layer_input_[epoch_][i]) +
        input_size_);
    layers_[i]->ForwardPass(layer_input_[epoch_][i], input, &hidden_, i *
        num_cells_);
    if (i < layers_.size() - 1) {
      auto start2 = begin(layer_input_[epoch_][i + 1]) + num_cells_ +
          input_size_;
      std::copy(start, start + num_cells_, start2);
    }
  }

#ifdef __AVX2__
  int hn = hidden_.size();
  const float* hp = &hidden_[0];
  float max_out = 0;
  for (unsigned int i = 0; i < output_size_; ++i) {
    float sum = AvxDotProduct(hp, &output_layer_[epoch_][i][0], hn);
    output_[epoch_][i] = sum;
    max_out = std::max(sum, max_out);
  }
#else
  float max_out = 0;
  for (unsigned int i = 0; i < output_size_; ++i) {
    float sum = 0;
    for (unsigned int j = 0; j < hidden_.size(); ++j) {
      sum += hidden_[j] * output_layer_[epoch_][i][j];
    }
    output_[epoch_][i] = sum;
    max_out = std::max(sum, max_out);
  }
#endif

  for (unsigned int i = 0; i < output_size_; ++i) {
    output_[epoch_][i] = exp(output_[epoch_][i] - max_out);
  }
  output_[epoch_] /= output_[epoch_].sum();
  int epoch = epoch_;
  ++epoch_;
  if (epoch_ == horizon_) epoch_ = 0;
  return output_[epoch];
}
