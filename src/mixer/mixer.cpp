#include "mixer.h"

#include "sigmoid.h"

#include <numeric>
#include <math.h>

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

// weights -= update * inputs
inline void AvxWeightUpdate(float* weights, const float* inputs, float update,
    int n) {
  __m256 vupdate = _mm256_set1_ps(update);
  int i = 0;
  for (; i + 7 < n; i += 8) {
    __m256 vw = _mm256_loadu_ps(weights + i);
    __m256 vi = _mm256_loadu_ps(inputs + i);
    vw = _mm256_fnmadd_ps(vupdate, vi, vw);
    _mm256_storeu_ps(weights + i, vw);
  }
  for (; i < n; ++i) {
    weights[i] -= update * inputs[i];
  }
}

// weights *= scale
inline void AvxScale(float* weights, float scale, int n) {
  __m256 vs = _mm256_set1_ps(scale);
  int i = 0;
  for (; i + 7 < n; i += 8) {
    __m256 vw = _mm256_loadu_ps(weights + i);
    vw = _mm256_mul_ps(vw, vs);
    _mm256_storeu_ps(weights + i, vw);
  }
  for (; i < n; ++i) {
    weights[i] *= scale;
  }
}

}  // namespace
#endif

Mixer::Mixer(const std::valarray<float>& inputs,
    const std::vector<float>& extra_inputs,
    const unsigned long long& context, float learning_rate,
    unsigned int extra_input_size) : inputs_(inputs),
    extra_inputs_vec_(extra_inputs), extra_inputs_(extra_input_size), p_(0.5),
    learning_rate_(learning_rate), context_(context), max_steps_(1), steps_(0)
    {}

ContextData* Mixer::GetContextData() {
  ContextData* data;
  unsigned long long limit = 10000;
  if (context_map_.size() >= limit && context_map_.find(context_) == context_map_.end()) {
    data = context_map_[0xDEADBEEF].get();
    if (data == nullptr) {
      context_map_[0xDEADBEEF] = std::unique_ptr<ContextData>(
          new ContextData(inputs_.size(), extra_inputs_.size()));
      data = context_map_[0xDEADBEEF].get();
    }
  } else {
    data = context_map_[context_].get();
    if (data == nullptr) {
      context_map_[context_] = std::unique_ptr<ContextData>(
          new ContextData(inputs_.size(), extra_inputs_.size()));
      data = context_map_[context_].get();
    }
  }

  return data;
}

float Mixer::Mix() {
  ContextData* data = GetContextData();

#ifdef __AVX2__
  int n = inputs_.size();
  p_ = AvxDotProduct(&inputs_[0], &data->weights[0], n);
#else
  float p = 0;
  for (int i = 0; i < inputs_.size(); ++i) {
    p += inputs_[i] * data->weights[i];
  }
  p_ = p;
#endif

  for (unsigned int i = 0; i < extra_inputs_.size(); ++i) {
    extra_inputs_[i] = extra_inputs_vec_[i];
  }

#ifdef __AVX2__
  int ne = extra_inputs_.size();
  if (ne > 0) {
    p_ += AvxDotProduct(&extra_inputs_[0], &data->extra_weights[0], ne);
  }
#else
  float e = 0;
  for (unsigned int i = 0; i < extra_inputs_.size(); ++i) {
    e += extra_inputs_[i] * data->extra_weights[i];
  }
  p_ += e;
#endif

  return p_;
}

void Mixer::Perceive(int bit) {
  ContextData* data = GetContextData();
  float decay = 0.9 / pow(0.0000001 * steps_ + 0.8, 0.8);
  decay *= 1.5 - ((1.0 * data->steps) / max_steps_);
  float update = decay * learning_rate_ * (Sigmoid::Logistic(p_) - bit);
  ++steps_;
  ++data->steps;
  if (data->steps > max_steps_) {
    max_steps_ = data->steps;
  }

#ifdef __AVX2__
  int n = inputs_.size();
  AvxWeightUpdate(&data->weights[0], &inputs_[0], update, n);

  int ne = extra_inputs_.size();
  if (ne > 0) {
    AvxWeightUpdate(&data->extra_weights[0], &extra_inputs_[0], update, ne);
  }

  if ((data->steps & 1023) == 0) {
    AvxScale(&data->weights[0], 1.0f - 3.0e-6f, n);
    if (ne > 0) {
      AvxScale(&data->extra_weights[0], 1.0f - 3.0e-6f, ne);
    }
  }
#else
  data->weights -= update * inputs_;
  data->extra_weights -= update * extra_inputs_;
  if ((data->steps & 1023) == 0) {
    data->weights *= 1.0f - 3.0e-6f;
    data->extra_weights *= 1.0f - 3.0e-6f;
  }
#endif
}
