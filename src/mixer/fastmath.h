#ifndef LSTM_FASTMATH_H
#define LSTM_FASTMATH_H

#include <valarray>
#include <cstdint>

// Float exp, degree-5 polynomial, single-constant range reduction:
// ~3.4e-6 max rel error for |x|<20 (the region gate logits occupy), up to
// ~7e-6 near the clamp boundary.
// Input clamped to [-87, 87]: exp(87) ~ 6.1e37 < FLT_MAX and n+127 <= 253,
// so the exponent reconstruction can never produce Inf; exp(-87) ~ 1.6e-38
// stays above FLT_MIN (no denormals on the low end either).
inline float fast_expf(float x) {
  x = x < -87.0f ? -87.0f : (x > 87.0f ? 87.0f : x);
  float z = x * 1.44269504f;
  int   n = (int)(z + (z >= 0.0f ? 0.5f : -0.5f));
  float r = x - (float)n * 0.693147182f;
  float p = 1.0f + r*(1.0f + r*(0.499999994f + r*(0.166668282f + r*(0.041656744f + r*0.008273227f))));
  union { float f; int32_t i; } u; u.i = (int32_t)((uint32_t)(n + 127) << 23);
  return p * u.f;
}

inline float fast_logistic(float x) {
  return 1.0f / (1.0f + fast_expf(-x));
}

// tanh via exp(2x); the clamp inside fast_expf saturates both tails correctly.
// Relative error blows up near 0 (cancellation) but absolute error stays
// ~1.5e-6, which is what matters for gate activations.
inline float fast_tanh(float x) {
  float e = fast_expf(2.0f * x);
  return (e - 1.0f) / (e + 1.0f);
}

inline std::valarray<float> fast_tanh_vec(const std::valarray<float>& v) {
  std::valarray<float> r(v.size());
  for (size_t i = 0; i < v.size(); ++i) r[i] = fast_tanh(v[i]);
  return r;
}

#endif
