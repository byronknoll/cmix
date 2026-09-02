#ifndef LSTM_QUANT_H
#define LSTM_QUANT_H

// -DLSTM_QUANT={8,16}: quantized LSTM gate-matvec forward pass.
//
// Scope (correctness-first design):
//   * ONLY the three gate matvecs in LstmLayer::ForwardPass(NeuronLayer&,..)
//     are quantized (the 160x526-shape hot loop). The one-hot embedding
//     column W[i][input_symbol] and the constant-1 bias column (last matvec
//     column; layer_input_[e][l][size-1] == 1 forever, set in the Lstm ctor
//     and never overwritten) stay float-exact.
//   * Weight STORAGE stays float; BPTT/Adam (BackwardPass) stays float.
//     The quantized weights are a derived cache, rebuilt lazily after Adam
//     runs (epoch==0 of the BPTT cycle => once per horizon steps; amortized
//     requant cost ~= 2 passes over the weights per horizon steps ~= 1.6%
//     of the matvec work at horizon=128).
//   * Weights: symmetric per-row scale, q in [-QMAX, +QMAX] (never -QMAX-1,
//     so |product| <= QMAX^2 and int16 madd pair-sums fit int32 exactly).
//   * Activations: symmetric per-step scale over the whole input vector.
//   * Accumulation is pure integer => order-independent => BIT-EXACT across
//     CPUs for a fixed binary (the determinism bonus: no FMA-vs-mul-add or
//     vector-width reassociation divergence in the accumulation; kernel
//     selection is compile-time only, no runtime dispatch).
//
// Variants:
//   LSTM_QUANT=8  : int8 weights+acts, int32 accumulator.
//                   x86: AVX512-VNNI (vpdpbusd, activations biased +128 into
//                   u8, corrected by 128*rowsum) when __AVX512VNNI__;
//                   else AVX2 (sign-extend to i16 + vpmaddwd; exact, no
//                   saturation -- vpmaddubsw is deliberately NOT used, its
//                   i16 pair-sum saturates at 2*255*127 > 32767).
//   LSTM_QUANT=16 : int16 weights+acts, int64 accumulator (int32 would
//                   overflow: n*32767^2 > 2^31 for n>2). x86: AVX2 vpmaddwd
//                   (pair-sums <= 2*32767^2 = 2147352578 < 2^31, exact) then
//                   widened to i64 every step. AVX512 vpdpwssd is NOT used:
//                   it accumulates in i32 and would wrap.
//   Every SIMD kernel computes the exact same integer as LstmQuantDotScalar
//   (integer addition is associative; no intermediate overflows by the
//   bounds above, audited in p13_quant/lstm_quant_selftest.cpp).
//
// Dequantization (the only float ops added, all IEEE mul/add/cvt):
//   f_i = W[i][input_symbol] + W[i][bias_col] + (float)acc_i * (sw_i * sx)
//
// -DLSTM_QUANT_VERIFY adds a runtime scalar crosscheck of the SIMD kernels
// (abort on mismatch; NOT assert -- NDEBUG kills asserts in this codebase).
// The x86 intrinsic paths compile under -target x86_64 locally but execute
// only on the VM: VM-verify.

#ifdef LSTM_QUANT

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#if LSTM_QUANT == 8
typedef int8_t  LstmQuantWeight;
typedef int8_t  LstmQuantAct;
typedef int32_t LstmQuantAcc;
#define LSTM_QUANT_QMAX 127
#elif LSTM_QUANT == 16
typedef int16_t LstmQuantWeight;
typedef int16_t LstmQuantAct;
typedef int64_t LstmQuantAcc;
#define LSTM_QUANT_QMAX 32767
#else
#error "LSTM_QUANT must be 8 (int8: VNNI, AVX2 fallback) or 16 (int16: AVX2)"
#endif

// Row stride padding (elements). 64 covers one AVX-512 int8 vector; padded
// weight elements are 0 so they contribute exactly 0 to the integer dot.
#define LSTM_QUANT_PAD 64u

// Compile-time kernel selection. No runtime dispatch: one binary always runs
// one kernel, so results are bit-exact across the CPUs that binary runs on.
#if LSTM_QUANT == 8 && defined(__AVX512VNNI__) && defined(__AVX512F__)
#define LSTM_QUANT_KERNEL_VNNI 1
#include <immintrin.h>
#elif defined(__AVX2__)
#define LSTM_QUANT_KERNEL_AVX2 1
#include <immintrin.h>
#else
#define LSTM_QUANT_KERNEL_SCALAR 1
#endif

// Round-to-nearest-even under the default FP environment on every target
// (x86 cvtps2si via MXCSR default; arm64 fcvtns). Matches what
// _mm256_cvtps_epi32 would do, so a future vectorized quantizer stays
// bit-compatible.
static inline int32_t LstmQuantRound(float v) {
  return (int32_t)lrintf(v);
}

// Quantize `count` floats from src into dst (padded to `stride` with zeros).
// Symmetric: q = clamp(round(x * QMAX/maxabs), -QMAX, QMAX). Returns the
// dequantization scale s = maxabs/QMAX (x ~= q*s); rowsum_out gets sum(q)
// (the VNNI +128-bias correction term; also cheap to keep in the =16 build).
// maxabs==0 => all-zero q, scale 1.0 (acc will be 0; avoids inf/NaN scales).
static inline float LstmQuantizeRow(const float* src, unsigned int count,
    LstmQuantWeight* dst, unsigned int stride, int32_t* rowsum_out) {
  float maxabs = 0.0f;
  for (unsigned int j = 0; j < count; ++j) {
    const float a = fabsf(src[j]);
    if (a > maxabs) maxabs = a;
  }
  if (maxabs == 0.0f) {
    memset(dst, 0, (size_t)stride * sizeof(LstmQuantWeight));
    *rowsum_out = 0;
    return 1.0f;
  }
  const float scale = maxabs * (1.0f / (float)LSTM_QUANT_QMAX);
  const float inv = (float)LSTM_QUANT_QMAX / maxabs;
  int32_t rowsum = 0;
  for (unsigned int j = 0; j < count; ++j) {
    int32_t q = LstmQuantRound(src[j] * inv);
    if (q > LSTM_QUANT_QMAX) q = LSTM_QUANT_QMAX;
    if (q < -LSTM_QUANT_QMAX) q = -LSTM_QUANT_QMAX;
    dst[j] = (LstmQuantWeight)q;
    rowsum += q;
  }
  for (unsigned int j = count; j < stride; ++j) dst[j] = 0;
  *rowsum_out = rowsum;
  return scale;
}

// Activation flavor: same math, no rowsum.
static inline float LstmQuantizeActs(const float* src, unsigned int count,
    LstmQuantAct* dst, unsigned int stride) {
  int32_t unused;
  return LstmQuantizeRow(src, count, dst, stride, &unused);
}


// THE determinism spec: exact integer dot product. Every SIMD kernel below
// must return exactly this value (auto-vectorization of this loop is safe:
// integer addition is associative, and by the range bounds no intermediate
// can overflow LstmQuantAcc).
static inline LstmQuantAcc LstmQuantDotScalar(const LstmQuantWeight* w,
    const LstmQuantAct* x, unsigned int n) {
  LstmQuantAcc acc = 0;
  for (unsigned int j = 0; j < n; ++j) {
    acc += (LstmQuantAcc)w[j] * (LstmQuantAcc)x[j];
  }
  return acc;
}

#ifdef LSTM_QUANT_KERNEL_VNNI
// int8 x int8 via u8 x s8 vpdpbusd (VM-verify: Cascade Lake AVX512_VNNI).
// xb[j] = (uint8)(x[j] + 128); sum((x+128)*w) - 128*sum(w) == sum(x*w),
// exact in int32: n <= 65536 gives |sum((x+128)*w)| <= n*255*127 and the
// bound check lives in the caller's ctor. vpdpbusd's inner 4-product sum
// (<= 4*255*127 = 129540) never saturates s32; accumulator adds wrap only
// past 2^31 which the n-bound excludes.
static inline int32_t LstmQuantDotVnni(const int8_t* w, const uint8_t* xb,
    unsigned int n, int32_t rowsum) {
  __m512i acc = _mm512_setzero_si512();
  for (unsigned int j = 0; j < n; j += 64) {
    const __m512i a = _mm512_loadu_si512((const void*)(xb + j));
    const __m512i b = _mm512_loadu_si512((const void*)(w + j));
    acc = _mm512_dpbusd_epi32(acc, a, b);
  }
  return _mm512_reduce_add_epi32(acc) - 128 * rowsum;
}
#endif

#ifdef LSTM_QUANT_KERNEL_AVX2
#if LSTM_QUANT == 8
// int8 AVX2 fallback: sign-extend both sides to i16, vpmaddwd, i32 lanes.
// Exact: pair-sums <= 2*127^2 = 32258; per-lane totals <= (n/16)*32258,
// safe to n ~ 1e6 (we use n <= ~1600). (VM-verify)
static inline int32_t LstmQuantDotAvx2(const int8_t* w, const int8_t* x,
    unsigned int n) {
  __m256i acc = _mm256_setzero_si256();
  for (unsigned int j = 0; j < n; j += 16) {
    const __m256i a =
        _mm256_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(x + j)));
    const __m256i b =
        _mm256_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(w + j)));
    acc = _mm256_add_epi32(acc, _mm256_madd_epi16(a, b));
  }
  __m128i lo = _mm256_castsi256_si128(acc);
  const __m128i hi = _mm256_extracti128_si256(acc, 1);
  lo = _mm_add_epi32(lo, hi);
  lo = _mm_add_epi32(lo, _mm_shuffle_epi32(lo, 0x4E));
  lo = _mm_add_epi32(lo, _mm_shuffle_epi32(lo, 0xB1));
  return _mm_cvtsi128_si32(lo);
}
#else
// int16 AVX2: vpmaddwd pair-sums are exact in i32 (<= 2*32767^2 < 2^31-1,
// and -32768 never occurs since quantization clamps to +/-32767); widened
// to i64 every iteration because n pair-sums can exceed 2^31. (VM-verify)
static inline int64_t LstmQuantDotAvx2(const int16_t* w, const int16_t* x,
    unsigned int n) {
  __m256i acc = _mm256_setzero_si256();  // 4 x i64
  for (unsigned int j = 0; j < n; j += 16) {
    const __m256i a = _mm256_loadu_si256((const __m256i*)(x + j));
    const __m256i b = _mm256_loadu_si256((const __m256i*)(w + j));
    const __m256i p = _mm256_madd_epi16(a, b);  // 8 x i32, exact
    acc = _mm256_add_epi64(acc,
        _mm256_cvtepi32_epi64(_mm256_castsi256_si128(p)));
    acc = _mm256_add_epi64(acc,
        _mm256_cvtepi32_epi64(_mm256_extracti128_si256(p, 1)));
  }
  __m128i lo = _mm256_castsi256_si128(acc);
  const __m128i hi = _mm256_extracti128_si256(acc, 1);
  lo = _mm_add_epi64(lo, hi);
  lo = _mm_add_epi64(lo, _mm_unpackhi_epi64(lo, lo));
  return (int64_t)_mm_cvtsi128_si64(lo);
}
#endif
#endif  // LSTM_QUANT_KERNEL_AVX2

// Single entry point. n must be a multiple of LSTM_QUANT_PAD (padded with
// zero weights). xb is the +128-biased activation buffer (VNNI builds only;
// pass nullptr elsewhere). rowsum is sum of the quantized weight row.
static inline LstmQuantAcc LstmQuantDot(const LstmQuantWeight* w,
    const LstmQuantAct* x, const uint8_t* xb, unsigned int n,
    int32_t rowsum) {
#if defined(LSTM_QUANT_KERNEL_SCALAR)
  (void)xb; (void)rowsum;
  return LstmQuantDotScalar(w, x, n);
#else
#if defined(LSTM_QUANT_KERNEL_VNNI)
  const LstmQuantAcc acc = LstmQuantDotVnni(w, xb, n, rowsum);
#else
  (void)xb; (void)rowsum;
  const LstmQuantAcc acc = LstmQuantDotAvx2(w, x, n);
#endif
#ifdef LSTM_QUANT_VERIFY
  const LstmQuantAcc ref = LstmQuantDotScalar(w, x, n);
  if (acc != ref) {
    fprintf(stderr,
        "LSTM_QUANT_VERIFY: SIMD kernel mismatch: acc=%lld ref=%lld n=%u\n",
        (long long)acc, (long long)ref, n);
    abort();
  }
#endif
  return acc;
#endif
}

#endif  // LSTM_QUANT
#endif  // LSTM_QUANT_H
