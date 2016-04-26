#include "sse.h"

SSE::SSE(const Logistic& logistic, const unsigned long long& byte_context,
    const unsigned int& bit_context, unsigned int num_buckets, float delta,
    unsigned long long table_size) : bit_context_(bit_context),
    byte_context_(byte_context), num_buckets_(num_buckets),
    bucket_(num_buckets / 2), divisor_(1 / delta),
    predictions_(table_size, std::vector<std::array<float, 256>>(
        num_buckets, std::array<float, 256>())) {
  for (unsigned int i = 0; i < table_size; ++i) {
    for (unsigned int j = 0; j < num_buckets; ++j) {
      float p = (j + 0.5) / num_buckets;
      predictions_[i][j].fill(logistic.Squash(20 * p - 10));
    }
  }
}

float SSE::Process(float input) {
  bucket_ = num_buckets_ * ((input + 10) / 20);
  if (bucket_ >= num_buckets_) bucket_ = num_buckets_ - 1;
  else if (bucket_ < 0) bucket_ = 0;
  return predictions_[byte_context_][bucket_][bit_context_];
}

void SSE::Perceive(int bit) {
  predictions_[byte_context_][bucket_][bit_context_] +=
      (bit - predictions_[byte_context_][bucket_][bit_context_]) * divisor_;
}

