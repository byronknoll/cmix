#include "direct.h"

Direct::Direct(const unsigned long long& byte_context,
    const unsigned int& bit_context, int limit, float delta, int size) :
    byte_context_(byte_context), bit_context_(bit_context), limit_(limit),
    delta_(delta), divisor_(1.0 / (limit + delta)),
    predictions_(size, std::array<float, 256>()),
    counts_(size, std::array<unsigned char, 256>()) {
  for (int i = 0; i < size; ++i) {
    predictions_[i].fill(0.5);
    counts_[i].fill(0);
  }
}

const std::valarray<float>& Direct::Predict() {
  outputs_[0] = predictions_[byte_context_][bit_context_];
  return outputs_;
}

void Direct::Perceive(int bit) {
  float divisor = divisor_;
  if (counts_[byte_context_][bit_context_] < limit_) {
    ++counts_[byte_context_][bit_context_];
    divisor = 1.0 / (counts_[byte_context_][bit_context_] + delta_);
  }
  predictions_[byte_context_][bit_context_] +=
      (bit - predictions_[byte_context_][bit_context_]) * divisor;
}
