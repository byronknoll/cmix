#include "direct.h"

Direct::Direct(const unsigned long long& byte_context,
    const unsigned int& bit_context, int limit, float delta, int size) :
    byte_context_(byte_context), bit_context_(bit_context), index_(0),
    limit_(limit), delta_(delta), divisor_(1.0 / (limit + delta)),
    predictions_(size, std::array<float, 256>()),
    counts_(size, std::array<unsigned char, 256>()) {
  for (int i = 0; i < size; ++i) {
    predictions_[i].fill(0.5);
  }
}

float Direct::Predict() {
  return predictions_[index_][bit_context_];
}

void Direct::Perceive(int bit) {
  float divisor = divisor_;
  if (counts_[index_][bit_context_] < limit_) {
    ++counts_[index_][bit_context_];
    divisor = 1.0 / (counts_[index_][bit_context_] + delta_);
  }
  predictions_[index_][bit_context_] +=
      (bit - predictions_[index_][bit_context_]) * divisor;
}

void Direct::ByteUpdate() {
  index_ = byte_context_ % predictions_.size();
}
