#include "direct-hash.h"

DirectHash::DirectHash(const unsigned long long& byte_context,
    const unsigned int& bit_context, int limit, float delta, int size) :
    byte_context_(byte_context), bit_context_(bit_context), index_(0),
    limit_(limit), delta_(delta), divisor_(1.0 / (limit + delta)),
    predictions_(size, std::array<float, 256>()),
    counts_(size, std::array<unsigned char, 256>()),
    checksums_(size, 0) {
  for (int i = 0; i < size; ++i) {
    predictions_[i].fill(0.5);
    counts_[i].fill(0);
  }
}

const std::valarray<float>& DirectHash::Predict() {
  outputs_[0] = predictions_[index_][bit_context_];
  return outputs_;
}

void DirectHash::Perceive(int bit) {
  float divisor = divisor_;
  if (counts_[index_][bit_context_] < limit_) {
    ++counts_[index_][bit_context_];
    divisor = 1.0 / (counts_[index_][bit_context_] + delta_);
  }
  predictions_[index_][bit_context_] +=
      (bit - predictions_[index_][bit_context_]) * divisor;
}

void DirectHash::ByteUpdate() {
  index_ = byte_context_ % predictions_.size();
  for (int i = 0; i < 20; ++i) {
    if (checksums_[index_] == 0) {
      checksums_[index_] = byte_context_;
      break;
    }
    if (checksums_[index_] == byte_context_) break;
    if (i == 19) {
      predictions_[index_].fill(0.5);
      counts_[index_].fill(0);
      checksums_[index_] = byte_context_;
      break;
    }
    ++index_;
    if (index_ == predictions_.size()) index_ = 0;
  }
}
