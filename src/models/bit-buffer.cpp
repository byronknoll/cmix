#include "bit-buffer.h"

BitBuffer::BitBuffer(int size, float delta, int limit) : bit_buffer_(size, 0),
    bit_pos_(0), delta_(delta), limit_(limit), divisor_(1.0 / (delta + limit)) {
  predictions_.fill(0.5);
}

float BitBuffer::Predict() {
  return predictions_[bit_buffer_[bit_pos_]];
}

void BitBuffer::Perceive(int bit) {
  bool c = bit_buffer_[bit_pos_];
  bit_buffer_[bit_pos_] = bit;
  ++bit_pos_;
  if (bit_pos_ == bit_buffer_.size()) bit_pos_ = 0;

  float divisor = divisor_;
  if (counts_[c] < limit_) {
    ++counts_[c];
    divisor = 1.0 / (counts_[c] + delta_);
  }
  predictions_[c] += (bit - predictions_[c]) * divisor;
}

