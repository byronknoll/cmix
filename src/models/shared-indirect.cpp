#include "shared-indirect.h"
#include <stdlib.h>

SharedIndirect::SharedIndirect(const State& state,
    const unsigned long long& byte_context,
    const unsigned int& bit_context, float delta,
    std::vector<unsigned char>& map) :  byte_context_(byte_context),
    bit_context_(bit_context), map_index_(0), map_offset_(0),
    divisor_(1.0 / delta), state_(state), map_(map) {
  map_offset_ = rand() % (map_.size() - 257);
  for (int i = 0; i < 256; ++i) {
    predictions_[i] = state_.InitProbability(i);
  }
}

float SharedIndirect::Predict() {
  map_index_ += bit_context_;
  return predictions_[map_[map_index_]];
}

void SharedIndirect::Perceive(int bit) {
  int state = map_[map_index_];
  predictions_[state] += (bit - predictions_[state]) * divisor_;
  map_[map_index_] = state_.Next(state, bit);
  map_index_ -= bit_context_;
}

void SharedIndirect::ByteUpdate() {
  map_index_ = (257 * byte_context_ + map_offset_) % (map_.size() - 257);
}
