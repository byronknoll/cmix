#include "indirect.h"

Indirect::Indirect(const State& state,
    const unsigned long long& byte_context,
    const unsigned int& bit_context, float delta,
    unsigned long long map_size) :  byte_context_(byte_context),
    bit_context_(bit_context), map_index_(0), divisor_(1.0 / delta),
    state_(state), map_(256 * map_size, 0) {
  for (int i = 0; i < 256; ++i) {
    predictions_[i] = state_.InitProbability(i);
  }
}

float Indirect::Predict() {
  map_index_ += bit_context_;
  return predictions_[map_[map_index_]];
}

void Indirect::Perceive(int bit) {
  int state = map_[map_index_];
  predictions_[state] += (bit - predictions_[state]) * divisor_;
  map_[map_index_] = state_.Next(state, bit);
  map_index_ -= bit_context_;
}

void Indirect::ByteUpdate() {
  map_index_ = (257 * byte_context_) % (map_.size() - 257);
}
