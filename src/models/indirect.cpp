#include "indirect.h"

Indirect::Indirect(const State& state,
    const unsigned long long& byte_context,
    const unsigned int& bit_context, float delta,
    unsigned long long map_size) :  byte_context_(byte_context),
    bit_context_(bit_context), map_index_(0), divisor_(1.0 / delta),
    state_(state), map_(map_size, std::array<unsigned char, 256>()) {
  for (int i = 0; i < 256; ++i) {
    predictions_[i] = state_.InitProbability(i);
  }
  for (unsigned int i = 0; i < map_.size(); ++i) {
    map_[i][0] = 255;
  }
}

float Indirect::Predict() {
  return predictions_[map_[map_index_][bit_context_]];
}

void Indirect::Perceive(int bit) {
  int state = map_[map_index_][bit_context_];
  predictions_[state] += (bit - predictions_[state]) * divisor_;
  map_[map_index_][bit_context_] = state_.Next(state, bit);
}

void Indirect::ByteUpdate() {
  map_index_ = byte_context_ % map_.size();
  unsigned char checksum = byte_context_ / map_.size();
  if (map_[map_index_][0] != 255 && map_[map_index_][0] != checksum) {
    map_index_ ^= 1;
    if (map_index_ == map_.size()) map_index_ = 0;
    if (map_[map_index_][0] != 255 && map_[map_index_][0] != checksum) {
      map_index_ ^= 2;
      if (map_index_ >= map_.size()) map_index_ = 0;
      if (map_[map_index_][0] != 255 && map_[map_index_][0] != checksum) {
        for (int i = 1; i < 256; ++i) map_[map_index_][i] = 0;
      }
    }
  }
  map_[map_index_][0] = checksum;
}
