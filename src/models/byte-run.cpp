#include "byte-run.h"

ByteRun::ByteRun(const unsigned long long& byte_context,
    const unsigned int& bit_context, float delta,
    unsigned long long map_size) :  byte_context_(byte_context),
    bit_context_(bit_context), byte_prediction_(0), run_length_(0),
    map_index_(0), bit_pos_(128), divisor_(1.0 / delta), map_(map_size, 0),
    counts_(map_size, 0) {
  for (int i = 0; i < 256; ++i) {
    predictions_[i] = 0.5 + ((i + 0.5) / 512);
  }
}

float ByteRun::Predict() {
  if (byte_prediction_ & bit_pos_) return predictions_[run_length_];
  return 1 - predictions_[run_length_];
}

void ByteRun::Perceive(int bit) {
  int match = 0;
  if (((byte_prediction_ & bit_pos_) != 0) == bit) match = 1;

  predictions_[run_length_] += (match - predictions_[run_length_]) * divisor_;

  bit_pos_ /= 2;
  if (match == 0) run_length_ = 0;

  if (bit_context_ >= 128) {
    unsigned char byte = (bit_context_ << 1) + bit;
    if (byte == byte_prediction_) {
      if (counts_[map_index_] < 255) ++counts_[map_index_];
    } else {
      map_[map_index_] = byte;
      counts_[map_index_] = 0;
    }
  }
}

void ByteRun::ByteUpdate() {
  map_index_ = byte_context_ % map_.size();
  byte_prediction_ = map_[map_index_];
  run_length_ = counts_[map_index_];
  bit_pos_ = 128;
}
