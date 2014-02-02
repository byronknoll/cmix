#include "match.h"

Match::Match(const std::vector<unsigned char>& history,
    const unsigned long long& byte_context, const unsigned int& bit_context,
    int limit, float delta, unsigned long long map_size) : history_(history),
    byte_context_(byte_context), bit_context_(bit_context), 
    history_pos_(0), cur_match_(0), cur_byte_(0),
    match_length_(0), limit_(limit), delta_(delta),
    divisor_(1.0 / (limit + delta)), map_(256 * map_size, 0) {
  for (int i = 0; i < 256; ++i) {
    predictions_[i] = 0.5 + (i + 0.5) / 512;
    unsigned char orig = i;
    unsigned char reverse = 0;
    for (int j = 0; j < 8; ++j) {
      reverse <<= 1;
      reverse += orig & 1;
      orig >>= 1;
    }
    reverse_table_[i] = reverse;
  }
}

float Match::Predict() {
  if (cur_byte_ & 1) return predictions_[match_length_];
  return 1 - predictions_[match_length_];
}

void Match::Perceive(int bit) {
  int match = 0;
  if (bit == (cur_byte_ & 1)) match = 1;

  float divisor = divisor_;
  if (counts_[match_length_] < limit_) {
    ++counts_[match_length_];
    divisor = 1.0 / (counts_[match_length_] + delta_);
  }
  predictions_[match_length_] +=
      (match - predictions_[match_length_]) * divisor;

  cur_byte_ >>= 1;
  if (match) {
    if (match_length_ < 255) ++match_length_;
  } else {
    match_length_ = 0;
  }

  if (bit_context_ >= 128) {
    map_[byte_context_ % map_.size()] = history_pos_;
    ++history_pos_;
  }
}

void Match::ByteUpdate() {
  if (match_length_ < 8) {
    cur_match_ = map_[byte_context_ % map_.size()];
  } else {
    ++cur_match_;
  }
  cur_byte_ = reverse_table_[history_.at(cur_match_)];
}
