#include "bracket.h"

Bracket::Bracket(const unsigned int& bit_context) : active_(-1), distance_(0),
    max_(200), byte_(bit_context),
    stats_(256, std::vector<std::pair<double, double>>(max_, {1, 256})) {
  brackets_ = {{'(',')'}, {'{','}'}, {'[',']'}, {'<','>'}};
}

void Bracket::ByteUpdate() {
  probs_ = 1./256;
  if (active_ == -1 || brackets_.find(byte_) != brackets_.end()) {
    if (brackets_.find(byte_) != brackets_.end()) {
      active_ = byte_;
      distance_ = 0;
      float p = stats_[active_][distance_].first /
          stats_[active_][distance_].second;
      probs_ = (1 - p) / 255;
      probs_[brackets_[active_]] = p;
    }
  } else {
    ++stats_[active_][distance_].second;
    if (brackets_[active_] == byte_) {
      ++stats_[active_][distance_].first;
    }
    if (brackets_[active_] == byte_ || distance_ >= max_ - 1) {
      active_ = -1;
    } else {
      ++distance_;
      float p = stats_[active_][distance_].first /
          stats_[active_][distance_].second;
      probs_ = (1 - p) / 255;
      probs_[brackets_[active_]] = p;
    }
  }
  ByteModel::ByteUpdate();
}
