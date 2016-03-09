#include "bracket.h"

Bracket::Bracket(const unsigned int& bit_context) : max_(200),
    byte_(bit_context), stats_(256, std::vector<std::pair<double, double>>
    (max_, {1, 256})) {
  brackets_ = {{'(',')'}, {'{','}'}, {'[',']'}, {'<','>'}, {'\'','\''},
      {'"','"'}};
}

void Bracket::ByteUpdate() {
  probs_ = 1./256;
  if (active_.empty() || brackets_.find(byte_) != brackets_.end()) {
    if (brackets_.find(byte_) != brackets_.end()) {
      active_.push_back(byte_);
      distance_.push_back(0);
      if (active_.size() > 10) {
        active_.erase(active_.begin());
        distance_.erase(distance_.begin());
      }
      float p = stats_[byte_][0].first / stats_[byte_][0].second;
      probs_ = (1 - p) / 255;
      probs_[brackets_[byte_]] = p;
    }
  } else {
    int active = active_[active_.size() - 1];
    int distance = distance_[distance_.size() - 1];
    ++stats_[active][distance].second;
    if (brackets_[active] == byte_) {
      ++stats_[active][distance].first;
    }
    if (brackets_[active] == byte_ || distance >= max_ - 1) {
      active_.pop_back();
      distance_.pop_back();
      if (!active_.empty()) {
        int active = active_[active_.size() - 1];
        int distance = distance_[distance_.size() - 1];
        float p = stats_[active][distance].first /
            stats_[active][distance].second;
        probs_ = (1 - p) / 255;
        probs_[brackets_[active]] = p;
      }
    } else {
      ++distance_[distance_.size() - 1];
      ++distance;
      float p = stats_[active][distance].first /
          stats_[active][distance].second;
      probs_ = (1 - p) / 255;
      probs_[brackets_[active]] = p;
    }
  }
  ByteModel::ByteUpdate();
}
