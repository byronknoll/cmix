#include "bracket.h"

Bracket::Bracket(const unsigned int& bit_context, int distance_limit,
    int stack_limit, int stats_limit, const std::vector<bool>& vocab) :
    ByteModel(vocab), distance_limit_(distance_limit),
    stack_limit_(stack_limit), stats_limit_(stats_limit), byte_(bit_context),
    stats_(256, std::vector<std::pair<unsigned int, unsigned int>>
    (distance_limit_, {1, 256})) {
  brackets_ = {{'(',')'}, {'{','}'}, {'[',']'}, {'<','>'}, {'\'','\''},
      {'"','"'}};
}

void Bracket::ByteUpdate() {
  probs_ = 1./256;
  if (active_.empty() || (brackets_.find(byte_) != brackets_.end() &&
      !(active_[active_.size() - 1] == byte_ && brackets_[byte_] == byte_))) {
    if (brackets_.find(byte_) != brackets_.end()) {
      active_.push_back(byte_);
      distance_.push_back(0);
      if (active_.size() > stack_limit_) {
        active_.erase(active_.begin());
        distance_.erase(distance_.begin());
      }
      float p = (1. * stats_[byte_][0].first) / stats_[byte_][0].second;
      probs_ = (1 - p) / 255;
      probs_[brackets_[byte_]] = p;
    }
  } else {
    unsigned int active = active_[active_.size() - 1];
    unsigned int distance = distance_[distance_.size() - 1];
    ++stats_[active][distance].second;
    if (brackets_[active] == byte_) {
      ++stats_[active][distance].first;
    }
    if (stats_[active][distance].second > stats_limit_) {
      stats_[active][distance].first /= 2;
      stats_[active][distance].second /= 2;
    }
    if (brackets_[active] == byte_ || distance >= distance_limit_ - 1) {
      active_.pop_back();
      distance_.pop_back();
      if (!active_.empty()) {
        int active = active_[active_.size() - 1];
        int distance = distance_[distance_.size() - 1];
        float p = (1. * stats_[active][distance].first) /
            stats_[active][distance].second;
        probs_ = (1 - p) / 255;
        probs_[brackets_[active]] = p;
      }
    } else {
      ++distance_[distance_.size() - 1];
      ++distance;
      float p = (1. * stats_[active][distance].first) /
          stats_[active][distance].second;
      probs_ = (1 - p) / 255;
      probs_[brackets_[active]] = p;
    }
  }
  ByteModel::ByteUpdate();
}
