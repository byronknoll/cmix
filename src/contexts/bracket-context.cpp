#include "bracket-context.h"

BracketContext::BracketContext(const unsigned int& bit_context,
    int distance_limit, int stack_limit) : byte_(bit_context),
    distance_limit_(distance_limit), stack_limit_(stack_limit) {
  brackets_ = {{'(',')'}, {'{','}'}, {'[',']'}, {'<','>'}};
  context_ = 0;
  size_ = 257 * distance_limit_;
}

void BracketContext::Update() {
  if (!active_.empty()) {
    if (brackets_[active_[active_.size() - 1]] == byte_ ||
        distance_[distance_.size() - 1] >= distance_limit_ - 1) {
      active_.pop_back();
      distance_.pop_back();
    } else {
      ++distance_[distance_.size() - 1];
    }
  }
  if (brackets_.find(byte_) != brackets_.end()) {
    active_.push_back(byte_);
    distance_.push_back(0);
    if (brackets_.size() > stack_limit_) {
      active_.erase(active_.begin());
      distance_.erase(distance_.begin());
    }
  }
  if (!active_.empty()) {
    context_ = distance_limit_ * (active_[active_.size() - 1] + 1) +
        distance_[distance_.size() - 1];
  } else {
    context_ = 0;
  }
}

bool BracketContext::IsEqual(Context* c) {
  BracketContext* p = dynamic_cast<BracketContext*>(c);
  if (!p) return false;
  if (distance_limit_ == p->distance_limit_ && stack_limit_ == p->stack_limit_)
    return true;
  return false;
}
