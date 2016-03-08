#include "bracket-context.h"

BracketContext::BracketContext(const unsigned int& bit_context) :
    byte_(bit_context), max_(200), active_(0), distance_(0) {
  brackets_ = {{'(',')'}, {'{','}'}, {'[',']'}, {'<','>'}};
  context_ = 0;
  size_ = 257 * max_;
}

void BracketContext::Update() {
  if (active_ != -1) {
    if (brackets_[active_] == byte_) {
      active_ = -1;
    } else if (distance_ < max_ - 1) {
      ++distance_;
    }
  }
  if (brackets_.find(byte_) != brackets_.end()) {
    active_ = byte_;
    distance_ = 0;
  }
  if (active_ != -1) {
    context_ = max_ * (active_ + 1) + distance_;
  } else {
    context_ = 0;
  }
}

bool BracketContext::IsEqual(Context* c) {
  BracketContext* p = dynamic_cast<BracketContext*>(c);
  if (!p) return false;
  return false;
}
