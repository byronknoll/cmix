#include "rubin-karp.h"

RubinKarp::RubinKarp(const unsigned int& bit_context, unsigned int order) :
    byte_(bit_context), order_(order), pow_(1), base_(29), rolling_(order, 0),
    pos_(0) {
  context_ = 0;
  size_ = 0;
  for (unsigned int i = 0; i < order_; ++i) {
    size_ = base_ * size_ + 255;
    pow_ *= base_;
  }
}

void RubinKarp::Update() {
  context_ = context_ * base_ + byte_;
  context_ -= pow_ * rolling_[pos_];
  rolling_[pos_] = byte_;
  ++pos_;
  if (pos_ == order_) pos_ = 0;
}

bool RubinKarp::IsEqual(Context* c) {
  RubinKarp* p = dynamic_cast<RubinKarp*>(c);
  if (!p) return false;
  if (order_ == p->order_) return true;
  return false;
}
