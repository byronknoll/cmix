#include "dmc.h"

DMC::DMC(float delta, unsigned int max_size) : cur_(0), delta_(delta),
    max_size_(max_size), threshold_(2), big_threshold_(2) {
  Reset();
}

void DMC::Reset() {
  nodes_.clear();
  nodes_.resize(256 * 256, Node());
  for (int i = 0; i < 256; ++i) {
    for (int j = 0; j < 256; ++j) {
      Node& n = nodes_[j * 256 + i];
      n.count[0] = 1.5;
      n.count[1] = 1.5;
      if (i < 127) {
        n.next[0] = 256 * j + 2 * i + 1;
        n.next[1] = 256 * j + 2 * i + 2;
      } else {
        n.next[0] = (i + 1) * 256;
        n.next[1] = (i - 127) * 256;
      }
    }
  }
  cur_ = 0;
}

float DMC::Predict() {
  return (delta_ + nodes_[cur_].count[1]) / (nodes_[cur_].count[0] +
      nodes_[cur_].count[1] + 2 * delta_);
}

void DMC::Perceive(int bit) {
  Node& cur = nodes_[cur_];
  Node& next = nodes_[cur.next[bit]];
  if (cur.count[bit] >= threshold_ &&
      next.count[0] + next.count[1] >= big_threshold_ + cur.count[bit]) {
    Node add;
    float r = cur.count[bit] / (next.count[0] + next.count[1]);
    add.count[0] = next.count[0] * r;
    next.count[0] -= add.count[0];
    add.count[1] = next.count[1] * r;
    next.count[1] -= add.count[1];
    add.next[0] = next.next[0];
    add.next[1] = next.next[1];
    cur.next[bit] = nodes_.size();
    nodes_.push_back(add);
    if (nodes_.size() > max_size_) Reset();
  }
  nodes_[cur_].count[bit]++;
  cur_ = nodes_[cur_].next[bit];
}
