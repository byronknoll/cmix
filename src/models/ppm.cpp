#include "ppm.h"

PPM::PPM(unsigned int order, const unsigned int& bit_context) :
    bit_context_(bit_context), cur_(NULL), cur_depth_(0), max_order_(order),
    top_(255), bot_(0) {
  cur_ = new Table();
  cur_->total = 1;
  cur_->links = new Links();
  probs_.fill(1.0 / 256);
}

float PPM::Predict() {
  float num = 0, denom = 0;
  int mid = bot_ + ((top_ - bot_) / 2);
  for (int i = bot_; i <= top_; ++i) {
    denom += probs_[i];
    if (i > mid) num += probs_[i];
  }
  return num / denom;
}

void PPM::Perceive(int bit) {
  int mid = bot_ + ((top_ - bot_) / 2);
  if (bit) {
    bot_ = mid + 1;
  } else {
    top_ = mid;
  }
}

Table* PPM::AddOrGetTable(Table* cur, unsigned int order, unsigned char byte) {
  if (cur->links->link[byte] != NULL) return cur->links->link[byte];
  Table* next = new Table();
  next->total = 1;
  cur->links->link[byte] = next;
  if (order < max_order_ - 1) next->links = new Links();
  if (order == 0) {
    next->lower_table = cur;
  } else {
    next->lower_table = AddOrGetTable(cur->lower_table, order - 1, byte);
  }
  return next;
}

void PPM::UpdateTable(Table* cur, unsigned char byte) {
  if (cur->counts[byte] == 255) {
    cur->total = 1;
    for (int i = 0; i < 256; ++i) {
      cur->counts[i] /= 2;
      cur->total += cur->counts[i];
    }
  }
  ++cur->counts[byte];
  ++cur->total;
  if (cur->lower_table != NULL) UpdateTable(cur->lower_table, byte);
}

void PPM::ByteUpdate() {
  UpdateTable(cur_, bit_context_);

  if (cur_depth_ == max_order_) {
    cur_ = cur_->lower_table;
    --cur_depth_;
  }
  cur_ = AddOrGetTable(cur_, cur_depth_, bit_context_);
  ++cur_depth_;

  Table* node = cur_;
  probs_.fill(0);
  double escape = 1;
  while (true) {
    int sum = 1;
    for (int i = 0; i < 256; ++i) {
      if (probs_[i] == 0) sum += node->counts[i];
    }
    for (int i = 0; i < 256; ++i) {
      if (probs_[i] == 0 && node->counts[i] > 0) {
        probs_[i] = (escape * node->counts[i]) / sum;
      }
    }
    bool done = true;
    for (int i = 0; i < 256; ++i) {
      if (probs_[i] == 0) {
        done = false;
        break;
      }
    }
    if (done) break;
    escape *= 1.0 / sum;
    node = node->lower_table;
    if (node == NULL) {
      for (int i = 0; i < 256; ++i) {
        if (probs_[i] == 0) {
          probs_[i] += escape / 256;
        }
      }
      break;
    }
  }
  top_ = 255;
  bot_ = 0;
}
