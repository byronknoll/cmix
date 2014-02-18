#include "ppm.h"

PPM::PPM(unsigned int order, const unsigned int& bit_context) :
    bit_context_(bit_context), cur_(NULL), cur_depth_(0), max_order_(order),
    top_(255), bot_(0) {
  cur_ = new Table();
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
  unsigned int index = 0;
  for (unsigned int i = 0; i < cur->entries.size(); ++i) {
    if (cur->entries[i].symbol == byte) {
      if (cur->links[i] != NULL) return cur->links[i];
      index = i;
      break;
    }
  }
  Table* next = new Table();
  cur->links[index] = next;
  if (order == 0) {
    next->lower_table = cur;
  } else {
    next->lower_table = AddOrGetTable(cur->lower_table, order - 1, byte);
  }
  return next;
}

void PPM::UpdateTable(Table* cur, unsigned int depth, unsigned char byte) {
  bool found = false;
  for (unsigned int i = 0; i < cur->entries.size(); ++i) {
    if (cur->entries[i].symbol == byte) {
      found = true;
      ++(cur->entries[i].count);
      if (cur->entries[i].count == 255) {
        for (unsigned int j = 0; j < cur->entries.size(); ++j) {
          cur->entries[j].count /= 2;
        }
      }
      break;
    }
  }
  if (!found) {
    Entry add;
    add.symbol = byte;
    add.count = 1;
    cur->entries.push_back(add);
    if (depth != max_order_) cur->links.push_back(NULL);
  }
  if (cur->lower_table != NULL) UpdateTable(cur->lower_table, depth - 1, byte);
}

void PPM::ByteUpdate() {
  UpdateTable(cur_, cur_depth_, bit_context_);

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
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      if (probs_[node->entries[i].symbol] == 0) sum += node->entries[i].count;
    }
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      unsigned char symbol = node->entries[i].symbol;
      if (probs_[symbol] == 0) {
        probs_[symbol] = (escape * node->entries[i].count) / sum;
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
