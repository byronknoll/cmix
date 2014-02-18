#include "ppm.h"

PPM::PPM(unsigned int order, const unsigned int& bit_context) :
    byte_(bit_context), cur_(NULL), cur_depth_(0), max_order_(order),
    top_(255), mid_(0), bot_(0) {
  cur_ = new Table();
  probs_.fill(1.0 / 256);
}

float PPM::Predict() {
  float num = 0, denom = 0;
  mid_ = bot_ + ((top_ - bot_) / 2);
  for (int i = bot_; i <= top_; ++i) {
    denom += probs_[i];
    if (i > mid_) num += probs_[i];
  }
  return num / denom;
}

void PPM::Perceive(int bit) {
  if (bit) {
    bot_ = mid_ + 1;
  } else {
    top_ = mid_;
  }
}

Table* PPM::AddOrGetTable(Table* cur, unsigned int order, unsigned char byte) {
  unsigned int index = 0;
  for (unsigned int i = 0; i < cur->entries.size(); ++i) {
    if (cur->entries[i].symbol == byte) {
      if (cur->entries[i].link != NULL) return cur->entries[i].link;
      index = i;
      break;
    }
  }
  Table* next = new Table();
  cur->entries[index].link = next;
  if (order == 0) {
    next->lower_table = cur;
  } else {
    next->lower_table = AddOrGetTable(cur->lower_table, order - 1, byte);
  }
  return next;
}

void PPM::UpdateTable(Table* cur, unsigned int depth, unsigned char byte) {
  for (unsigned int i = 0; i < cur->entries.size(); ++i) {
    if (cur->entries[i].symbol == byte) {
      ++(cur->entries[i].count);
      if (cur->entries[i].count == 255) {
        for (unsigned int j = 0; j < cur->entries.size(); ++j) {
          cur->entries[j].count /= 2;
        }
      }
      return;
    }
  }
  Entry add;
  add.symbol = byte;
  add.count = 1;
  add.link = NULL;
  cur->entries.push_back(add);
  if (cur->lower_table != NULL) UpdateTable(cur->lower_table, depth - 1, byte);
}

void PPM::ByteUpdate() {
  UpdateTable(cur_, cur_depth_, byte_);

  if (cur_depth_ == max_order_) {
    cur_ = cur_->lower_table;
    --cur_depth_;
  }
  cur_ = AddOrGetTable(cur_, cur_depth_, byte_);
  ++cur_depth_;

  Table* node = cur_;
  probs_.fill(0);
  float escape = 1;
  while (true) {
    int sum = 1;
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      if (probs_[node->entries[i].symbol] == 0) sum += node->entries[i].count;
    }
    float factor = escape / sum;
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      unsigned char symbol = node->entries[i].symbol;
      if (probs_[symbol] == 0) {
        probs_[symbol] = factor * node->entries[i].count;
      }
    }
    escape *= 1.0 / sum;
    node = node->lower_table;
    if (node == NULL) {
      factor = escape / 256;
      for (int i = 0; i < 256; ++i) {
        if (probs_[i] == 0) {
          probs_[i] = factor;
        }
      }
      break;
    }
  }
  top_ = 255;
  bot_ = 0;
}
