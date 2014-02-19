#include "ppm.h"

PPM::PPM(unsigned int order, const unsigned int& bit_context,
    unsigned int max_size) : byte_(bit_context), cur_(0), cur_depth_(0),
    max_order_(order), max_size_(max_size), top_(255), mid_(0), bot_(0) {
  Reset();
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

int PPM::AddOrGetTable(int table_index, unsigned int order,
    unsigned char byte) {
  Table& cur = tables_[table_index];
  unsigned int index = 0;
  for (unsigned int i = 0; i < cur.entries.size(); ++i) {
    if (cur.entries[i].symbol == byte) {
      if (cur.entries[i].link != -1) return cur.entries[i].link;
      index = i;
      break;
    }
  }
  tables_.push_back(Table());
  int next = tables_.size() - 1;
  tables_[table_index].entries[index].link = next;
  tables_[next].escape = 1;
  if (order == 0) {
    tables_[next].lower_table = 0;
  } else {
    int temp = AddOrGetTable(tables_[table_index].lower_table, order - 1, byte);
    tables_[next].lower_table = temp;
  }
  return next;
}

void PPM::UpdateTable(int table_index, unsigned int depth, unsigned char byte) {
  Table& cur = tables_[table_index];
  for (unsigned int i = 0; i < cur.entries.size(); ++i) {
    if (cur.entries[i].symbol == byte) {
      ++(cur.entries[i].count);
      if (cur.entries[i].count == 50) {
        for (unsigned int j = 0; j < cur.entries.size(); ++j) {
          cur.entries[j].count /= 2;
        }
        cur.escape /= 2;
        if (depth == 0) {
          for (unsigned int j = 0; j < cur.entries.size(); ++j) {
            if (cur.entries[j].count == 0) cur.entries[j].count = 1;
          }
        } else if (cur.escape == 0) cur.escape = 1;
      }
      return;
    }
  }
  Entry add;
  add.symbol = byte;
  add.count = 1;
  add.link = -1;
  cur.entries.push_back(add);
  ++cur.escape;
  if (cur.lower_table != -1) UpdateTable(cur.lower_table, depth - 1, byte);
}

void PPM::Reset() {
  tables_.clear();
  Table t;
  t.lower_table = -1;
  t.escape = 0;
  tables_.push_back(t);
  cur_ = 0;
  cur_depth_ = 0;
}

void PPM::ByteUpdate() {
  if (tables_.size() > max_size_) Reset();

  UpdateTable(cur_, cur_depth_, byte_);

  if (cur_depth_ == max_order_) {
    cur_ = tables_[cur_].lower_table;
    --cur_depth_;
  }
  cur_ = AddOrGetTable(cur_, cur_depth_, byte_);
  ++cur_depth_;

  Table* node = &tables_[cur_];
  probs_.fill(0);
  float escape = 1;
  while (true) {
    float sum = node->escape;
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      if (probs_[node->entries[i].symbol] == 0) sum += node->entries[i].count;
    }
    float factor = escape / sum;
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      if (probs_[node->entries[i].symbol] == 0) {
        probs_[node->entries[i].symbol] = factor * node->entries[i].count;
      }
    }
    if (node->lower_table == -1) break;
    escape *= (1.0 * node->escape) / sum;
    node = &tables_[node->lower_table];
  }
  top_ = 255;
  bot_ = 0;
}
