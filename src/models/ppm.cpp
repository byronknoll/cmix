#include "ppm.h"

PPM::PPM(unsigned int order, const unsigned int& bit_context, float delta,
    unsigned int max_size) : byte_(bit_context), divisor_(1.0 / delta),
    escape_map_(256 * 256, 0), cur_(0), cur_depth_(0), max_order_(order),
    max_size_(max_size), top_(255), mid_(0), bot_(0) {
  Reset();
  probs_.fill(1.0 / 256);
  for (int escape = 1; escape < 128; ++escape) {
    for (int sum = 0; sum < 256; ++sum) {
      for (unsigned int order = 0; order <= max_order_; ++order) {
        float prob = (1.0 * escape) / (escape + sum);
        int context = EscapeContext(escape, sum, order);
        escape_map_[context] = prob;
      }
    }
  }
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
  int sum = 0;
  for (unsigned int i = 0; i < cur.entries.size(); ++i) {
    sum += cur.entries[i].count;
    if (cur.entries[i].symbol == byte) {
      cur.entries[i].count += 2;
      if (cur.entries[i].count >= 254) {
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
      for (unsigned int j = i + 1; j < cur.entries.size(); ++j) {
        sum += cur.entries[j].count;
      }
      int escape_context = EscapeContext(cur.escape, sum, depth);
      escape_map_[escape_context] -= escape_map_[escape_context] * divisor_;

      if (cur.lower_table != -1) {
        Table& low = tables_[cur.lower_table];
        for (unsigned int j = 0; j < low.entries.size(); ++j) {
          if (low.entries[j].symbol == byte) {
            ++(low.entries[j].count);
          }
        }
      }

      return;
    }
  }
  Entry add;
  add.symbol = byte;
  add.count = 1;
  add.link = -1;
  cur.entries.push_back(add);
  int escape_context = EscapeContext(cur.escape, sum, depth);
  escape_map_[escape_context] += (1 - escape_map_[escape_context]) * divisor_;
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

int PPM::EscapeContext(int escape_count, int match_count, int order) {
  if (match_count > 255) match_count = 255;
  if (escape_count > 127) escape_count = 127;
  if (order > 3) order = 1;
  else order = 0;
  return order + escape_count * 2 + match_count * 256;
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
  int order = cur_depth_;
  while (true) {
    int sum = 0;
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      if (probs_[node->entries[i].symbol] == 0) sum += node->entries[i].count;
    }
    float escape_prob = escape_map_[EscapeContext(node->escape, sum, order)];
    float factor = (1 - escape_prob) * escape / sum;
    for (unsigned int i = 0; i < node->entries.size(); ++i) {
      if (probs_[node->entries[i].symbol] == 0) {
        probs_[node->entries[i].symbol] = factor * node->entries[i].count;
      }
    }
    if (node->lower_table == -1) break;
    escape *= escape_prob;
    node = &tables_[node->lower_table];
    --order;
  }
  top_ = 255;
  bot_ = 0;
}
