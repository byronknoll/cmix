#include "ppm2.h"
#include <math.h>

PPM2::PPM2(unsigned int order, const unsigned int& bit_context,
    unsigned int max_size) : byte_(bit_context), recent_bytes_(order, 0),
    recent_bytes_pos_(0), max_order_(order), max_size_(max_size), top_(255),
    mid_(0), bot_(0) {
  Reset();
  probs_.fill(1.0 / 256);
}

float PPM2::Predict() {
  float num = 0, denom = 0;
  mid_ = bot_ + ((top_ - bot_) / 2);
  for (int i = bot_; i <= top_; ++i) {
    denom += probs_[i];
    if (i > mid_) num += probs_[i];
  }
  return num / denom;
}

void PPM2::Perceive(int bit) {
  if (bit) {
    bot_ = mid_ + 1;
  } else {
    top_ = mid_;
  }
}

void PPM2::Reset() {
  tables_.clear();
  tables_.push_back(Table2());
  for (int i = 0; i < 256; ++i) {
    Entry2 e;
    e.count = 1;
    e.symbol = i;
    e.link = tables_.size();
    tables_[0].entries.push_back(e);
    tables_.push_back(Table2());
  }
}

void PPM2::UpdateTables() {
  unsigned int cur_context = 0;
  for (unsigned int order = 0; order < max_order_; ++order) {
    unsigned int next_context = 0;
    bool found = false;
    Table2& cur = tables_[cur_context];
    unsigned char recent_byte = GetRecentByte(order);
    for (unsigned int i = 0; i < cur.entries.size(); ++i) {
      if (cur.entries[i].symbol == byte_) {
        found = true;
        ++cur.entries[i].count;
        if (cur.entries[i].count == 255) {
          for (unsigned int j = 0; j < cur.entries.size(); ++j) {
            if (cur.entries[j].count > 1) cur.entries[j].count /= 2;
          }
        }
        if (next_context != 0) break;
      }
      if (cur.entries[i].symbol == recent_byte) {
        next_context = cur.entries[i].link;
        if (found) break;
      }
    }
    if (!found) {
      Entry2 e;
      e.count = 1;
      e.symbol = byte_;
      if (byte_ == recent_byte) {
        next_context = tables_.size();
      }
      e.link = tables_.size();
      cur.entries.push_back(e);
      tables_.push_back(Table2());
    }
    if (next_context == 0) {
      Entry2 e;
      e.count = 0;
      e.symbol = recent_byte;
      e.link = tables_.size();
      next_context = tables_.size();
      tables_[cur_context].entries.push_back(e);
      tables_.push_back(Table2());
    }
    cur_context = next_context;
  }
}

void PPM2::UpdatePrediction() {
  for (int i = 0; i < 256; ++i) {
    probs_[i] = 0;
  }

  matches_.clear();
  unsigned int cur_context = 0;  
  for (unsigned int order = 0; order < max_order_; ++order) {
    unsigned int next_context = 0;
    Table2& cur = tables_[cur_context];
    unsigned char recent_byte = GetRecentByte(order);
    for (unsigned int i = 0; i < cur.entries.size(); ++i) {
      if (cur.entries[i].symbol == recent_byte) {
        next_context = cur.entries[i].link;
      }
    }
    matches_.push_back(cur_context);
    if (next_context == 0) break;
    cur_context = next_context;
  }

  float w1 = 1;
  float param1 = 4.4, param2 = 0.4;
  for (int order = matches_.size() - 1; order >= 0; --order) {
    Table2& cur = tables_[matches_[order]];
    float a = 0;
    float b = 0;
    for (unsigned int i = 0; i < cur.entries.size(); ++i) {
      if (cur.entries[i].count > 0) {
        a += cur.entries[i].count;
        ++b;
      }
    }
    float w2 = w1 * pow((a / (a + b)), param1);
    for (unsigned int i = 0; i < cur.entries.size(); ++i) {
      if (cur.entries[i].count > 0) {
        float p = (cur.entries[i].count / a) * w2;
        probs_[cur.entries[i].symbol] += p;
      }
    }
    if (a != 0) w1 = param2 * w1 + (1 - param2) * ((w1 * b) / (a + b));
  }
}

void PPM2::SetRecentByte() {
  recent_bytes_[recent_bytes_pos_] = byte_;
  ++recent_bytes_pos_;
  if (recent_bytes_pos_ == max_order_) {
    recent_bytes_pos_ = 0;
  }
}

unsigned char PPM2::GetRecentByte(int age) {
  int index = recent_bytes_pos_ - age - 1;
  if (index < 0) index += max_order_;
  return recent_bytes_[index];
}

void PPM2::ByteUpdate() {
  if (tables_.size() > max_size_) Reset();

  UpdateTables();
  SetRecentByte();
  UpdatePrediction();

  top_ = 255;
  bot_ = 0;
}
