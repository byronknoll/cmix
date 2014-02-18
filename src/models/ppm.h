#ifndef PPM_H
#define PPM_H

#include "model.h"

#include <vector>
#include <array>

struct Table;

struct Entry {
  unsigned char count;
  unsigned char symbol;
  Table* link;
};

struct Table {
  Table() : entries(0, Entry()) {}
  std::vector<Entry> entries;
  Table* lower_table;
};

class PPM : public Model {
 public:
  PPM(unsigned int order, const unsigned int& bit_context);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate();

 private:
  Table* AddOrGetTable(Table* cur, unsigned int depth, unsigned char byte);
  void UpdateTable(Table* cur, unsigned int depth, unsigned char byte);

  const unsigned int& byte_;
  Table* cur_;
  unsigned int cur_depth_;
  unsigned int max_order_;
  std::array<float, 256> probs_;
  int top_;
  int mid_;
  int bot_;
};

#endif
