#ifndef PPM_H
#define PPM_H

#include "model.h"

#include <vector>
#include <array>

struct Links;

struct Table {
  std::array<unsigned char, 256> counts;
  unsigned short total;
  Links* links;
  Table* lower_table;
};

struct Links {
  std::array<Table*, 256> link;
};

class PPM : public Model {
 public:
  PPM(unsigned int order, const unsigned int& bit_context);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate();

 private:
  Table* AddOrGetTable(Table* cur, unsigned int order, unsigned char byte);
  void UpdateTable(Table* cur, unsigned char byte);

  const unsigned int& bit_context_;
  Table* cur_;
  unsigned int cur_depth_;
  unsigned int max_order_;
  std::array<double, 256> probs_;
  int top_;
  int bot_;
};

#endif
