#ifndef PPM_H
#define PPM_H

#include "byte-model.h"

#include <vector>
#include <array>

struct Table;

struct Entry {
  unsigned char count;
  unsigned char symbol;
  int link;
};

struct Table {
  Table() : entries(0, Entry()) {}
  std::vector<Entry> entries;
  int lower_table;
  unsigned char escape;
};

class PPM : public ByteModel {
 public:
  PPM(unsigned int order, const unsigned int& bit_context, float delta,
      unsigned int max_size);
  void ByteUpdate();

 private:
  int AddOrGetTable(int table_index, unsigned int depth, unsigned char byte);
  void UpdateTable(int table_index, unsigned int depth, unsigned char byte);
  void Reset();
  int EscapeContext(int escape_count, int match_count, int order);

  const unsigned int& byte_;
  float divisor_;
  std::vector<Table> tables_;
  std::vector<float> escape_map_;
  unsigned int cur_, cur_depth_, max_order_, max_size_;
};

#endif
