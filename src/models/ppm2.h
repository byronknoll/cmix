#ifndef PPM2_H
#define PPM2_H

#include "byte-model.h"

#include <vector>
#include <array>

struct Entry2 {
  unsigned char count;
  unsigned char symbol;
  unsigned int link;
};

struct Table2 {
  std::vector<Entry2> entries;
};

class PPM2 : public ByteModel {
 public:
  PPM2(unsigned int order, const unsigned int& bit_context,
      unsigned int max_size);
  void ByteUpdate();

 private:
  void Reset();
  void UpdateTables();
  void UpdatePrediction();
  void SetRecentByte();
  unsigned char GetRecentByte(int age);

  const unsigned int& byte_;
  std::vector<Table2> tables_;
  std::vector<unsigned char> recent_bytes_;
  unsigned int recent_bytes_pos_, max_order_, max_size_;
  std::vector<unsigned int> matches_;
};

#endif
