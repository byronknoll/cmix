#ifndef PPMD_H
#define PPMD_H

#include "byte-model.h"

class PPMD : public ByteModel {
 public:
  PPMD(int order, int memory, const unsigned int& bit_context,
      const std::vector<bool>& vocab);
  void ByteUpdate();
 private:
  const unsigned int& byte_;
};

#endif
