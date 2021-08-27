#ifndef PPMD_H
#define PPMD_H

#include "byte-model.h"

#include <memory>

namespace PPMD {

struct ppmd_Model;

class PPMD : public ByteModel {
 public:
  PPMD(int order, int memory, const unsigned int& bit_context,
      const std::vector<bool>& vocab);
  void ByteUpdate();
 private:
  const unsigned int& byte_;
  std::unique_ptr<ppmd_Model> ppmd_model_;
};

} // namespace PPMD

#endif
