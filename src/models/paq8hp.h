#ifndef PAQ8HP_H
#define PAQ8HP_H

#include "model.h"
#include <vector>

class PAQ8HP : public Model {
 public:
  PAQ8HP(int memory);
  const std::valarray<float>& Predict();
  unsigned int NumOutputs();
  void Perceive(int bit);
  void ByteUpdate() {};
};

#endif
