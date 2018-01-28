#ifndef PAQ8_H
#define PAQ8_H

#include "model.h"
#include <vector>

class PAQ8 : public Model {
 public:
  PAQ8(int memory);
  const std::valarray<float>& Predict();
  unsigned int NumOutputs();
  void Perceive(int bit);
  void ByteUpdate() {};
};

#endif
