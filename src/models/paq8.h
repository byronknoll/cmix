#ifndef PAQ8_H
#define PAQ8_H

#include "model.h"
#include <vector>

class PAQ8 : public Model {
 public:
  PAQ8(int memory);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {};
  const std::vector<float>& ModelPredictions();
};

#endif
