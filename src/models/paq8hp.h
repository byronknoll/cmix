#ifndef PAQ8HP_H
#define PAQ8HP_H

#include "model.h"
#include <vector>

class PAQ8HP : public Model {
 public:
  PAQ8HP(int memory);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {};
  const std::vector<float>& ModelPredictions();
};

#endif
