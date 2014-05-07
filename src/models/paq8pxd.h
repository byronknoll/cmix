#ifndef PAQ8PXD_H
#define PAQ8PXD_H

#include "model.h"
#include <vector>

class PAQ8PXD : public Model {
 public:
  PAQ8PXD(int memory);
  float Predict();
  void Perceive(int bit);
  void ByteUpdate() {};
  const std::vector<float>& ModelPredictions();
};

#endif
