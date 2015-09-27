#ifndef BYTE_MODEL_H
#define BYTE_MODEL_H

#include "model.h"

#include <array>

class ByteModel : Model {
 public:
  virtual ~ByteModel() {}
  ByteModel();
  std::array<float, 256> BytePredict();
  float Predict();
  void Perceive(int bit);
  virtual void ByteUpdate();

 protected:
  std::array<float, 256> probs_;
  int top_, mid_, bot_;
};

#endif
