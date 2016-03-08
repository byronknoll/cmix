#ifndef BYTE_MODEL_H
#define BYTE_MODEL_H

#include "model.h"

#include <valarray>

class ByteModel : public Model {
 public:
  virtual ~ByteModel() {}
  ByteModel();
  const std::valarray<float>& BytePredict();
  float Predict();
  void Perceive(int bit);
  virtual void ByteUpdate();

 protected:
  int top_, mid_, bot_;
  std::valarray<float> probs_;
};

#endif
