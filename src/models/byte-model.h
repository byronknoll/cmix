#ifndef BYTE_MODEL_H
#define BYTE_MODEL_H

#include "model.h"

#include <valarray>
#include <vector>

class ByteModel : public Model {
 public:
  virtual ~ByteModel() {}
  ByteModel(const std::vector<bool>& vocab);
  const std::valarray<float>& BytePredict();
  const std::valarray<float>& Predict();
  void Perceive(int bit);
  virtual void ByteUpdate();

 protected:
  int top_, mid_, bot_;
  const std::vector<bool>& vocab_;
  std::valarray<float> probs_;
};

#endif
