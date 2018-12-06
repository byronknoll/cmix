#ifndef PAQ8_H
#define PAQ8_H

#include "model.h"
#include <vector>
#include <memory>

namespace paq8 {
  class Predictor;
}

class PAQ8 : public Model {
 public:
  PAQ8(int memory);
  const std::valarray<float>& Predict();
  unsigned int NumOutputs();
  void Perceive(int bit);
  void ByteUpdate() {};

 private:
  std::unique_ptr<paq8::Predictor> predictor_;
};

#endif
