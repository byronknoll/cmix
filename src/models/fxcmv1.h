#ifndef FXCM_H
#define FXCM_H

#include "model.h"
#include <vector>
#include <memory>

namespace fxcmv1 {
  class Predictor;
}

class FXCM : public Model {
 public:
  FXCM(int memory);
  const std::valarray<float>& Predict();
  unsigned int NumOutputs();
  void Perceive(int bit);
  void ByteUpdate() {};

 private:
  std::unique_ptr<fxcmv1::Predictor> predictor_;
};

#endif
