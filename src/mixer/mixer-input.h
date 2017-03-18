#ifndef MIXER_INPUT_H
#define MIXER_INPUT_H

#include "logistic.h"

#include <valarray>

class MixerInput {
 public:
  MixerInput(const Logistic& logistic, float eps);
  void SetNumModels(int num_models);
  void SetInput(int index, float p);
  void SetStretchedInput(int index, float p) { inputs_[index] = p; }
  const std::valarray<float>& Inputs() const { return inputs_; }

 private:
  std::valarray<float> inputs_;
  const Logistic& logistic_;
  float min_, max_;
};

#endif
