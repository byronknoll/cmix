#ifndef MIXER_INPUT_H
#define MIXER_INPUT_H

#include "logistic.h"

#include <valarray>

class MixerInput {
 public:
  MixerInput(const Logistic& logistic, float eps);
  void SetNumModels(int num_models);
  void SetInput(int index, float p);

  std::valarray<float> inputs_;

 private:
  const Logistic& logistic_;
  float min_;
  float max_;
};

#endif
