#ifndef MIXER_INPUT_H
#define MIXER_INPUT_H

#include "sigmoid.h"

#include <valarray>
#include <vector>

class MixerInput {
 public:
  MixerInput(const Sigmoid& sigmoid, float eps);
  void SetNumModels(int num_models);
  void SetInput(int index, float p);
  void SetStretchedInput(int index, float p);
  void SetExtraInput(float p);
  void ClearExtraInputs() { extra_inputs_.clear(); }
  const std::valarray<float>& Inputs() const { return inputs_; }
  const std::vector<float>& ExtraInputs() const { return extra_inputs_; }

 private:
  std::valarray<float> inputs_;
  std::vector<float> extra_inputs_;
  const Sigmoid& sigmoid_;
  float min_, max_, stretched_min_, stretched_max_;
};

#endif
