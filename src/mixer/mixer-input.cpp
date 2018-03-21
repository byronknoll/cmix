#include "mixer-input.h"

MixerInput::MixerInput(const Sigmoid& sigmoid, float eps) :
    inputs_(0.5, 1), sigmoid_(sigmoid), min_(eps), max_(1 - eps) {}

void MixerInput::SetNumModels(int num_models) {
  inputs_.resize(num_models, 0.5);
}

void MixerInput::SetInput(int index, float p) {
  if (p < min_) p = min_;
  else if (p > max_) p = max_;
  inputs_[index] = sigmoid_.Logit(p);
}
