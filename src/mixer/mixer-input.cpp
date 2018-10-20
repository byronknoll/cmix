#include "mixer-input.h"

MixerInput::MixerInput(const Sigmoid& sigmoid, float eps) :
    inputs_(0.5, 1), sigmoid_(sigmoid), min_(eps), max_(1 - eps),
    stretched_min_(sigmoid.Logit(0)), stretched_max_(sigmoid.Logit(1)) {}

void MixerInput::SetNumModels(int num_models) {
  inputs_.resize(num_models, 0.5);
}

void MixerInput::SetInput(int index, float p) {
  if (p < min_) p = min_;
  else if (p > max_) p = max_;
  inputs_[index] = sigmoid_.Logit(p);
}

void MixerInput::SetStretchedInput(int index, float p) {
  if (p > stretched_max_) p = stretched_max_;
  else if (p < stretched_min_) p = stretched_min_;
  inputs_[index] = p;
}

void MixerInput::SetExtraInput(float p) {
  if (p > stretched_max_) p = stretched_max_;
  else if (p < stretched_min_) p = stretched_min_;
  extra_inputs_.push_back(p);
}
