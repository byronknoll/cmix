#include "mixer-input.h"

#include <math.h>

MixerInput::MixerInput(float eps, int stretch_size, int squash_size) :
    inputs_(1, 0.5), min_(eps), max_(1 - eps), stretch_size_(stretch_size),
    squash_size_(squash_size), stretch_table_(stretch_size, 0),
    squash_table_(squash_size, 0) {
  for (int i = 0; i < stretch_size_; ++i) {
    stretch_table_[i] = SlowStretch((i + 0.5) / stretch_size_);
  }
  for (int i = 0; i < squash_size_; ++i) {
    squash_table_[i] = SlowSquash(100 * ((i + 0.5) / squash_size_) - 50);
  }
}

void MixerInput::SetNumModels(int num_models) {
  inputs_.resize(num_models, 0.5);
}

void MixerInput::SetInput(int index, float p) {
  if (p < min_) p = min_;
  else if (p > max_) p = max_;
  inputs_[index] = Stretch(p);
}

float MixerInput::Stretch(float p) const {
  int index = p * stretch_size_;
  return stretch_table_[index];
}

float MixerInput::Squash(float p) const {
  int index = ((p + 50) / 100) * squash_size_;
  if (index >= squash_size_) index = squash_size_ - 1;
  else if (index < 0) index = 0;
  return squash_table_[index];
}

float MixerInput::SlowStretch(float p) {
  return log(p / (1 - p));
}

float MixerInput::SlowSquash(float p) {
  return 1 / (1 + exp(-p));
}
