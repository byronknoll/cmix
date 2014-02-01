#include "mixer.h"

Mixer::Mixer(const std::vector<float>& inputs, const Logistic& logistic,
    const unsigned long long& byte_context, const unsigned int& bit_context,
    float learning_rate, unsigned long long table_size) :
    inputs_(inputs), logistic_(logistic), p_(0.5),
    learning_rate_(learning_rate), bit_context_(bit_context),
    byte_context_(byte_context),
    weights_(table_size, std::array<std::vector<float>, 256>()) {}

void Mixer::SetNumModels(int num_models) {
  for (unsigned int i = 0; i < weights_.size(); ++i) {
    for (int j = 0; j < 256; ++j) {
      weights_[i][j].resize(num_models, 0);
    }
  }
}

float Mixer::Mix() {
  p_ = 0;
  for (unsigned int i = 0; i < inputs_.size(); ++i) {
    p_ += weights_[byte_context_][bit_context_][i] * inputs_[i];
  }
  p_ = logistic_.Squash(p_);
  return p_;
}

void Mixer::Perceive(int bit) {
  float update = learning_rate_ * (bit - p_);
  for (unsigned int i = 0; i < inputs_.size(); ++i) {
    weights_[byte_context_][bit_context_][i] += update * inputs_[i];
  }
}
