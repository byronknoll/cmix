#include "mixer.h"

Mixer::Mixer(const std::vector<float>& inputs, const Logistic& logistic,
    const unsigned long long& context, float learning_rate,
    unsigned long long context_size, unsigned long long input_size) :
    inputs_(inputs), logistic_(logistic), p_(0.5),
    learning_rate_(learning_rate), context_(context),
    weights_(context_size, std::vector<float>(input_size, 0)) {}

float Mixer::Mix() {
  p_ = 0;
  for (unsigned int i = 0; i < inputs_.size(); ++i) {
    p_ += weights_[context_][i] * inputs_[i];
  }
  p_ = logistic_.Squash(p_);
  return p_;
}

void Mixer::Perceive(int bit) {
  float update = learning_rate_ * (bit - p_);
  for (unsigned int i = 0; i < inputs_.size(); ++i) {
    weights_[context_][i] += update * inputs_[i];
  }
}

unsigned long long Mixer::GetNumNeurons() {
  return weights_.size();
}
