#include "mixer.h"

#include <numeric>

Mixer::Mixer(const std::valarray<float>& inputs, const Logistic& logistic,
    const unsigned long long& context, float learning_rate,
    unsigned long long context_size, unsigned long long input_size) :
    inputs_(inputs), logistic_(logistic), p_(0.5),
    learning_rate_(learning_rate), context_(context),
    weights_(context_size, std::valarray<float>(0.0, input_size)) {}

float Mixer::Mix() {
  p_ = logistic_.Squash(std::inner_product(&inputs_[0],
      &inputs_[inputs_.size()], &weights_[context_][0], 0.0));
  return p_;
}

void Mixer::Perceive(int bit) {
  float update = learning_rate_ * (bit - p_);
  weights_[context_] += update * inputs_;
}

unsigned long long Mixer::GetNumNeurons() {
  return weights_.size();
}
