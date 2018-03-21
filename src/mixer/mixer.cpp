#include "mixer.h"

#include "sigmoid.h"

#include <numeric>
#include <math.h>

Mixer::Mixer(const std::valarray<float>& inputs,
    const unsigned long long& context, float learning_rate,
    unsigned long long input_size) : inputs_(inputs), p_(0.5),
    learning_rate_(learning_rate), context_(context), max_steps_(1), steps_(0),
    input_size_(input_size) {}

float Mixer::Mix() {
  ContextData* data = context_map_[context_].get();
  if (data == nullptr) {
    context_map_[context_] = std::unique_ptr<ContextData>(
        new ContextData(input_size_));
    data = context_map_[context_].get();
  }
  p_ = (inputs_ * data->weights).sum();
  return p_;
}

void Mixer::Perceive(int bit) {
  ContextData* data = context_map_[context_].get();
  if (data == nullptr) {
    context_map_[context_] = std::unique_ptr<ContextData>(
        new ContextData(input_size_));
    data = context_map_[context_].get();
  }
  float decay = 0.9 / pow(0.0000001 * steps_ + 0.8, 0.8);
  decay *= 1.5 - ((1.0 * data->steps) / max_steps_);
  float update = decay * learning_rate_ * (bit - Sigmoid::Logistic(p_));
  ++steps_;
  ++data->steps;
  if (data->steps > max_steps_) {
    max_steps_ = data->steps;
  }
  data->weights += update * inputs_;
}

