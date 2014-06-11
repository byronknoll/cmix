#ifndef MIXER_H
#define MIXER_H

#include "logistic.h"

#include <vector>
#include <valarray>

class Mixer {
 public:
  Mixer(const std::valarray<float>& inputs, const Logistic& logistic,
      const unsigned long long& context, float learning_rate,
      unsigned long long context_size, unsigned long long input_size);
  float Mix();
  void Perceive(int bit);
  unsigned long long GetNumNeurons();

 private:
  const std::valarray<float>& inputs_;
  const Logistic& logistic_;
  float p_;
  float learning_rate_;
  const unsigned long long& context_;
  std::vector<std::valarray<float>> weights_;
};

#endif
