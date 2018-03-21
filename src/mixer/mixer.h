#ifndef MIXER_H
#define MIXER_H

#include <vector>
#include <valarray>
#include <unordered_map>
#include <memory>

struct ContextData {
  ContextData(unsigned long long input_size) : steps(0),
      weights(0.0, input_size) {};
  unsigned long long steps;
  std::valarray<float> weights;
};

class Mixer {
 public:
  Mixer(const std::valarray<float>& inputs, const unsigned long long& context,
      float learning_rate, unsigned long long input_size);
  float Mix();
  void Perceive(int bit);

 private:
  const std::valarray<float>& inputs_;
  float p_, learning_rate_;
  const unsigned long long& context_;
  unsigned long long max_steps_, steps_, input_size_;
  std::unordered_map<unsigned int, std::unique_ptr<ContextData>> context_map_;
};

#endif
