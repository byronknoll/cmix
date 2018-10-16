#ifndef MIXER_H
#define MIXER_H

#include <vector>
#include <valarray>
#include <unordered_map>
#include <memory>

struct ContextData {
  ContextData(unsigned long long input_size,
      unsigned long long extra_input_size) : steps(0), weights(input_size),
      extra_weights(extra_input_size) {};
  unsigned long long steps;
  std::valarray<float> weights, extra_weights;
};

class Mixer {
 public:
  Mixer(const std::valarray<float>& inputs,
      const std::vector<float>& extra_inputs, const unsigned long long& context,
      float learning_rate, unsigned int extra_input_size);
  float Mix();
  void Perceive(int bit);

 private:
  ContextData* GetContextData();
  const std::valarray<float>& inputs_;
  const std::vector<float>& extra_inputs_vec_;
  std::valarray<float> extra_inputs_;
  float p_, learning_rate_;
  const unsigned long long& context_;
  unsigned long long max_steps_, steps_;
  std::unordered_map<unsigned int, std::unique_ptr<ContextData>> context_map_;
};

#endif
