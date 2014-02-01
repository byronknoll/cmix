#ifndef MIXER_H
#define MIXER_H

#include "logistic.h"

#include <vector>
#include <array>

class Mixer {
 public:
  Mixer(const std::vector<float>& inputs, const Logistic& logistic,
      const unsigned long long& byte_context, const unsigned int& bit_context,
      float learning_rate, unsigned long long table_size);
  void SetNumModels(int num_models);
  float Mix();
  void Perceive(int bit);

 private:
  const std::vector<float>& inputs_;
  const Logistic& logistic_;
  float p_;
  float learning_rate_;
  const unsigned int& bit_context_;
  const unsigned long long& byte_context_;
  std::vector<std::array<std::vector<float>, 256>> weights_;
};

#endif
