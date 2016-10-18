#ifndef LAYER_H
#define LAYER_H

#include <valarray>
#include <stdlib.h>
#include <math.h>

class Layer {
 public:
  Layer(unsigned int input_size, unsigned int auxiliary_input_size,
      unsigned int num_cells, int horizon, float learning_rate);
  const std::valarray<float>& ForwardPass(const std::valarray<float>& input);
  const std::valarray<float>& BackwardPass(const std::valarray<float>& input,
      const std::valarray<float>& hidden_error, int epoch);
  static inline float Rand() {
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  }
  static inline float Logistic(float val) { return 1 / (1 + exp(-val)); }

 private:
  std::valarray<float> state_, hidden_, hidden_error_, output_gate_error_,
      state_error_, input_node_error_, input_gate_error_, forget_gate_error_,
      stored_error_;
  std::valarray<std::valarray<float>> tanh_state_, output_gate_state_,
      input_node_state_, input_gate_state_, forget_gate_state_, last_state_,
      forget_gate_, input_node_, input_gate_, output_gate_, forget_gate_update_,
      input_node_update_, input_gate_update_, output_gate_update_;
  float learning_rate_;
  unsigned int num_cells_, epoch_, horizon_, input_size_;
};

#endif
