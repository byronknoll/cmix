#ifndef LSTM_LAYER_H
#define LSTM_LAYER_H

#include <valarray>
#include <stdlib.h>
#include <math.h>

class LstmLayer {
 public:
  LstmLayer(unsigned int input_size, unsigned int auxiliary_input_size,
      unsigned int output_size, unsigned int num_cells, int horizon,
      float gradient_clip);
  void ForwardPass(const std::valarray<float>& input, int input_symbol,
      std::valarray<float>* hidden, int hidden_start);
  void BackwardPass(const std::valarray<float>& input, int epoch,
      int layer, int input_symbol, std::valarray<float>* hidden_error);
  static inline float Rand() {
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  }

 private:
  std::valarray<float> state_, output_gate_error_, state_error_,
      input_node_error_, forget_gate_error_, stored_error_;
  std::valarray<std::valarray<float>> tanh_state_, output_gate_state_,
      input_node_state_, input_gate_state_, forget_gate_state_, last_state_,
      forget_gate_, input_node_, output_gate_, forget_gate_update_,
      input_node_update_, output_gate_update_, forget_gate_m_, input_node_m_,
      output_gate_m_, forget_gate_v_, input_node_v_, output_gate_v_;
  float gradient_clip_;
  unsigned int num_cells_, epoch_, horizon_, input_size_, output_size_;
  unsigned long long update_steps_ = 0;

  void ClipGradients(std::valarray<float>* arr);
};

#endif
