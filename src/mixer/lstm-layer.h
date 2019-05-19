#ifndef LSTM_LAYER_H
#define LSTM_LAYER_H

#include <valarray>
#include <stdlib.h>
#include <math.h>

struct NeuronLayer {
  NeuronLayer(unsigned int input_size, unsigned int num_cells, int horizon,
    int offset) : error_(num_cells),
    weights_(std::valarray<float>(input_size), num_cells),
    state_(std::valarray<float>(num_cells), horizon),
    update_(std::valarray<float>(input_size), num_cells),
    m_(std::valarray<float>(input_size), num_cells),
    v_(std::valarray<float>(input_size), num_cells),
    transpose_(std::valarray<float>(num_cells), input_size - offset) {};

  std::valarray<float> error_;
  std::valarray<std::valarray<float>> weights_, state_, update_, m_, v_,
      transpose_;
};

class LstmLayer {
 public:
  LstmLayer(unsigned int input_size, unsigned int auxiliary_input_size,
      unsigned int output_size, unsigned int num_cells, int horizon,
      float gradient_clip, float learning_rate);
  void ForwardPass(const std::valarray<float>& input, int input_symbol,
      std::valarray<float>* hidden, int hidden_start);
  void BackwardPass(const std::valarray<float>& input, int epoch,
      int layer, int input_symbol, std::valarray<float>* hidden_error);
  static inline float Rand() {
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  }

 private:
  std::valarray<float> state_, state_error_, stored_error_, cache_, cache2_;
  std::valarray<std::valarray<float>> tanh_state_, input_gate_state_,
      last_state_;
  float gradient_clip_, learning_rate_;
  unsigned int num_cells_, epoch_, horizon_, input_size_, output_size_;
  unsigned long long update_steps_ = 0;
  NeuronLayer forget_gate_, input_node_, output_gate_;

  void ClipGradients(std::valarray<float>* arr);
  void BackwardPass(NeuronLayer& neurons, const std::valarray<float>&input,
      int epoch, int layer, int input_symbol,
      std::valarray<float>* hidden_error);
};

#endif
