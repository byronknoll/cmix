#ifndef LSTM_COMPRESS_H
#define LSTM_COMPRESS_H

#include <valarray>
#include <vector>
#include <memory>

#include "layer.h"

class Lstm {
 public:
  Lstm(unsigned int input_size, unsigned int num_cells, unsigned int num_layers,
      int horizon, float learning_rate);
  std::valarray<float>& Perceive(unsigned char input);
  std::valarray<float>& Predict(unsigned char input);
  void SetInput(int index, float val);

 private:
  std::vector<std::unique_ptr<Layer>> layers_;
  std::vector<unsigned char> input_history_;
  std::valarray<float> hidden_, hidden_error_;
  std::valarray<std::valarray<std::valarray<float>>> layer_input_,
      output_layer_;
  std::valarray<std::valarray<float>> output_;
  float learning_rate_;
  unsigned int num_cells_, epoch_, horizon_, input_size_;
};

#endif

