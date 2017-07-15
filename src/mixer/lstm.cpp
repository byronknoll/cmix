#include "lstm.h"

#include <numeric>

Lstm::Lstm(unsigned int input_size, unsigned int output_size, unsigned int
    num_cells, unsigned int num_layers, int horizon, float learning_rate) :
    input_history_(horizon), hidden_(num_cells * num_layers + 1),
    hidden_error_(num_cells),
    layer_input_(std::valarray<std::valarray<float>>(std::valarray<float>
    (input_size + output_size + 1 + num_cells * 2), num_layers), horizon),
    output_layer_(std::valarray<std::valarray<float>>(std::valarray<float>
    (num_cells * num_layers + 1), output_size), horizon),
    output_(std::valarray<float>(1.0 / output_size, output_size), horizon),
    learning_rate_(learning_rate), num_cells_(num_cells), epoch_(0),
    horizon_(horizon), input_size_(input_size), output_size_(output_size) {
  hidden_[hidden_.size() - 1] = 1;
  for (int epoch = 0; epoch < horizon; ++epoch) {
    layer_input_[epoch][0].resize(output_size + 1 + num_cells + input_size);
    for (unsigned int i = 0; i < num_layers; ++i) {
      layer_input_[epoch][i][layer_input_[epoch][i].size() - 1] = 1;
    }
  }
  for (unsigned int i = 0; i < num_layers; ++i) {
    layers_.push_back(std::unique_ptr<Layer>(new Layer(layer_input_[0][i].
        size(), input_size_, output_size_, num_cells, horizon, learning_rate)));
  }
}

void Lstm::SetInput(int index, float val) {
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    layer_input_[epoch_][i][output_size_ + index] = val;
  }
}

std::valarray<float>& Lstm::Perceive(unsigned int input) {
  int last_epoch = epoch_ - 1;
  if (last_epoch == -1) last_epoch = horizon_ - 1;
  input_history_[last_epoch] = input;
  if (epoch_ == 0) {
    for (int epoch = horizon_ - 1; epoch >= 0; --epoch) {
      for (int layer = layers_.size() - 1; layer >= 0; --layer) {
        int offset = layer * num_cells_;
        for (unsigned int i = 0; i < output_size_; ++i) {
          float error = 0;
          if (i == input_history_[epoch]) error = (1 - output_[epoch][i]);
          else error = -output_[epoch][i];
          for (unsigned int j = 0; j < hidden_error_.size(); ++j) {
            hidden_error_[j] += output_layer_[epoch][i][j + offset] * error;
          }
        }
        hidden_error_ = layers_[layer]->BackwardPass(layer_input_[epoch][layer],
            hidden_error_, epoch);
      }
    }
  }

  output_layer_[epoch_] = output_layer_[last_epoch];
  for (unsigned int i = 0; i < output_size_; ++i) {
    float error = 0;
    if (i == input) error = (1 - output_[last_epoch][i]);
    else error = -output_[last_epoch][i];
    output_layer_[epoch_][i] += learning_rate_ * error * hidden_;
  }
  return Predict(input);
}

std::valarray<float>& Lstm::Predict(unsigned int input) {
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    std::fill_n(begin(layer_input_[epoch_][i]), output_size_, 0);
    layer_input_[epoch_][i][input] = 1;
    auto start = begin(hidden_) + i * num_cells_;
    std::copy(start, start + num_cells_, begin(layer_input_[epoch_][i]) +
        output_size_ + input_size_);
    const auto& hidden = layers_[i]->ForwardPass(layer_input_[epoch_][i]);
    std::copy(begin(hidden), end(hidden), start);
    if (i < layers_.size() - 1) {
      start = begin(layer_input_[epoch_][i + 1]) + output_size_ + num_cells_ +
          input_size_;
      std::copy(begin(hidden), end(hidden), start);
    }
  }
  for (unsigned int i = 0; i < output_size_; ++i) {
    output_[epoch_][i] = exp(std::inner_product(&hidden_[0],
        &hidden_[hidden_.size()], &output_layer_[epoch_][i][0], 0.0));
  }
  double sum = 0;
  for (unsigned int i = 0; i < output_size_; ++i) {
    sum += output_[epoch_][i];
  }
  output_[epoch_] /= sum;
  int epoch = epoch_;
  ++epoch_;
  if (epoch_ == horizon_) epoch_ = 0;
  return output_[epoch];
}
