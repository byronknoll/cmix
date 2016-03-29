#include "byte-mixer.h"

#include <stdlib.h>
#include <math.h>

ByteMixer::ByteMixer(const Logistic& logistic, int input_neurons,
    int hidden_neurons, const unsigned int& bit_context, float learning_rate) :
    byte_(bit_context), logistic_(logistic), learning_rate_(-learning_rate),
    steps_(1), outputs_(0.0, 256) {
  int output_neurons = 256;
  input_neurons += hidden_neurons;
  weights_.resize(2);
  srand(0xDEADBEEF);
  for (size_t layer = 0; layer < weights_.size(); ++layer) {
    if (layer == weights_.size() - 1) {
      weights_[layer].resize(output_neurons);
    } else {
      weights_[layer].resize(hidden_neurons);
    }
    for (size_t neuron = 0; neuron < weights_[layer].size(); ++neuron) {
      if (layer == 0) {
        weights_[layer][neuron].resize(input_neurons + 1);
      } else {
        weights_[layer][neuron].resize(weights_[layer - 1].size() + 1);
      }
      int num_in = weights_[layer][neuron].size() - 1;
      int num_out = weights_[layer].size() - 1;
      float low = -4 * (sqrt(6.0 / (num_in + num_out)));
      float range = -2 * low;
      for (size_t weight = 0; weight < weights_[layer][neuron].size();
          ++weight) {
        float w = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        weights_[layer][neuron][weight] = low + w * range;
      }
    }
  }

  states_.resize(weights_.size() + 1);
  for (size_t layer = 0; layer < states_.size(); ++layer) {
    if (layer == 0) {
      states_[layer].resize(input_neurons + 1);
    } else if (layer == states_.size() - 1) {
      states_[layer].resize(output_neurons);
    } else {
      states_[layer].resize(weights_[layer - 1].size() + 1);
    }
    states_[layer][states_[layer].size() - 1] = 1;
  }

  errors_.resize(weights_.size());
  for (size_t layer = 0; layer < errors_.size(); ++layer) {
    if (layer == errors_.size() - 1) {
      errors_[layer].resize(output_neurons);
    } else {
      errors_[layer].resize(weights_[layer].size() + 1);
    }
  }
}

void ByteMixer::SetInput(int index, float val) {
  states_[0][index] = val;
}

void ByteMixer::Train() {
  outputs_ = 0;
  outputs_[byte_] = 1;
  errors_[errors_.size() - 1] = states_[states_.size() - 1] - outputs_;
  for (int layer = errors_.size() - 2; layer >= 0; --layer) {
    int next = layer + 1;
    errors_[layer] = weights_[next][0] * errors_[next][0];
    for (size_t neuron = 1; neuron < weights_[next].size(); ++neuron) {
      errors_[layer] += weights_[next][neuron] * errors_[next][neuron];
    }
    errors_[layer] *= states_[next] * (1.0f - states_[next]);
  }

  float rate = (1.0 / sqrt(0.0000001 * steps_ + 0.8)) * learning_rate_;
  for (size_t layer = 0; layer < weights_.size(); ++layer) {
    for (size_t neuron = 0; neuron < weights_[layer].size(); ++neuron) {
      weights_[layer][neuron] += (rate * errors_[layer][neuron]) *
          states_[layer];
    }
  }
  ++steps_;
}

void ByteMixer::ByteUpdate() {
  for (size_t index = 0; index < states_[1].size() - 1; ++index) {
    int input_index = states_[0].size() - states_[1].size() + index;
    states_[0][input_index] = states_[1][index];
  }
  for (size_t layer = 0; layer < weights_.size(); ++layer) {
    int offset = layer + 1;
    for (size_t neuron = 0; neuron < weights_[layer].size(); ++neuron) {
      float x1 = 0, x2 = 0, x3 = 0, x4 = 0;
      size_t weight = 0;
      for (; weight < states_[layer].size() - 3; weight += 4) {
        x1 += states_[layer][weight] * weights_[layer][neuron][weight];
        x2 += states_[layer][weight+1] * weights_[layer][neuron][weight+1];
        x3 += states_[layer][weight+2] * weights_[layer][neuron][weight+2];
        x4 += states_[layer][weight+3] * weights_[layer][neuron][weight+3];
      }
      states_[offset][neuron] = x1 + x2 + x3 + x4;
      for (; weight < states_[layer].size(); ++weight) {
        states_[offset][neuron] += states_[layer][weight] *
            weights_[layer][neuron][weight];
      }
      states_[offset][neuron] = logistic_.Squash(states_[offset][neuron]);
    }
  }
  probs_ = states_[states_.size() - 1];
  ByteModel::ByteUpdate();
}
