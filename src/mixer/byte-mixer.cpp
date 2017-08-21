#include "byte-mixer.h"

ByteMixer::ByteMixer(unsigned int num_models, unsigned int num_cells,
    unsigned int num_layers, int horizon, float learning_rate,
    const unsigned int& bit_context) : byte_(bit_context), lstm_(
    256, 256, num_cells, num_layers, horizon, learning_rate),
    inputs_(0.0, 256), num_models_(num_models) {}

void ByteMixer::SetInput(int index, float val) {
  inputs_[index] += val;
}

void ByteMixer::ByteUpdate() {
  for (int i = 0; i < 256; ++i) {
    lstm_.SetInput(i, 2*inputs_[i] / num_models_);
  }
  inputs_ = 0;
  probs_ = lstm_.Perceive(byte_);
  ByteModel::ByteUpdate();
}
