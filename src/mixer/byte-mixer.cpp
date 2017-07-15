#include "byte-mixer.h"

ByteMixer::ByteMixer(unsigned int input_size, unsigned int num_cells,
    unsigned int num_layers, int horizon, float learning_rate,
    const unsigned int& bit_context) : byte_(bit_context),
    lstm_(input_size, 256, num_cells, num_layers, horizon, learning_rate) {}

void ByteMixer::SetInput(int index, float val) {
  lstm_.SetInput(index, val);
}

void ByteMixer::ByteUpdate() {
  probs_ = lstm_.Perceive(byte_);
  ByteModel::ByteUpdate();
}
