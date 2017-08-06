#include "byte-mixer.h"

ByteMixer::ByteMixer(unsigned int num_models, unsigned int num_cells,
    unsigned int num_layers, int horizon, float learning_rate,
    const unsigned int& bit_context, const std::vector<bool>& vocab,
    unsigned int vocab_size) : ByteModel(vocab), byte_(bit_context), lstm_(
    vocab_size, vocab_size, num_cells, num_layers, horizon, learning_rate),
    inputs_(vocab_size, 0), byte_map_(256, 0),
    outputs_(1.0 / vocab_size, vocab_size), num_models_(num_models),
    vocab_size_(vocab_size), offset_(0) {
  for (int i = 0; i < 256; ++i) {
    byte_map_[i] = offset_;
    if (vocab_[i]) ++offset_;
  }
  offset_ = 0;
}

void ByteMixer::SetInput(int index, float val) {
  if (!vocab_[index]) return;
  inputs_[offset_] += val;
  ++offset_;
  if (offset_ == vocab_size_) offset_ = 0;
}

void ByteMixer::ByteUpdate() {
  for (unsigned int i = 0; i < vocab_size_; ++i) {
    lstm_.SetInput(i, 2*inputs_[i]/num_models_);
    inputs_[i] = 0;
  }
  outputs_ = lstm_.Perceive(byte_map_[byte_]);
  offset_ = 0;
  for (int i = 0; i < 256; ++i) {
    if (vocab_[i]) {
      probs_[i] = outputs_[offset_];
      ++offset_;
    } else {
      probs_[i] = 0;
    }
  }
  offset_ = 0;
  ByteModel::ByteUpdate();
}
