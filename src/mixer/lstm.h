#ifndef LSTM_COMPRESS_H
#define LSTM_COMPRESS_H

#include <valarray>
#include <vector>
#include <memory>
#include <string>

#include "lstm-layer.h"

#include "../ds/emhash_set.hpp"
#include "../ds/SmallVector.h"


class Lstm {
 public:
  Lstm(unsigned int input_size, unsigned int output_size, unsigned int
      num_cells, unsigned int num_layers, int horizon, float learning_rate,
      float gradient_clip);
  ~Lstm();
  std::valarray<float>& Perceive(unsigned int input);
  std::valarray<float>& Predict(unsigned int input);
  void SetInput(const std::valarray<float>& input);
  //void SaveToDisk(const std::string& path);
  //void LoadFromDisk(const std::string& path);

 private:
  llvm::SmallVector<LstmLayer, 1> layers_;
  std::vector<uint8_t> input_history_; // horizon
  std::valarray<float> hidden_, hidden_error_;
  std::valarray<std::valarray<std::valarray<float>>> layer_input_;
#ifdef P1_MODE
  std::valarray<std::valarray<float>> output_w_;
#if P1_MODE == 2
  std::valarray<std::valarray<float>> work_;
#endif
#else
  std::valarray<std::valarray<std::valarray<float>>> output_layer_;
#endif
  std::valarray<std::valarray<float>> output_;
  float learning_rate_;
#ifdef OLR_RAMP
  // position-scheduled output LR (code-dive 2026-07-11): switches from the early
  // value (OUTPUT_LR) to OLR_LATE once olr_bytes_ >= OLR_RAMP. ADDS MEMBERS =>
  // the define MUST arrive via CFLAGS_DEFINES (all three flag groups; ODR rule).
  unsigned long long olr_bytes_ = 0;
  float olr_cur_ = 0;
#endif
  unsigned int num_cells_, epoch_, horizon_, input_size_, output_size_;
  int last_input_ = -1;
};
#include "lstm.hpp"
#endif


