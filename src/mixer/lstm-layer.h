#ifndef LSTM_LAYER_H
#define LSTM_LAYER_H

#include <valarray>
#include <vector>
#include <stdlib.h>
#include <math.h>

#ifdef LSTM_QUANT
#include "lstm-quant.h"
#endif

struct NeuronLayer {
  NeuronLayer(unsigned int input_size, unsigned int num_cells, int horizon,
    int offset) : error_(num_cells), ivar_(horizon), gamma_(1.0, num_cells),
    gamma_u_(num_cells), gamma_m_(num_cells), gamma_v_(num_cells),
    beta_(num_cells), beta_u_(num_cells), beta_m_(num_cells),
    beta_v_(num_cells), weights_(std::valarray<float>(input_size), num_cells),
    state_(std::valarray<float>(num_cells), horizon),
    update_(std::valarray<float>(input_size), num_cells),
    m_(std::valarray<float>(input_size), num_cells),
    v_(std::valarray<float>(input_size), num_cells),
    transpose_(std::valarray<float>(num_cells), input_size - offset),
    norm_(std::valarray<float>(num_cells), horizon) {
#ifdef P3_BATCH
    error_hist_.resize(horizon);
    for (int e = 0; e < horizon; ++e) error_hist_[e].resize(num_cells);
#endif
  };

  std::valarray<float> error_, ivar_, gamma_, gamma_u_, gamma_m_, gamma_v_,
      beta_, beta_u_, beta_m_, beta_v_;
  std::valarray<std::valarray<float>> weights_, state_, update_, m_, v_,
      transpose_, norm_;
#ifdef P3_BATCH
  std::valarray<std::valarray<float>> error_hist_;
#endif
#ifdef LSTM_QUANT
  // Quantized shadow of the matvec block of weights_ (columns
  // [output_size, row_size-1); the embed columns and the trailing bias
  // column stay float). Rebuilt lazily after each Adam step (once per
  // horizon cycle). Sized by the owning LstmLayer's ctor.
  std::vector<LstmQuantWeight> qweights_;  // num_cells x qstride_, row-major
  std::vector<float> qscale_;              // per-row dequantization scale
  std::vector<int32_t> qrowsum_;           // per-row sum(q) (VNNI correction)
#endif
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
  std::vector<std::valarray<std::valarray<float>>*> Weights();
#ifdef LSTM_GRADCHECK
  // Test-only harness hook: expose the gate NeuronLayers so a
  // finite-difference gradcheck can read update_/gamma_u_/beta_u_ and
  // perturb weights_/gamma_/beta_. Methods only, never in a make build.
  std::vector<NeuronLayer*> Gates() {
    return {&forget_gate_, &input_node_, &output_gate_};
  }
#endif

 private:
  std::valarray<float> state_, state_error_, stored_error_;
  std::valarray<std::valarray<float>> tanh_state_, input_gate_state_,
      last_state_;
  float gradient_clip_, learning_rate_;
  unsigned int num_cells_, epoch_, horizon_, input_size_, output_size_;
  unsigned long long update_steps_ = 0;
  NeuronLayer forget_gate_, input_node_, output_gate_;
#ifdef P3_BATCH
  std::vector<const std::valarray<float>*> input_hist_;
  std::vector<int> sym_hist_;
#endif
#ifdef LSTM_QUANT
  bool qdirty_ = true;         // set after Adam updates weights_ (epoch 0)
  unsigned int qcols_ = 0;     // quantized matvec columns (excl. bias col)
  unsigned int qstride_ = 0;   // qcols_ padded to LSTM_QUANT_PAD
  float qact_scale_ = 1.0f;    // per-step activation dequantization scale
  std::vector<LstmQuantAct> qinput_;   // quantized activations (padded)
#ifdef LSTM_QUANT_KERNEL_VNNI
  std::vector<uint8_t> qinput_biased_;  // qinput_ + 128 (vpdpbusd operand)
#endif

  void QuantRequantize();
  void QuantRequantizeGate(NeuronLayer& neurons);
  void QuantizeInput(const std::valarray<float>& input);
  LstmQuantAcc QuantDotRow(const NeuronLayer& neurons, unsigned int i) const;
#endif

  void ClipGradients(std::valarray<float>* arr);
  void ForwardPass(NeuronLayer& neurons, const std::valarray<float>& input,
      int input_symbol);
  void BackwardPass(NeuronLayer& neurons, const std::valarray<float>&input,
      int epoch, int layer, int input_symbol,
      std::valarray<float>* hidden_error);
};
#include "lstm-layer.hpp"

#endif

