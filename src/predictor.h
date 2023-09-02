#ifndef PREDICTOR_H
#define PREDICTOR_H

#include "mixer/sigmoid.h"
#include "mixer/mixer-input.h"
#include "mixer/mixer.h"
#include "mixer/byte-mixer.h"
#include "mixer/sse.h"
#include "models/model.h"
#include "models/byte-model.h"
#include "context-manager.h"

#include <vector>
#include <set>
#include <memory>

class Predictor {
 public:
  Predictor(const std::vector<bool>& vocab);
  float Predict();
  void Perceive(int bit);
  void Pretrain(int bit);

 private:
  unsigned long long GetNumModels();
  void AddModel(Model* model);
  void AddByteModel(ByteModel* model);
  void AddMixer(int layer, const unsigned long long& context,
      float learning_rate);
  void AddByteMixer(ByteMixer* byte_mixer);
  void AddAuxiliary();
  void AddPAQ8();
  void AddFXCM();
  void AddPPMD();
  void AddBracket();
  void AddWord();
  void AddDirect();
  void AddMatch();
  void AddDoubleIndirect();
  void AddMixers();

  std::vector<std::unique_ptr<Model>> models_;
  std::vector<std::unique_ptr<ByteModel>> byte_models_;
  SSE sse_;
  std::vector<std::unique_ptr<MixerInput>> layers_;
  std::vector<std::vector<std::unique_ptr<Mixer>>> mixers_;
  std::vector<unsigned int> auxiliary_;
  ContextManager manager_;
  Sigmoid sigmoid_;
  std::vector<std::unique_ptr<ByteMixer>> byte_mixers_;
  std::vector<bool> vocab_;
};

#endif
