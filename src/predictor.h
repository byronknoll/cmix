#ifndef PREDICTOR_H
#define PREDICTOR_H

#include "mixer/logistic.h"
#include "mixer/mixer-input.h"
#include "mixer/mixer.h"
#include "mixer/byte-mixer.h"
#include "mixer/sse.h"
#include "models/model.h"
#include "models/byte-model.h"
#include "context-manager.h"

#include <vector>
#include <memory>

class Predictor {
 public:
  Predictor(const std::vector<bool>& vocab);
  float Predict();
  void Perceive(int bit);
  void Pretrain(int bit);

 private:
  void PrintStats();
  unsigned long long GetNumModels();
  unsigned long long GetNumNeurons();
  unsigned long long GetNumConnections();
  void AddModel(Model* model);
  void AddByteModel(ByteModel* model);
  void AddMixer(int layer, Mixer* mixer);
  void AddByteMixer(ByteMixer* byte_mixer);
  void AddPAQ8L();
  void AddPAQ8HP();
  void AddPPM();
  void AddPPMD();
  void AddBracket();
  void AddDMC();
  void AddByteRun();
  void AddNonstationary();
  void AddEnglish();
  void AddSparse();
  void AddDirect();
  void AddRunMap();
  void AddMatch();
  void AddDoubleIndirect();
  void AddInterval();
  void AddMixers();

  std::vector<std::unique_ptr<Model>> models_;
  std::vector<std::unique_ptr<ByteModel>> byte_models_;
  SSE sse_;
  std::vector<std::unique_ptr<MixerInput>> layers_;
  std::vector<std::vector<std::unique_ptr<Mixer>>> mixers_;
  std::vector<unsigned int> auxiliary_;
  ContextManager manager_;
  Logistic logistic_;
  std::vector<std::unique_ptr<ByteMixer>> byte_mixers_;
  std::vector<bool> vocab_;
};

#endif
