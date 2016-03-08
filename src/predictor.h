#ifndef PREDICTOR_H
#define PREDICTOR_H

#include "mixer/logistic.h"
#include "mixer/mixer-input.h"
#include "mixer/mixer.h"
#include "mixer/byte-mixer.h"
#include "sse.h"
#include "models/model.h"
#include "models/byte-model.h"
#include "manager.h"

#include <vector>
#include <memory>

class Predictor {
 public:
  Predictor();
  float Predict();
  void Perceive(int bit);

 private:
  void PrintStats();
  unsigned long long GetNumModels();
  unsigned long long GetNumNeurons();
  unsigned long long GetNumConnections();
  void Add(Model* model);
  void AddByteModel(ByteModel* model);
  void Add(int layer, Mixer* mixer);
  void AddPAQ8L();
  void AddPAQ8HP();
  void AddPPM();
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
  void AddSSE();
  void AddMixers();

  std::vector<std::unique_ptr<Model>> models_;
  std::vector<std::unique_ptr<ByteModel>> byte_models_;
  std::vector<std::unique_ptr<SSE>> sse_;
  std::vector<std::unique_ptr<MixerInput>> layers_;
  std::vector<std::vector<std::unique_ptr<Mixer>>> mixers_;
  std::vector<unsigned int> auxiliary_;
  Manager manager_;
  Logistic logistic_;
  std::unique_ptr<ByteMixer> byte_mixer_;
};

#endif
