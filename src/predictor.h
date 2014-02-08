#ifndef PREDICTOR_H
#define PREDICTOR_H

#include "mixer/logistic.h"
#include "mixer/mixer-input.h"
#include "mixer/mixer.h"
#include "sse.h"
#include "models/model.h"
#include "manager.h"

#include <vector>
#include <memory>

class Predictor {
 public:
  Predictor(unsigned long long file_size);
  float Predict();
  void Perceive(int bit);

 private:
  void Add(Model* model);
  void Add(int layer, Mixer* mixer);
  void AddByteRun();
  void AddNonstationary();
  void AddEnglish();
  void AddSparse();
  void AddDirect();
  void AddRunMap();
  void AddMatch();
  void AddPic();
  void AddSSE();
  void AddMixers();

  std::vector<std::unique_ptr<Model>> models_;
  std::vector<std::unique_ptr<SSE>> sse_;
  std::vector<std::unique_ptr<MixerInput>> layers_;
  std::vector<std::vector<std::unique_ptr<Mixer>>> mixers_;
  Manager manager_;
  Logistic logistic_;
};

#endif
