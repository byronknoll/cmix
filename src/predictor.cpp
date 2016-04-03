#include "predictor.h"
#include "models/direct.h"
#include "models/direct-hash.h"
#include "models/indirect.h"
#include "models/byte-run.h"
#include "models/match.h"
#include "models/dmc.h"
#include "models/ppm.h"
#include "models/bracket.h"
#include "models/paq8l.h"
#include "models/paq8hp.h"
#include "contexts/context-hash.h"
#include "contexts/bracket-context.h"
#include "contexts/sparse.h"
#include "contexts/indirect-hash.h"
#include "contexts/bit-context.h"
#include "models/facade.h"

#include <vector>
#include <stdio.h>

Predictor::Predictor() : manager_(), logistic_(10000, 1000) {
  AddBracket();
  AddPAQ8HP();
  AddPAQ8L();
  AddPPM();
  AddDMC();
  AddByteRun();
  AddNonstationary();
  AddEnglish();
  AddSparse();
  AddDirect();
  AddRunMap();
  AddMatch();
  AddDoubleIndirect();

  AddMixers();
  AddSSE();

  // PrintStats();
}

void Predictor::PrintStats() {
  printf("Number of models: %llu\n", GetNumModels());
  printf("Number of neurons: %llu\n", GetNumNeurons());
  printf("Number of connections: %llu\n", GetNumConnections());
}

unsigned long long Predictor::GetNumModels() {
  return models_.size() + byte_models_.size() + 1;
}

unsigned long long Predictor::GetNumNeurons() {
  unsigned long long neurons = GetNumModels();
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    for (const auto& mixer : mixers_[i]) {
      neurons += mixer->GetNumNeurons();
    }
  }
  return neurons;
}

unsigned long long Predictor::GetNumConnections() {
  unsigned long long connections = 0;
  for (unsigned int i = 0; i < layers_.size(); ++i) {
    for (const auto& mixer : mixers_[i]) {
      connections += mixer->GetNumConnections();
    }
  }
  return connections;
}

void Predictor::Add(Model* model) {
  models_.push_back(std::unique_ptr<Model>(model));
}

void Predictor::AddByteModel(ByteModel* model) {
  byte_models_.push_back(std::unique_ptr<ByteModel>(model));
}

void Predictor::Add(int layer, Mixer* mixer) {
  mixers_[layer].push_back(std::unique_ptr<Mixer>(mixer));
}

void Predictor::AddPAQ8HP() {
  auxiliary_.push_back(models_.size());
  PAQ8HP* paq = new PAQ8HP(11);
  Add(paq);
  const std::vector<float>& predictions = paq->ModelPredictions();
  for (unsigned int i = 0; i < predictions.size(); ++i) {
    Add(new Facade(predictions[i]));
  }
}

void Predictor::AddPAQ8L() {
  auxiliary_.push_back(models_.size());
  PAQ8L* paq = new PAQ8L(10);
  Add(paq);
  const std::vector<float>& predictions = paq->ModelPredictions();
  for (unsigned int i = 0; i < predictions.size(); ++i) {
    Add(new Facade(predictions[i]));
  }
}

void Predictor::AddBracket() {
  Add(new Bracket(manager_.bit_context_, 200, 10, 100000));
  const Context& context = manager_.AddContext(std::unique_ptr<Context>(
      new BracketContext(manager_.bit_context_, 256, 15)));
  Add(new Direct(context.context_, manager_.bit_context_, 30, 0,
      context.size_));
  Add(new Indirect(manager_.nonstationary_, context.context_,
      manager_.bit_context_, 300, context.size_));
}

void Predictor::AddPPM() {
  AddByteModel(new PPM(7, manager_.bit_context_, 10000, 80000000));
  AddByteModel(new PPM(5, manager_.bit_context_, 10000, 20000000));
  AddByteModel(new PPM(3, manager_.bit_context_, 10000, 5000000));
}

void Predictor::AddDMC() {
  Add(new DMC(0.02, 70000000));
}

void Predictor::AddByteRun() {
  unsigned long long max_size = 10000000;
  float delta = 200;
  std::vector<std::vector<unsigned int>> model_params = {{0, 8}, {1, 5}, {1, 8},
      {2, 8}};

  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    Add(new ByteRun(context.context_, manager_.bit_context_, delta,
        std::min(max_size, context.size_)));
  }
}

void Predictor::AddNonstationary() {
  unsigned long long max_size = 1000000;
  float delta = 500;
  std::vector<std::vector<unsigned int>> model_params = {{0, 8}, {2, 8}, {4, 7},
      {8, 3}, {12, 1}, {16, 1}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    Add(new Indirect(manager_.nonstationary_, context.context_,
        manager_.bit_context_, delta, std::min(max_size, context.size_)));
  }
}

void Predictor::AddEnglish() {
  float delta = 200;
  unsigned long long max_size = 200000;
  std::vector<std::vector<unsigned int>> model_params = {{0}, {0, 1}, {7, 2},
      {7}, {1}, {1, 2}, {1, 2, 3}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {3, 4},
      {1, 2, 4}, {1, 2, 3, 4}, {2, 3, 4}, {2}, {1, 2, 3, 4, 5},
      {1, 2, 3, 4, 5, 6}};
  for (const auto& params : model_params) {
    std::unique_ptr<Context> hash(new Sparse(manager_.words_, params));
    const Context& context = manager_.AddContext(std::move(hash));
    Add(new Indirect(manager_.nonstationary_, context.context_,
        manager_.bit_context_, delta, max_size));
  }

  std::vector<std::vector<unsigned int>> model_params2 = {{0}, {1}, {7},
      {1, 3}, {1, 2, 3}, {7, 2}};
  for (const auto& params : model_params2) {
    std::unique_ptr<Context> hash(new Sparse(manager_.words_, params));
    const Context& context = manager_.AddContext(std::move(hash));
    Add(new Match(manager_.history_, context.context_, manager_.bit_context_,
        200, 0.5, 10000000, &(manager_.longest_match_)));
    Add(new ByteRun(context.context_, manager_.bit_context_, 100, 10000000));
    if (params[0] == 1 && params.size() == 1) {
      Add(new Indirect(manager_.run_map_, context.context_,
          manager_.bit_context_, delta, 500000));
      Add(new DirectHash(context.context_, manager_.bit_context_, 30, 0,
          500000));
    }
  }
}

void Predictor::AddSparse() {
  float delta = 300;
  unsigned long long max_size = 256;
  std::vector<std::vector<unsigned int>> model_params = {{1}, {2}, {3}, {4},
      {5}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {1, 2},
      {1, 3}, {2, 3}, {2, 5}, {3, 4}, {3, 5}, {3, 7}};
  for (const auto& params : model_params) {
    if (params.size() > 1) max_size = 256 * 256;
    std::unique_ptr<Context> hash(new Sparse(manager_.recent_bytes_, params));
    const Context& context = manager_.AddContext(std::move(hash));
    Add(new Indirect(manager_.nonstationary_, context.context_,
        manager_.bit_context_, delta, max_size));
  }
  std::vector<std::vector<unsigned int>> model_params2 = {{1}, {0, 2}, {0, 4},
      {1, 2}, {2, 3}, {3, 4}, {3, 7}};
  for (const auto& params : model_params2) {
    std::unique_ptr<Context> hash(new Sparse(manager_.recent_bytes_, params));
    const Context& context = manager_.AddContext(std::move(hash));
    Add(new Match(manager_.history_, context.context_, manager_.bit_context_,
        200, 0.5, 10000000, &(manager_.longest_match_)));
    Add(new ByteRun(context.context_, manager_.bit_context_, 100, 10000000));
  }
}

void Predictor::AddDirect() {
  float delta = 0;
  int limit = 30;
  std::vector<std::vector<int>> model_params = {{0, 8}, {1, 8}, {2, 8}, {3, 8}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_,params[0], params[1])));
    if (params[0] < 3) {
      Add(new Direct(context.context_, manager_.bit_context_, limit, delta,
          context.size_));
    } else {
      Add(new DirectHash(context.context_, manager_.bit_context_, limit, delta,
          100000));
    }
  }
}

void Predictor::AddRunMap() {
  unsigned long long max_size = 5000;
  float delta = 200;
  std::vector<std::vector<unsigned int>> model_params = {{0, 8}, {1, 5}, {1, 7},
      {1, 8}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    Add(new Indirect(manager_.run_map_, context.context_, manager_.bit_context_,
        delta, std::min(max_size, context.size_)));
  }
}

void Predictor::AddMatch() {
  float delta = 0.5;
  int limit = 200;
  unsigned long long max_size = 20000000;
  std::vector<std::vector<int>> model_params = {{0, 8}, {1, 8}, {2, 8}, {7, 4},
      {11, 3}, {13, 2}, {15, 2}, {17, 2}, {20, 1}, {25, 1}};

  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_,params[0], params[1])));
    Add(new Match(manager_.history_, context.context_, manager_.bit_context_,
        limit, delta, std::min(max_size, context.size_),
        &(manager_.longest_match_)));
  }
}

void Predictor::AddDoubleIndirect() {
  unsigned long long max_size = 100000;
  float delta = 400;
  std::vector<std::vector<unsigned int>> model_params = {{1, 8, 1, 8},
      {2, 8, 1, 8}, {1, 8, 2, 8}, {2, 8, 2, 8}, {1, 8, 3, 8}, {3, 8, 1, 8},
      {4, 6, 4, 8}, {5, 5, 5, 5}, {1, 8, 4, 8}, {1, 8, 5, 6}, {6, 4, 6, 4}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new IndirectHash(manager_.bit_context_, params[0], params[1],
        params[2], params[3])));
    Add(new Indirect(manager_.nonstationary_, context.context_,
        manager_.bit_context_, delta, std::min(max_size, context.size_)));
  }
}

void Predictor::AddSSE() {
  std::vector<std::vector<unsigned int>> model_params = {{0, 8, 2000, 1000},
      {1, 8, 500, 100}, {2, 4, 500, 100}, {3, 3, 200, 100}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    sse_.push_back(std::unique_ptr<SSE>(new SSE(logistic_, context.context_,
        manager_.bit_context_, params[2], params[3], context.size_)));
  }
}

void Predictor::AddMixers() {
  byte_mixer_.reset(new ByteMixer(logistic_, 256 * byte_models_.size(), 128,
      manager_.bit_context_, 0.3));

  for (int i = 0; i < 3; ++i) {
    layers_.push_back(std::unique_ptr<MixerInput>(new MixerInput(logistic_,
        1.0e-4)));
    mixers_.push_back(std::vector<std::unique_ptr<Mixer>>());
  }

  unsigned long long input_size = models_.size() + byte_models_.size() + 1;
  layers_[0]->SetNumModels(input_size);
  std::vector<std::vector<double>> model_params = {{0, 8, 0.005},
      {0, 8, 0.0005}, {1, 8, 0.005}, {1, 8, 0.0005}, {2, 4, 0.005},
      {3, 2, 0.002}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    const BitContext& bit_context = manager_.AddBitContext(std::unique_ptr
        <BitContext>(new BitContext(manager_.long_bit_context_,
        context.context_, context.size_)));
    Add(0, new Mixer(layers_[0]->inputs_, logistic_, bit_context.context_,
        params[2], bit_context.size_, input_size));
  }

  model_params = {{0, 0.001}, {2, 0.002}, {3, 0.005}};
  for (const auto& params : model_params) {
    Add(0, new Mixer(layers_[0]->inputs_, logistic_,
        manager_.recent_bytes_[params[0]], params[1], 256, input_size));
  }
  Add(0, new Mixer(layers_[0]->inputs_, logistic_, manager_.zero_context_,
      0.00005, 1, input_size));
  Add(0, new Mixer(layers_[0]->inputs_, logistic_, manager_.line_break_,
      0.0007, 100, input_size));
  Add(0, new Mixer(layers_[0]->inputs_, logistic_, manager_.longest_match_,
      0.0005, 8, input_size));

  const BitContext& bit_context1 = manager_.AddBitContext(std::unique_ptr
      <BitContext>(new BitContext(manager_.long_bit_context_,
      manager_.recent_bytes_[1], 256)));
  Add(0, new Mixer(layers_[0]->inputs_, logistic_, bit_context1.context_,
      0.005, bit_context1.size_, input_size));

  const BitContext& bit_context2 = manager_.AddBitContext(std::unique_ptr
      <BitContext>(new BitContext(manager_.recent_bytes_[1],
      manager_.recent_bytes_[0], 256)));
  Add(0, new Mixer(layers_[0]->inputs_, logistic_, bit_context2.context_, 0.005,
      bit_context2.size_, input_size));

  const BitContext& bit_context3 = manager_.AddBitContext(std::unique_ptr
      <BitContext>(new BitContext(manager_.recent_bytes_[2],
      manager_.recent_bytes_[1], 256)));
  Add(0, new Mixer(layers_[0]->inputs_, logistic_, bit_context3.context_, 0.003,
      bit_context3.size_, input_size));

  input_size = mixers_[0].size() + auxiliary_.size();
  layers_[1]->SetNumModels(input_size);

  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.zero_context_,
      0.005, 1, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.zero_context_,
      0.0005, 1, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.long_bit_context_,
      0.005, 256, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.long_bit_context_,
      0.0005, 256, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.long_bit_context_,
      0.00001, 256, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.recent_bytes_[0],
      0.005, 256, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.recent_bytes_[1],
      0.005, 256, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.recent_bytes_[2],
      0.005, 256, input_size));
  Add(1, new Mixer(layers_[1]->inputs_, logistic_, manager_.longest_match_,
      0.0005, 8, input_size));

  input_size = mixers_[1].size() + auxiliary_.size();
  layers_[2]->SetNumModels(input_size);
  Add(2, new Mixer(layers_[2]->inputs_, logistic_, manager_.zero_context_,
      0.0003, 1, input_size));
}

float Predictor::Predict() {
  // return models_[0]->Predict();
  for (unsigned int i = 0; i < models_.size(); ++i) {
    float p = models_[i]->Predict();
    layers_[0]->SetInput(i, p);
  }
  for (unsigned int i = 0; i < byte_models_.size(); ++i) {
    float p = byte_models_[i]->Predict();
    layers_[0]->SetInput(models_.size() + i, p);
  }
  float byte_mixer_p = byte_mixer_->Predict();
  layers_[0]->SetInput(models_.size() + byte_models_.size(), byte_mixer_p);
  for (unsigned int layer = 1; layer <= 2; ++layer) {
    for (unsigned int i = 0; i < mixers_[layer - 1].size(); ++i) {
      layers_[layer]->SetInput(i, mixers_[layer - 1][i]->Mix());
    }
    for (unsigned int i = 0; i < auxiliary_.size(); ++i) {
      layers_[layer]->inputs_[mixers_[layer - 1].size() + i] =
          layers_[0]->inputs_[auxiliary_[i]];
    }
  }
  float mixer_output = mixers_[2][0]->Mix();
  float p = logistic_.Stretch(mixer_output);
  p = (mixer_output + 3 * sse_[0]->Process(p) + sse_[1]->Process(p) +
      sse_[2]->Process(p) + sse_[3]->Process(p)) / 7;
  return p;
}

void Predictor::Perceive(int bit) {
  for (const auto& model : models_) {
    model->Perceive(bit);
  }
  for (const auto& model : byte_models_) {
    model->Perceive(bit);
  }
  byte_mixer_->Perceive(bit);
  for (unsigned int i = 0; i < mixers_.size(); ++i) {
    for (const auto& mixer : mixers_[i]) {
      mixer->Perceive(bit);
    }
  }
  for (const auto& sse : sse_) {
    sse->Perceive(bit);
  }

  bool byte_update = false;
  if (manager_.bit_context_ >= 128) byte_update = true;

  manager_.Perceive(bit);
  if (byte_update) {
    for (const auto& model : models_) {
      model->ByteUpdate();
    }
    for (const auto& model : byte_models_) {
      model->ByteUpdate();
    }
    byte_mixer_->Train();
    for (unsigned int i = 0; i < byte_models_.size(); ++i) {
      const std::valarray<float>& p = byte_models_[i]->BytePredict();
      for (unsigned int j = 0; j < 256; ++j) {
        byte_mixer_->SetInput(i*256 + j, p[j]);
      }
    }
    byte_mixer_->ByteUpdate();
    manager_.bit_context_ = 1;
  }
}
