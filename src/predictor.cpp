#include "predictor.h"
#include "models/direct.h"
#include "models/indirect.h"
#include "models/match.h"
#include "models/bit-buffer.h"
#include "contexts/context-hash.h"
#include "contexts/rubin-karp.h"
#include "contexts/sparse.h"

#include <vector>

Predictor::Predictor(unsigned long long file_size) : manager_(file_size) {
  AddNonstationary();
  AddEnglish();
  AddSparse();
  AddDirect();
  AddRunMap();
  AddMatch();
  AddPic();
  AddSSE();
  AddMixers();
}

int Predictor::Discretize(float p) {
  return 1 + 65534 * p;
}

void Predictor::Add(Model* model) {
  models_.push_back(std::unique_ptr<Model>(model));
}

void Predictor::AddNonstationary() {
  unsigned long long max_size = 1000000;
  float delta = 500;
  std::vector<std::vector<unsigned int>> model_params = {{0, 8}, {1, 5},
      {1, 8}, {2, 6}, {2, 8}, {3, 7}, {4, 7}, {5, 4}, {6, 4}, {7, 3}, {8, 3},
      {9, 2}, {10,2}, {11, 2}, {12, 1}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    Add(new Indirect(manager_.nonstationary_, context.context_,
        manager_.bit_context_, delta, std::min(max_size, context.size_)));
  }

  max_size = 200000;
  const Context& context = manager_.AddContext(std::unique_ptr<Context>(
      new RubinKarp(manager_.bit_context_, 5)));
  Add(new Indirect(manager_.nonstationary_, context.context_,
      manager_.bit_context_, delta, std::min(max_size, context.size_ / 5)));
}

void Predictor::AddEnglish() {
  float delta = 200;
  unsigned long long max_size = 1500000;
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

  std::unique_ptr<Context> hash(new Sparse(manager_.words_,
      std::vector<unsigned int>(1, 1)));
  const Context& context = manager_.AddContext(std::move(hash));
  Add(new Indirect(manager_.nonstationary_, context.context_,
      manager_.bit_context_, delta, 1000));
  Add(new Indirect(manager_.run_map_, context.context_, manager_.bit_context_,
      delta, 1000000));
  Add(new Direct(context.context_, manager_.bit_context_, 30, 0, 1000000));
  Add(new Match(manager_.history_, context.context_, manager_.bit_context_, 200,
      0.5, 10000000));
}

void Predictor::AddSparse() {
  float delta = 300;
  unsigned long long max_size = 256;
  std::vector<std::vector<unsigned int>> model_params = {{1}, {2}, {3}, {4},
      {5}, {6}, {7}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {1, 2},
      {1, 3}, {2, 3}, {2, 5}, {3, 4}, {3, 5}, {3, 7}};
  for (const auto& params : model_params) {
    if (params.size() > 1) max_size = 256 * 256;
    std::unique_ptr<Context> hash(new Sparse(manager_.recent_bytes_, params));
    const Context& context = manager_.AddContext(std::move(hash));
    Add(new Indirect(manager_.nonstationary_, context.context_,
        manager_.bit_context_, delta, max_size));
  }
}

void Predictor::AddDirect() {
  float delta = 0;
  int limit = 30;
  std::vector<std::vector<int>> model_params = {{0, 8}, {1, 8}, {2, 8}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_,params[0], params[1])));
    Add(new Direct(context.context_, manager_.bit_context_, limit, delta,
        context.size_));
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
  std::vector<std::vector<int>> model_params = {{0, 8}, {1, 8}, {2, 8}, {3, 8},
    {5, 6}, {7, 5}, {9, 4}};

  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_,params[0], params[1])));
    Add(new Match(manager_.history_, context.context_, manager_.bit_context_,
        limit, delta, std::min(max_size, context.size_)));
  }
}

void Predictor::AddPic() {
  float delta = 0;
  int limit = 200;
  std::vector<int> params = {216 * 8 + 1, 216 * 8, 216 * 8 - 1, 216 * 8 - 2,
      216 * 8 - 3, 1, 2, 3, 4, 5, 6, 7};
  for (int size : params) {
    Add(new BitBuffer(size, delta, limit));
  }
}

void Predictor::AddSSE() {
  std::vector<std::vector<unsigned int>> model_params = {{0, 8, 2000, 1000},
      {1, 8, 500, 100}, {2, 4, 500, 100}, {3, 3, 200, 100}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    sse_.push_back(std::unique_ptr<SSE>(new SSE(context.context_,
        manager_.bit_context_, params[2], params[3], context.size_)));
  }
}

void Predictor::AddMixers() {
  for (int i = 0; i < 3; ++i) {
    layers_.push_back(std::unique_ptr<MixerInput>(new MixerInput(1.0e-4,
        10000, 1000)));
    mixers_.push_back(std::vector<std::unique_ptr<Mixer>>());
  }

  std::vector<std::vector<double>> model_params = {{0, 8, 0.005},
      {0, 8, 0.0005}, {1, 8, 0.005}, {1, 8, 0.0005}, {3, 3, 0.002}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    mixers_[0].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[0]),
        context.context_, manager_.bit_context_, params[2], context.size_)));
  }

  model_params = {{0, 0.001}, {1, 0.01}, {2, 0.002}, {3, 0.005}};
  const Context& context = manager_.AddContext(std::unique_ptr<Context>(
      new ContextHash(manager_.bit_context_, 0, 8)));
  for (const auto& params : model_params) {
    mixers_[0].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[0]),
        context.context_, manager_.recent_bytes2_[params[0]], params[1],
        context.size_)));
  }
  mixers_[0].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[0]),
      context.context_, manager_.zero_context_, 0.00005, context.size_)));
  mixers_[0].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[0]),
      manager_.recent_bytes_[1], manager_.bit_context_, 0.005, 256)));
  mixers_[0].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[0]),
      manager_.recent_bytes_[1], manager_.recent_bytes2_[0], 0.005, 256)));
  mixers_[0].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[0]),
      manager_.recent_bytes_[2], manager_.recent_bytes2_[1], 0.003, 256)));


  mixers_[1].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[1]),
      context.context_, manager_.zero_context_, 0.005, context.size_)));
  mixers_[1].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[1]),
      context.context_, manager_.zero_context_, 0.0005, context.size_)));
  mixers_[1].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[1]),
      context.context_, manager_.bit_context_, 0.005, context.size_)));
  mixers_[1].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[1]),
      context.context_, manager_.bit_context_, 0.0005, context.size_)));
  mixers_[1].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[1]),
      context.context_, manager_.recent_bytes2_[0], 0.005, context.size_)));
  mixers_[1].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[1]),
      context.context_, manager_.recent_bytes2_[1], 0.005, context.size_)));
  mixers_[1].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[1]),
      context.context_, manager_.recent_bytes2_[2], 0.005, context.size_)));

  mixers_[2].push_back(std::unique_ptr<Mixer>(new Mixer(*(layers_[2]),
      context.context_, manager_.zero_context_, 0.0003, context.size_)));

  layers_[0]->SetNumModels(models_.size());
  for (const auto& mixer : mixers_[0]) {
    mixer->SetNumModels(models_.size());
  }
  layers_[1]->SetNumModels(mixers_[0].size());
  for (const auto& mixer : mixers_[1]) {
    mixer->SetNumModels(mixers_[0].size());
  }
  layers_[2]->SetNumModels(mixers_[1].size());
  for (const auto& mixer : mixers_[2]) {
    mixer->SetNumModels(mixers_[1].size());
  }
}

int Predictor::Predict() {
  //return Discretize(models_[0]->Predict());
  for (unsigned int i = 0; i < models_.size(); ++i) {
    float p = models_[i]->Predict();
    layers_[0]->SetInput(i, p);
  }
  for (unsigned int i = 0; i < mixers_[0].size(); ++i) {
    layers_[1]->SetInput(i, mixers_[0][i]->Mix());
  }
  for (unsigned int i = 0; i < mixers_[1].size(); ++i) {
    layers_[2]->SetInput(i, mixers_[1][i]->Mix());
  }
  float p = mixers_[2][0]->Mix();
  p = (p + 2 * sse_[0]->Process(p) + sse_[1]->Process(p) + sse_[2]->Process(p) +
      sse_[3]->Process(p)) / 6;
  return Discretize(p);
}

void Predictor::Perceive(int bit) {
  for (const auto& model : models_) {
    model->Perceive(bit);
  }
  for (unsigned int i = 0; i < mixers_.size(); ++i) {
    for (const auto& mixer : mixers_[i]) {
      mixer->Perceive(bit);
    }
  }
  for (const auto& sse : sse_) {
    sse->Perceive(bit);
  }
  manager_.Perceive(bit);
  if (manager_.bit_context_ == 1) {
    for (const auto& model : models_) {
      model->ByteUpdate();
    }
  }
}
