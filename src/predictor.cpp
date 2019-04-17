#include "predictor.h"
#include "models/direct.h"
#include "models/direct-hash.h"
#include "models/indirect.h"
#include "models/byte-run.h"
#include "models/match.h"
#include "models/ppmd.h"
#include "models/bracket.h"
#include "models/paq8.h"
#include "models/paq8hp.h"
#include "contexts/context-hash.h"
#include "contexts/bracket-context.h"
#include "contexts/sparse.h"
#include "contexts/indirect-hash.h"
#include "contexts/interval.h"
#include "contexts/interval-hash.h"
#include "contexts/bit-context.h"
#include "contexts/combined-context.h"

#include <vector>
#include <stdlib.h>
#include <stdio.h>

Predictor::Predictor(const std::vector<bool>& vocab) : manager_(),
    sigmoid_(100001), vocab_(vocab) {
  srand(0xDEADBEEF);

  AddBracket();
  AddPAQ8HP();
  AddPAQ8();
  AddPPMD();
  AddWord();
  AddDirect();
  AddMatch();
  AddDoubleIndirect();
  AddMixers();
}

unsigned long long Predictor::GetNumModels() {
  unsigned long long num = 0;
  for (unsigned int i = 0; i < models_.size(); ++i) {
    num += models_[i]->NumOutputs();
  }
  for (unsigned int i = 0; i < byte_models_.size(); ++i) {
    num += byte_models_[i]->NumOutputs();
  }
  for (unsigned int i = 0; i < byte_mixers_.size(); ++i) {
    num += byte_mixers_[i]->NumOutputs();
  }
  return num;
}

void Predictor::AddModel(Model* model) {
  models_.push_back(std::unique_ptr<Model>(model));
}

void Predictor::AddByteModel(ByteModel* model) {
  byte_models_.push_back(std::unique_ptr<ByteModel>(model));
}

void Predictor::AddMixer(int layer, const unsigned long long& context,
    float learning_rate) {
  mixers_[layer].push_back(std::unique_ptr<Mixer>(new Mixer(
      layers_[layer]->Inputs(), layers_[layer]->ExtraInputs(), context,
      learning_rate, mixers_[layer].size())));
}

void Predictor::AddByteMixer(ByteMixer* byte_mixer) {
  byte_mixers_.push_back(std::unique_ptr<ByteMixer>(byte_mixer));
}

void Predictor::AddAuxiliary() {
  unsigned int index = GetNumModels() - 1;
  auxiliary_.push_back(index);
}

void Predictor::AddPAQ8HP() {
  PAQ8HP* paq = new PAQ8HP(11);
  AddModel(paq);
  AddAuxiliary();
}

void Predictor::AddPAQ8() {
  PAQ8* paq = new PAQ8(11);
  AddModel(paq);
  AddAuxiliary();
}

void Predictor::AddBracket() {
  AddModel(new Bracket(manager_.bit_context_, 200, 10, 100000, vocab_));
  const Context& context = manager_.AddContext(std::unique_ptr<Context>(
      new BracketContext(manager_.bit_context_, 256, 15)));
  AddModel(new Direct(context.GetContext(), manager_.bit_context_, 30, 0,
      context.Size()));
  AddModel(new Indirect(manager_.nonstationary_, context.GetContext(),
      manager_.bit_context_, 300, manager_.shared_map_));
}

void Predictor::AddPPMD() {
  AddByteModel(new PPMD::PPMD(6, 1200, manager_.bit_context_, vocab_));
  AddByteModel(new PPMD::PPMD(16, 1200, manager_.bit_context_, vocab_));
}

void Predictor::AddWord() {
  float delta = 200;
  std::vector<std::vector<unsigned int>> model_params = {{0}, {0, 1}, {7, 2},
      {7}, {1}, {1, 2}, {1, 2, 3}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {3, 4},
      {1, 2, 4}, {1, 2, 3, 4}, {2, 3, 4}, {2}, {1, 2, 3, 4, 5},
      {1, 2, 3, 4, 5, 6}};
  for (const auto& params : model_params) {
    std::unique_ptr<Context> hash(new Sparse(manager_.words_, params));
    const Context& context = manager_.AddContext(std::move(hash));
    AddModel(new Indirect(manager_.nonstationary_, context.GetContext(),
        manager_.bit_context_, delta, manager_.shared_map_));
  }

  std::vector<std::vector<unsigned int>> model_params2 = {{0}, {1}, {7},
      {1, 3}, {1, 2, 3}, {7, 2}};
  for (const auto& params : model_params2) {
    std::unique_ptr<Context> hash(new Sparse(manager_.words_, params));
    const Context& context = manager_.AddContext(std::move(hash));
    AddModel(new Match(manager_.history_, context.GetContext(),
        manager_.bit_context_, 200, 0.5, 10000000, &(manager_.longest_match_)));
    AddModel(new ByteRun(context.GetContext(), manager_.bit_context_, 100,
        10000000));
    if (params[0] == 1 && params.size() == 1) {
      AddModel(new Indirect(manager_.run_map_, context.GetContext(),
          manager_.bit_context_, delta, manager_.shared_map_));
      AddModel(new DirectHash(context.GetContext(), manager_.bit_context_, 30,
          0, 500000));
    }
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
      AddModel(new Direct(context.GetContext(), manager_.bit_context_, limit,
          delta, context.Size()));
    } else {
      AddModel(new DirectHash(context.GetContext(), manager_.bit_context_,
          limit, delta, 100000));
    }
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
    AddModel(new Match(manager_.history_, context.GetContext(),
        manager_.bit_context_, limit, delta, std::min(max_size, context.Size()),
        &(manager_.longest_match_)));
  }
}

void Predictor::AddDoubleIndirect() {
  float delta = 400;
  std::vector<std::vector<unsigned int>> model_params = {{1, 8, 1, 8},
      {2, 8, 1, 8}, {1, 8, 2, 8}, {2, 8, 2, 8}, {1, 8, 3, 8}, {3, 8, 1, 8},
      {4, 6, 4, 8}, {5, 5, 5, 5}, {1, 8, 4, 8}, {1, 8, 5, 6}, {6, 4, 6, 4}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new IndirectHash(manager_.bit_context_, params[0], params[1],
        params[2], params[3])));
    AddModel(new Indirect(manager_.nonstationary_, context.GetContext(),
        manager_.bit_context_, delta, manager_.shared_map_));
  }
}

void Predictor::AddMixers() {
  unsigned int vocab_size = 0;
  for (unsigned int i = 0; i < vocab_.size(); ++i) {
    if (vocab_[i]) ++vocab_size;
  }
  AddByteMixer(new ByteMixer(byte_models_.size(), 200, 2, 40, 0.035, 10,
      manager_.bit_context_, vocab_, vocab_size));
  AddAuxiliary();

  for (int i = 0; i < 3; ++i) {
    layers_.push_back(std::unique_ptr<MixerInput>(new MixerInput(sigmoid_,
        1.0e-4)));
    mixers_.push_back(std::vector<std::unique_ptr<Mixer>>());
  }

  unsigned long long input_size = GetNumModels();
  layers_[0]->SetNumModels(input_size);
  std::vector<std::vector<double>> model_params = {{0, 8, 0.005},
      {0, 8, 0.0005}, {1, 8, 0.005}, {1, 8, 0.0005}, {2, 4, 0.005},
      {3, 2, 0.002}};
  for (const auto& params : model_params) {
    const Context& context = manager_.AddContext(std::unique_ptr<Context>(
        new ContextHash(manager_.bit_context_, params[0], params[1])));
    const BitContext& bit_context = manager_.AddBitContext(std::unique_ptr
        <BitContext>(new BitContext(manager_.long_bit_context_,
        context.GetContext(), context.Size())));
    AddMixer(0, bit_context.GetContext(), params[2]);
  }

  model_params = {{0, 0.001}, {2, 0.002}, {3, 0.005}};
  for (const auto& params : model_params) {
    AddMixer(0, manager_.recent_bytes_[params[0]], params[1]);
  }
  AddMixer(0, manager_.zero_context_, 0.00005);
  AddMixer(0, manager_.line_break_, 0.0007);
  AddMixer(0, manager_.longest_match_, 0.0005);
  AddMixer(0, manager_.wrt_context_, 0.002);
  AddMixer(0, manager_.auxiliary_context_, 0.0005);

  std::vector<int> map(256, 0);
  for (int i = 0; i < 256; ++i) {
    map[i] = (i < 1) + (i < 32) + (i < 64) + (i < 128) + (i < 255) +
      (i < 142) + (i < 138) + (i < 140) + (i < 137) + (i < 97);
  }
  const Context& interval1 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 8)));
  AddMixer(0, interval1.GetContext(), 0.001);

  for (int i = 0; i < 256; ++i) {
    map[i] = (i < 41) + (i < 92) + (i < 124) + (i < 58) +
        (i < 11) + (i < 46) + (i < 36) + (i < 47) +
        (i < 64) + (i < 4) + (i < 61) + (i < 97) +
        (i < 125) + (i < 45) + (i < 48);
  }
  const Context& interval2 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 8)));
  AddMixer(0, interval2.GetContext(), 0.001);

  for (int i = 0; i < 256; ++i) map[i] = 0;
  for (int i = 'a'; i <= 'z'; ++i) map[i] = 1;
  for (int i = 'A'; i <= 'Z'; ++i) map[i] = 1;
  for (int i = '0'; i <= '9'; ++i) map[i] = 1;
  for (int i = 0x80; i < 256; ++i) map[i] = 1;
  const Context& interval3 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 7)));
  AddMixer(0, interval3.GetContext(), 0.001);
  const BitContext& bit_context5 = manager_.AddBitContext(std::unique_ptr
      <BitContext>(new BitContext(manager_.long_bit_context_,
      interval3.GetContext(), interval3.Size())));
  AddMixer(0, bit_context5.GetContext(), 0.005);

  for (int i = 0; i < 256; ++i) map[i] = 0;
  for (int i = 0x30; i < 0x60; ++i) map[i] = 1;
  for (int i = 0x60; i < 0xD0; ++i) map[i] = 2;
  for (int i = 0xD0; i < 256; ++i) map[i] = 3;
  const Context& interval4 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 10)));
  AddMixer(0, interval4.GetContext(), 0.001);
  const Context& interval5 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 15)));
  AddMixer(0, interval5.GetContext(), 0.001);
  const Context& interval8 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 7)));
  const BitContext& bit_context4 = manager_.AddBitContext(std::unique_ptr
      <BitContext>(new BitContext(manager_.long_bit_context_,
      interval8.GetContext(), interval8.Size())));
  AddMixer(0, bit_context4.GetContext(), 0.005);

  for (int i = 0; i < 256; ++i) map[i] = 0;
  for (int i = 0x20; i <= 0x7E; ++i) map[i] = 1;
  for (int i = 'a'; i <= 'z'; ++i) map[i] = 2;
  for (int i = '0'; i <= '9'; ++i) map[i] = 3;
  for (int i = 0x80; i <= 0xCF; ++i) map[i] = 4;
  for (int i = 0xD0; i <= 0xEF; ++i) map[i] = 5;
  for (int i = 0xF0; i <= 0xFF; ++i) map[i] = 6;
  map[' '] = 7;
  const Context& interval6 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 9)));
  AddMixer(0, interval6.GetContext(), 0.001);
  const Context& interval7 = manager_.AddContext(std::unique_ptr<Context>(
      new IntervalHash(manager_.bit_context_, map, 8, 7, 2)));
  AddMixer(0, interval7.GetContext(), 0.001);
  const Context& interval9 = manager_.AddContext(std::unique_ptr<Context>(
      new Interval(manager_.bit_context_, map, 7)));
  const BitContext& bit_context6 = manager_.AddBitContext(std::unique_ptr
      <BitContext>(new BitContext(manager_.long_bit_context_,
      interval9.GetContext(), interval9.Size())));
  AddMixer(0, bit_context6.GetContext(), 0.005);

  const BitContext& bit_context1 = manager_.AddBitContext(std::unique_ptr
      <BitContext>(new BitContext(manager_.long_bit_context_,
      manager_.recent_bytes_[1], 256)));
  AddMixer(0, bit_context1.GetContext(), 0.005);

  const Context& combined1 = manager_.AddContext(std::unique_ptr
      <Context>(new CombinedContext(manager_.recent_bytes_[1],
      manager_.recent_bytes_[0], 256, 256)));
  AddMixer(0, combined1.GetContext(), 0.005);

  const Context& combined2 = manager_.AddContext(std::unique_ptr
      <Context>(new CombinedContext(manager_.recent_bytes_[2],
      manager_.recent_bytes_[1], 256, 256)));
  AddMixer(0, combined2.GetContext(), 0.003);

  input_size = mixers_[0].size() + auxiliary_.size();
  layers_[1]->SetNumModels(input_size);

  AddMixer(1, manager_.zero_context_, 0.005);
  AddMixer(1, manager_.zero_context_, 0.0005);
  AddMixer(1, manager_.long_bit_context_, 0.005);
  AddMixer(1, manager_.long_bit_context_, 0.0005);
  AddMixer(1, manager_.long_bit_context_, 0.00001);
  AddMixer(1, manager_.recent_bytes_[0], 0.005);
  AddMixer(1, manager_.recent_bytes_[1], 0.005);
  AddMixer(1, manager_.recent_bytes_[2], 0.005);
  AddMixer(1, manager_.longest_match_, 0.0005);
  AddMixer(1, manager_.wrt_context_, 0.002);
  AddMixer(1, interval1.GetContext(), 0.001);
  AddMixer(1, interval2.GetContext(), 0.001);
  AddMixer(1, interval3.GetContext(), 0.001);
  AddMixer(1, interval4.GetContext(), 0.001);
  AddMixer(1, interval5.GetContext(), 0.001);
  AddMixer(1, interval6.GetContext(), 0.001);
  AddMixer(1, interval7.GetContext(), 0.001);
  AddMixer(1, bit_context4.GetContext(), 0.001);
  AddMixer(1, bit_context5.GetContext(), 0.001);
  AddMixer(1, bit_context6.GetContext(), 0.001);

  input_size = mixers_[0].size() + mixers_[1].size() + auxiliary_.size();
  layers_[2]->SetNumModels(input_size);
  AddMixer(2, manager_.zero_context_, 0.0003);
}

float Predictor::Predict() {
  unsigned int input_index = 0;
  for (unsigned int i = 0; i < models_.size(); ++i) {
    const std::valarray<float>& outputs = models_[i]->Predict();
    for (unsigned int j = 0; j < outputs.size(); ++j) {
      layers_[0]->SetInput(input_index, outputs[j]);
      ++input_index;
    }
  }

  for (unsigned int i = 0; i < byte_models_.size(); ++i) {
    const std::valarray<float>& outputs = byte_models_[i]->Predict();
    for (unsigned int j = 0; j < outputs.size(); ++j) {
      layers_[0]->SetInput(input_index, outputs[j]);
      ++input_index;
    }
  }
  float byte_mixer_override = -1;
  for (unsigned int i = 0; i < byte_mixers_.size(); ++i) {
    const std::valarray<float>& outputs = byte_mixers_[i]->Predict();
    for (unsigned int j = 0; j < outputs.size(); ++j) {
      float p = outputs[j];
      if (p == 0 || p == 1) byte_mixer_override = p;
      layers_[0]->SetInput(input_index, p);
      ++input_index;
    }
  }
  float auxiliary_average = 0;
  for (unsigned int i = 0; i < auxiliary_.size(); ++i) {
    auxiliary_average += Sigmoid::Logistic(layers_[0]->Inputs()[auxiliary_[i]]);
  }
  auxiliary_average /= auxiliary_.size();
  manager_.auxiliary_context_ = auxiliary_average * 15;

  for (unsigned int i = 0; i < mixers_[0].size(); ++i) {
    float p = mixers_[0][i]->Mix();
    layers_[0]->SetExtraInput(p);
    layers_[1]->SetStretchedInput(i, p);
    layers_[2]->SetStretchedInput(i, p);
  }
  layers_[0]->ClearExtraInputs();
  for (unsigned int i = 0; i < auxiliary_.size(); ++i) {
    float p = layers_[0]->Inputs()[auxiliary_[i]];
    layers_[1]->SetStretchedInput(mixers_[0].size() + i, p);
    layers_[2]->SetStretchedInput(mixers_[0].size() + mixers_[1].size() + i, p);
  }
  for (unsigned int i = 0; i < mixers_[1].size(); ++i) {
    float p = mixers_[1][i]->Mix();
    layers_[1]->SetExtraInput(p);
    layers_[2]->SetStretchedInput(mixers_[0].size() + i, p);
  }
  layers_[1]->ClearExtraInputs();
  float p = Sigmoid::Logistic(mixers_[2][0]->Mix());
  p = sse_.Predict(p);
  if (byte_mixer_override >= 0) {
    return byte_mixer_override;
  }
  return p;
}

void Predictor::Perceive(int bit) {
  for (const auto& model : models_) {
    model->Perceive(bit);
  }
  for (const auto& model : byte_models_) {
    model->Perceive(bit);
  }
  for (const auto& byte_mixer : byte_mixers_) {
    byte_mixer->Perceive(bit);
  }
  for (unsigned int i = 0; i < mixers_.size(); ++i) {
    for (const auto& mixer : mixers_[i]) {
      mixer->Perceive(bit);
    }
  }
  sse_.Perceive(bit);

  bool byte_update = false;
  if (manager_.bit_context_ >= 128) byte_update = true;

  manager_.UpdateContexts(bit);
  if (byte_update) {
    for (const auto& model : models_) {
      model->ByteUpdate();
    }
    for (const auto& model : byte_models_) {
      model->ByteUpdate();
    }
    for (unsigned int i = 0; i < byte_models_.size(); ++i) {
      const std::valarray<float>& p = byte_models_[i]->BytePredict();
      for (const auto& byte_mixer : byte_mixers_) {
        for (unsigned int j = 0; j < 256; ++j) {
          byte_mixer->SetInput(j, p[j]);
        }
      }
    }
    for (const auto& byte_mixer : byte_mixers_) {
      byte_mixer->ByteUpdate();
    }
    manager_.bit_context_ = 1;
  }
}

void Predictor::Pretrain(int bit) {
  for (const auto& model : models_) {
    model->Predict();
  }
  for (const auto& model : models_) {
    model->Perceive(bit);
  }
  bool byte_update = false;
  if (manager_.bit_context_ >= 128) byte_update = true;
  manager_.UpdateContexts(bit);
  if (byte_update) {
    for (const auto& model : models_) {
      model->ByteUpdate();
    }
    manager_.bit_context_ = 1;
  }
}
