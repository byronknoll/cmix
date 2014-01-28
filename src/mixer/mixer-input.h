#ifndef MIXER_INPUT_H
#define MIXER_INPUT_H

#include <vector>
#include <array>

class MixerInput {
 public:
  MixerInput(float eps, int stretch_size, int squash_size);
  void SetNumModels(int num_models);
  void SetInput(int index, float p);
  float Stretch(float p) const;
  float Squash(float p) const;

  std::vector<float> inputs_;

 private:
  float SlowStretch(float p);
  float SlowSquash(float p);

  float min_;
  float max_;
  int stretch_size_;
  int squash_size_;
  std::vector<float> stretch_table_;
  std::vector<float> squash_table_;
};

#endif
