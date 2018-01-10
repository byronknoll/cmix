#include "interval.h"

Interval::Interval(const unsigned int& bit_context,
    const std::vector<int>& map, unsigned int num_bits) : byte_(bit_context),
    map_(map) {
  context_ = 0;
  int max_value = 0;
  for (unsigned int i = 0; i < map.size(); ++i) {
    if (map[i] > max_value) max_value = map[i];
  }
  shift_ = 1;
  while ((1 << shift_) <= max_value) ++shift_;
  size_ = 1 << num_bits;
  mask_ = size_ - 1;
}

void Interval::Update() {
  context_ = mask_ & ((context_ << shift_) + map_[byte_]);
}

bool Interval::IsEqual(Context* c) {
  Interval* p = dynamic_cast<Interval*>(c);
  if (!p) return false;
  if (p->size_ != size_) return false;
  for (int i = 0; i < 256; ++i) {
    if (p->map_[i] != map_[i]) return false;
  }
  return true;
}
