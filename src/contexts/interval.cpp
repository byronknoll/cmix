#include "interval.h"

Interval::Interval(const unsigned int& bit_context,
    const std::vector<int>& map) : byte_(bit_context), map_(map) {
  context_ = 0;
  int max_value = 0;
  for (unsigned int i = 0; i < map.size(); ++i) {
    if (map[i] > max_value) max_value = map[i];
  }
  size_ = (max_value << 4) + max_value + 1;
}

void Interval::Update() {
  context_ = 0xFF & ((context_ << 4) + map_[byte_]);
}

bool Interval::IsEqual(Context* c) {
  Interval* p = dynamic_cast<Interval*>(c);
  if (!p) return false;
  for (int i = 0; i < 256; ++i) {
    if (p->map_[i] != map_[i]) return false;
  }
  return true;
}
