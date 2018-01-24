#include "interval-hash.h"

IntervalHash::IntervalHash(const unsigned int& bit_context,
    const std::vector<int>& map, unsigned int num_bits, unsigned int order,
    unsigned int hash_size) : byte_(bit_context), map_(map),
    hash_size_(hash_size), interval_(0) {
  context_ = 0;
  int max_value = 0;
  for (unsigned int i = 0; i < map.size(); ++i) {
    if (map[i] > max_value) max_value = map[i];
  }
  shift_ = 1;
  while ((1 << shift_) <= max_value) ++shift_;
  mask_ = (1 << num_bits) - 1;
  size_ = (unsigned long long)1 << (hash_size * order);
}

void IntervalHash::Update() {
  interval_ = mask_ & ((interval_ << shift_) + map_[byte_]);
  context_ = (context_ * (1 << hash_size_) + interval_) % size_;
}

bool IntervalHash::IsEqual(Context* c) {
  IntervalHash* p = dynamic_cast<IntervalHash*>(c);
  if (!p) return false;
  for (int i = 0; i < 256; ++i) {
    if (p->map_[i] != map_[i]) return false;
  }
  if (size_ != p->size_ || hash_size_ != p->hash_size_ ||
    mask_ != p->mask_) return false;
  return true;
}
