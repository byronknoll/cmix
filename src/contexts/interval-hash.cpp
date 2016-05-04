#include "interval-hash.h"

IntervalHash::IntervalHash(const unsigned int& bit_context,
    const std::vector<int>& map, unsigned int order, unsigned int hash_size) :
    byte_(bit_context), map_(map), hash_size_(hash_size), interval_(0) {
  context_ = 0;
  size_ = (unsigned long long)1 << (hash_size * order);
}

void IntervalHash::Update() {
  interval_ = 0xFF & ((interval_ << 4) + map_[byte_]);
  context_ = (context_ * (1 << hash_size_) + interval_) % size_;
}

bool IntervalHash::IsEqual(Context* c) {
  IntervalHash* p = dynamic_cast<IntervalHash*>(c);
  if (!p) return false;
  for (int i = 0; i < 256; ++i) {
    if (p->map_[i] != map_[i]) return false;
  }
  if (size_ != p->size_ || hash_size_ != p->hash_size_) return false;
  return true;
}
