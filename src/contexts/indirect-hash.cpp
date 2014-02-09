#include "indirect-hash.h"

IndirectHash::IndirectHash(const unsigned int& bit_context, unsigned int order1,
    unsigned int hash_size1, unsigned int order2, unsigned int hash_size2) :
    byte_(bit_context), context1_(0), hash_size1_(hash_size1),
    hash_size2_(hash_size2), size1_(0), hashes_(0, 0) {
  context_ = 0;
  size1_ = (unsigned long long)1 << (hash_size1 * order1);
  size_ = (unsigned long long)1 << (hash_size2 * order2);
  hashes_.resize(size1_, 0);
}

void IndirectHash::Update() {
  hashes_[context1_] = (context_ * (1 << hash_size2_) + byte_) % size_;
  context1_ = (context1_ * (1 << hash_size1_) + byte_) % size1_;
  context_ = hashes_[context1_];
}

bool IndirectHash::IsEqual(Context* c) {
  IndirectHash* p = dynamic_cast<IndirectHash*>(c);
  if (!p) return false;
  if (size_ == p->size_ && size1_ == p->size1_ &&
      hash_size1_ == p->hash_size1_ && hash_size2_ == p->hash_size2_) {
    return true;
  }
  return false;
}
