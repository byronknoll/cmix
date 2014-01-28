#include "context-hash.h"

ContextHash::ContextHash(const unsigned int& bit_context, unsigned int order,
    unsigned int hash_size) : byte_(bit_context), hash_size_(hash_size) {
  context_ = 0;
  size_ = (unsigned long long)1 << (hash_size * order);
}

void ContextHash::Update() {
  context_ = (context_ * (1 << hash_size_) + byte_) % size_;
}

bool ContextHash::IsEqual(Context* c) {
  ContextHash* p = dynamic_cast<ContextHash*>(c);
  if (!p) return false;
  if (size_ == p->size_ && hash_size_ == p->hash_size_) return true;
  return false;
}
