#include "bit-context.h"

BitContext::BitContext(const unsigned long long& bit_context,
    const unsigned long long& byte_context,
    unsigned long long byte_context_size) : bit_context_(bit_context),
    byte_context_(byte_context) {
  context_ = 0;
  size_ = 256 * byte_context_size;
}

void BitContext::Update() {
  context_ = (byte_context_<<8) + bit_context_;
}

bool BitContext::IsEqual(Context* c) {
  BitContext* p = dynamic_cast<BitContext*>(c);
  if (!p) return false;
  if (&(p->byte_context_) != &byte_context_) return false;
  if (&(p->bit_context_) != &bit_context_) return false;
  return true;
}
