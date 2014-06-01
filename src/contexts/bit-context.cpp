#include "bit-context.h"

BitContext::BitContext(const unsigned long long& bit_context,
    const unsigned long long& byte_context,
    unsigned long long byte_context_size) : context_(0),
    size_(256 * byte_context_size), bit_context_(bit_context),
    byte_context_(byte_context) {}

void BitContext::Update() {
  context_ = (byte_context_<<8) + bit_context_;
}
