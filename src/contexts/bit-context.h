#ifndef BIT_CONTEXT_H
#define BIT_CONTEXT_H

class BitContext {
 public:
  BitContext(const unsigned long long& bit_context, const unsigned long long&
      byte_context, unsigned long long byte_context_size);
  void Update();

  unsigned long long context_;
  unsigned long long size_;

 private:
  const unsigned long long& bit_context_;
  const unsigned long long& byte_context_;
};

#endif
