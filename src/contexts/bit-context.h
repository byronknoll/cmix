#ifndef BIT_CONTEXT_H
#define BIT_CONTEXT_H

class BitContext {
 public:
  BitContext(const unsigned long long& bit_context, const unsigned long long&
      byte_context, unsigned long long byte_context_size);
  void Update();
  const unsigned long long& GetContext() const { return context_; }
  unsigned long long Size() const { return size_; }

 private:
  unsigned long long context_, size_;
  const unsigned long long& bit_context_;
  const unsigned long long& byte_context_;
};

#endif
