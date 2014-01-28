#ifndef CONTEXT_H
#define CONTEXT_H

class Context {
 public:
  virtual ~Context() {}
  virtual void Update() {}
  virtual bool IsEqual(Context* c) {return false;}

  unsigned long long context_;
  unsigned long long size_;
};

#endif
