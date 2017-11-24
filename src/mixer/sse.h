#ifndef SSE_H
#define SSE_H

namespace SSE_sh {
  struct M_T1;
}

class SSE {
 public:
  SSE();
  ~SSE();
  float Predict(float input);
  void Perceive(int bit);

 private:
  SSE_sh::M_T1* sse_;
};

#endif
