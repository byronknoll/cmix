#ifndef SSE2_H
#define SSE2_H

namespace SSE_sh {
  class M_T1;
}

class SSE2 {
 public:
  SSE2();
  ~SSE2();
  float Process(float input);
  void Perceive(int bit);

 private:
  SSE_sh::M_T1* sse_;
};

#endif
