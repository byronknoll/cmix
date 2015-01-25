// This code is a hybrid of paq8l and paq8pxd_11 (released by Kaido Orav).

/*
    Copyright (C) 2006 Matt Mahoney, Serge Osnach, Alexander Ratushnyak,
    Bill Pettis, Przemyslaw Skibinski, Matthew Fite, wowtiger, Andrew Paterson,

    LICENSE

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details at
    Visit <http://www.gnu.org/copyleft/gpl.html>.
*/

#include "paq8l.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#define NDEBUG  // remove for debugging (turns on Array bound checks)

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#endif

#ifdef WINDOWS
#include <windows.h>
#endif

#ifndef DEFAULT_OPTION
#define DEFAULT_OPTION 9
#endif

#ifndef NOASM
#define NOASM
#endif

namespace {
typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;

#ifndef WINDOWS
inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}
#endif

void quit(const char* message=0) {
  throw message;
}

//////////////////////////// Array ////////////////////////////

template <class T, int ALIGN=0> class Array {
private:
  U32 n;
  U32 reserved;
  char *ptr;
  T* data;
  void create(U32 i);
public:
  explicit Array(U32 i=0) {create(i);}
  ~Array();
  T& operator[](U32 i) {
    return data[i];
  }
  const T& operator[](U32 i) const {
    return data[i];
  }
  U32 size() const {return n;}
  void resize(U32 i);
  void pop_back() {if (n>0) --n;}
  void push_back(const T& x);
private:
  Array(const Array&);
  Array& operator=(const Array&);
};

template<class T, int ALIGN> void Array<T, ALIGN>::resize(U32 i) {
  if (i<=reserved) {
    n=i;
    return;
  }
  char *saveptr=ptr;
  T *savedata=data;
  int saven=n;
  create(i);
  if (saveptr) {
    if (savedata) {
      memcpy(data, savedata, sizeof(T)*min(i, saven));
    }
    free(saveptr);
  }
}

template<class T, int ALIGN> void Array<T, ALIGN>::create(U32 i) {
  n=reserved=i;
  if (i<=0) {
    data=0;
    ptr=0;
    return;
  }
  const U32 sz=ALIGN+n*sizeof(T);
  ptr = (char*)calloc(sz, 1);
  if (!ptr) quit("Out of memory");
  data = (ALIGN ? (T*)(ptr+ALIGN-(((long long)ptr)&(ALIGN-1))) : (T*)ptr);
}

template<class T, int ALIGN> Array<T, ALIGN>::~Array() {
  free(ptr);
}

template<class T, int ALIGN> void Array<T, ALIGN>::push_back(const T& x) {
  if (n==reserved) {
    int saven=n;
    resize(max(1, n*2));
    n=saven;
  }
  data[n++]=x;
}

/////////////////////////// String /////////////////////////////

class String: public Array<char> {
public:
  const char* c_str() const {return &(*this)[0];}
  void operator=(const char* s) {
    resize(strlen(s)+1);
    strcpy(&(*this)[0], s);
  }
  void operator+=(const char* s) {
    pop_back();
    while (*s) push_back(*s++);
    push_back(0);
  }
  String(const char* s=""): Array<char>(1) {
    (*this)+=s;
  }
};


//////////////////////////// rnd ///////////////////////////////

class Random{
  Array<U32> table;
  int i;
public:
  Random(): table(64) {
    table[0]=123456789;
    table[1]=987654321;
    for(int j=0; j<62; j++) table[j+2]=table[j+1]*11+table[j]*23/16;
    i=0;
  }
  U32 operator()() {
    return ++i, table[i&63]=table[(i-24)&63]^table[(i-55)&63];
  }
} rnd;

////////////////////////////// Buf /////////////////////////////

int pos;

class Buf {
  Array<U8> b;
public:
  Buf(U32 i=0): b(i) {}
  void setsize(U32 i) {
    if (!i) return;
    b.resize(i);
  }
  U8& operator[](U32 i) {
    return b[i&(b.size()-1)];
  }
  int operator()(U32 i) const {
    return b[(pos-i)&(b.size()-1)];
  }
  U32 size() const {
    return b.size();
  }
};

/////////////////////// Global context /////////////////////////

int level=DEFAULT_OPTION;  // Compression level 0 to 9
#define MEM (0x10000<<level)
int y=0;  // Last bit, 0 or 1, set by encoder

int c0=1;
U32 c4=0;
int bpos=0;
int blpos=0;
Buf buf;

///////////////////////////// ilog //////////////////////////////

class Ilog {
  Array<U8> t;
public:
  int operator()(U16 x) const {return t[x];}
  Ilog();
} ilog;

Ilog::Ilog(): t(65536) {
  U32 x=14155776;
  for (int i=2; i<65536; ++i) {
    x+=774541002/(i*2-1);
    t[i]=x>>24;
  }
}

inline int llog(U32 x) {
  if (x>=0x1000000)
    return 256+ilog(x>>16);
  else if (x>=0x10000)
    return 128+ilog(x>>8);
  else
    return ilog(x);
}

///////////////////////// state table ////////////////////////

static const U8 State_table[256][4]={
  {  1,  2, 0, 0},{  3,  5, 1, 0},{  4,  6, 0, 1},{  7, 10, 2, 0},
  {  8, 12, 1, 1},{  9, 13, 1, 1},{ 11, 14, 0, 2},{ 15, 19, 3, 0},
  { 16, 23, 2, 1},{ 17, 24, 2, 1},{ 18, 25, 2, 1},{ 20, 27, 1, 2},
  { 21, 28, 1, 2},{ 22, 29, 1, 2},{ 26, 30, 0, 3},{ 31, 33, 4, 0},
  { 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1},
  { 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2},
  { 34, 37, 2, 2},{ 34, 37, 2, 2},{ 36, 39, 1, 3},{ 36, 39, 1, 3},
  { 36, 39, 1, 3},{ 36, 39, 1, 3},{ 38, 40, 0, 4},{ 41, 43, 5, 0},
  { 42, 45, 4, 1},{ 42, 45, 4, 1},{ 44, 47, 3, 2},{ 44, 47, 3, 2},
  { 46, 49, 2, 3},{ 46, 49, 2, 3},{ 48, 51, 1, 4},{ 48, 51, 1, 4},
  { 50, 52, 0, 5},{ 53, 43, 6, 0},{ 54, 57, 5, 1},{ 54, 57, 5, 1},
  { 56, 59, 4, 2},{ 56, 59, 4, 2},{ 58, 61, 3, 3},{ 58, 61, 3, 3},
  { 60, 63, 2, 4},{ 60, 63, 2, 4},{ 62, 65, 1, 5},{ 62, 65, 1, 5},
  { 50, 66, 0, 6},{ 67, 55, 7, 0},{ 68, 57, 6, 1},{ 68, 57, 6, 1},
  { 70, 73, 5, 2},{ 70, 73, 5, 2},{ 72, 75, 4, 3},{ 72, 75, 4, 3},
  { 74, 77, 3, 4},{ 74, 77, 3, 4},{ 76, 79, 2, 5},{ 76, 79, 2, 5},
  { 62, 81, 1, 6},{ 62, 81, 1, 6},{ 64, 82, 0, 7},{ 83, 69, 8, 0},
  { 84, 71, 7, 1},{ 84, 71, 7, 1},{ 86, 73, 6, 2},{ 86, 73, 6, 2},
  { 44, 59, 5, 3},{ 44, 59, 5, 3},{ 58, 61, 4, 4},{ 58, 61, 4, 4},
  { 60, 49, 3, 5},{ 60, 49, 3, 5},{ 76, 89, 2, 6},{ 76, 89, 2, 6},
  { 78, 91, 1, 7},{ 78, 91, 1, 7},{ 80, 92, 0, 8},{ 93, 69, 9, 0},
  { 94, 87, 8, 1},{ 94, 87, 8, 1},{ 96, 45, 7, 2},{ 96, 45, 7, 2},
  { 48, 99, 2, 7},{ 48, 99, 2, 7},{ 88,101, 1, 8},{ 88,101, 1, 8},
  { 80,102, 0, 9},{103, 69,10, 0},{104, 87, 9, 1},{104, 87, 9, 1},
  {106, 57, 8, 2},{106, 57, 8, 2},{ 62,109, 2, 8},{ 62,109, 2, 8},
  { 88,111, 1, 9},{ 88,111, 1, 9},{ 80,112, 0,10},{113, 85,11, 0},
  {114, 87,10, 1},{114, 87,10, 1},{116, 57, 9, 2},{116, 57, 9, 2},
  { 62,119, 2, 9},{ 62,119, 2, 9},{ 88,121, 1,10},{ 88,121, 1,10},
  { 90,122, 0,11},{123, 85,12, 0},{124, 97,11, 1},{124, 97,11, 1},
  {126, 57,10, 2},{126, 57,10, 2},{ 62,129, 2,10},{ 62,129, 2,10},
  { 98,131, 1,11},{ 98,131, 1,11},{ 90,132, 0,12},{133, 85,13, 0},
  {134, 97,12, 1},{134, 97,12, 1},{136, 57,11, 2},{136, 57,11, 2},
  { 62,139, 2,11},{ 62,139, 2,11},{ 98,141, 1,12},{ 98,141, 1,12},
  { 90,142, 0,13},{143, 95,14, 0},{144, 97,13, 1},{144, 97,13, 1},
  { 68, 57,12, 2},{ 68, 57,12, 2},{ 62, 81, 2,12},{ 62, 81, 2,12},
  { 98,147, 1,13},{ 98,147, 1,13},{100,148, 0,14},{149, 95,15, 0},
  {150,107,14, 1},{150,107,14, 1},{108,151, 1,14},{108,151, 1,14},
  {100,152, 0,15},{153, 95,16, 0},{154,107,15, 1},{108,155, 1,15},
  {100,156, 0,16},{157, 95,17, 0},{158,107,16, 1},{108,159, 1,16},
  {100,160, 0,17},{161,105,18, 0},{162,107,17, 1},{108,163, 1,17},
  {110,164, 0,18},{165,105,19, 0},{166,117,18, 1},{118,167, 1,18},
  {110,168, 0,19},{169,105,20, 0},{170,117,19, 1},{118,171, 1,19},
  {110,172, 0,20},{173,105,21, 0},{174,117,20, 1},{118,175, 1,20},
  {110,176, 0,21},{177,105,22, 0},{178,117,21, 1},{118,179, 1,21},
  {110,180, 0,22},{181,115,23, 0},{182,117,22, 1},{118,183, 1,22},
  {120,184, 0,23},{185,115,24, 0},{186,127,23, 1},{128,187, 1,23},
  {120,188, 0,24},{189,115,25, 0},{190,127,24, 1},{128,191, 1,24},
  {120,192, 0,25},{193,115,26, 0},{194,127,25, 1},{128,195, 1,25},
  {120,196, 0,26},{197,115,27, 0},{198,127,26, 1},{128,199, 1,26},
  {120,200, 0,27},{201,115,28, 0},{202,127,27, 1},{128,203, 1,27},
  {120,204, 0,28},{205,115,29, 0},{206,127,28, 1},{128,207, 1,28},
  {120,208, 0,29},{209,125,30, 0},{210,127,29, 1},{128,211, 1,29},
  {130,212, 0,30},{213,125,31, 0},{214,137,30, 1},{138,215, 1,30},
  {130,216, 0,31},{217,125,32, 0},{218,137,31, 1},{138,219, 1,31},
  {130,220, 0,32},{221,125,33, 0},{222,137,32, 1},{138,223, 1,32},
  {130,224, 0,33},{225,125,34, 0},{226,137,33, 1},{138,227, 1,33},
  {130,228, 0,34},{229,125,35, 0},{230,137,34, 1},{138,231, 1,34},
  {130,232, 0,35},{233,125,36, 0},{234,137,35, 1},{138,235, 1,35},
  {130,236, 0,36},{237,125,37, 0},{238,137,36, 1},{138,239, 1,36},
  {130,240, 0,37},{241,125,38, 0},{242,137,37, 1},{138,243, 1,37},
  {130,244, 0,38},{245,135,39, 0},{246,137,38, 1},{138,247, 1,38},
  {140,248, 0,39},{249,135,40, 0},{250, 69,39, 1},{ 80,251, 1,39},
  {140,252, 0,40},{249,135,41, 0},{250, 69,40, 1},{ 80,251, 1,40},
  {140,252, 0,41}};  // 252, 253-255 are reserved

#define nex(state,sel) State_table[state][sel]

///////////////////////////// Squash //////////////////////////////

class Squash {
  Array<U16> t;
public:
  Squash();
  int operator()(int p) const {
    if (p>2047) return 4095;
  if (p<-2047) return 0;
    return t[p+2048];
  }
} squash;

Squash::Squash(): t(4096) {
                  int ts[33]={
    1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
    1546,2047,2549,2994,3348,3607,3785,3901,3975,4022,
    4050,4068,4079,4085,4089,4092,4093,4094};
    int w,d;
  for (int i=-2047; i<=2047; ++i){
    w=i&127;
  d=(i>>7)+16;
  t[i+2048]=(ts[d]*(128-w)+ts[(d+1)]*w+64) >> 7;
    }
}

//////////////////////////// Stretch ///////////////////////////////

class Stretch {
  Array<short> t;
public:
  Stretch();
  int operator()(int p) const {
    return t[p];
  }
} stretch;

Stretch::Stretch(): t(4096) {
  int pi=0;
  for (int x=-2047; x<=2047; ++x) {  // invert squash()
    int i=squash(x);
    for (int j=pi; j<=i; ++j)
      t[j]=x;
    pi=i+1;
  }
  t[4095]=2047;
}

//////////////////////////// Mixer /////////////////////////////

#if !defined(__GNUC__)

#if (2 == _M_IX86_FP) // 2 if /arch:SSE2 was used.
# define __SSE2__
#elif (1 == _M_IX86_FP) // 1 if /arch:SSE was used.
# define __SSE__
#endif

#endif /* __GNUC__ */

static int dot_product (const short* const t, const short* const w, int n);

static void train (const short* const t, short* const w, int n, const int e);
#if defined(__AVX2__) // fast
#include<immintrin.h>
#define OPTIMIZE "AVX2-"
static int dot_product (const short* const t, const short* const w, int n) {
  __m256i sum = _mm256_setzero_si256 ();
  while ((n -= 16) >= 0) {
    __m256i tmp = _mm256_madd_epi16 (*(__m256i *) &t[n], *(__m256i *) &w[n]);
    tmp = _mm256_srai_epi32 (tmp, 8);
    sum = _mm256_add_epi32 (sum, tmp);
  } 
  __m128i low = _mm_add_epi32 (_mm256_extracti128_si256(sum,0),_mm256_extracti128_si256(sum,1)); 
  low = _mm_add_epi32 (low, _mm_srli_si128 (low, 8));
  low = _mm_add_epi32 (low, _mm_srli_si128 (low, 4));
  return _mm_cvtsi128_si32 (low);
}

static void train (const short* const t, short* const w, int n, const int e) {
  if (e) {
    const __m256i one = _mm256_set1_epi16 (1);
    const __m256i err = _mm256_set1_epi16 (short(e));
    while ((n -= 16) >= 0) {
      __m256i tmp = _mm256_adds_epi16 (*(__m256i *) &t[n], *(__m256i *) &t[n]);
      tmp = _mm256_mulhi_epi16 (tmp, err);
      tmp = _mm256_adds_epi16 (tmp, one);
      tmp = _mm256_srai_epi16 (tmp, 1);
      tmp = _mm256_adds_epi16 (tmp, *(__m256i *) &w[n]);
      *(__m256i *) &w[n] = tmp;
    }
  }
}

#elif defined(__SSE2__) // faster
#include <emmintrin.h>
#define OPTIMIZE "SSE2-"

static int dot_product (const short* const t, const short* const w, int n) {
  __m128i sum = _mm_setzero_si128 ();
  while ((n -= 8) >= 0) {
    __m128i tmp = _mm_madd_epi16 (*(__m128i *) &t[n], *(__m128i *) &w[n]);
    tmp = _mm_srai_epi32 (tmp, 8);
    sum = _mm_add_epi32 (sum, tmp);
  }
  sum = _mm_add_epi32 (sum, _mm_srli_si128 (sum, 8));
  sum = _mm_add_epi32 (sum, _mm_srli_si128 (sum, 4));
  return _mm_cvtsi128_si32 (sum);
}

static void train (const short* const t, short* const w, int n, const int e) {
  if (e) {
    const __m128i one = _mm_set1_epi16 (1);
    const __m128i err = _mm_set1_epi16 (short(e));
    while ((n -= 8) >= 0) {
      __m128i tmp = _mm_adds_epi16 (*(__m128i *) &t[n], *(__m128i *) &t[n]);
      tmp = _mm_mulhi_epi16 (tmp, err);
      tmp = _mm_adds_epi16 (tmp, one);
      tmp = _mm_srai_epi16 (tmp, 1);
      tmp = _mm_adds_epi16 (tmp, *(__m128i *) &w[n]);
      *(__m128i *) &w[n] = tmp;
    }
  }
}

#elif defined(__SSE__) // fast
#include <xmmintrin.h>
#define OPTIMIZE "SSE-"

static int dot_product (const short* const t, const short* const w, int n) {
  __m64 sum = _mm_setzero_si64 ();
  while ((n -= 8) >= 0) {
    __m64 tmp = _mm_madd_pi16 (*(__m64 *) &t[n], *(__m64 *) &w[n]);
    tmp = _mm_srai_pi32 (tmp, 8);
    sum = _mm_add_pi32 (sum, tmp);

    tmp = _mm_madd_pi16 (*(__m64 *) &t[n + 4], *(__m64 *) &w[n + 4]);
    tmp = _mm_srai_pi32 (tmp, 8);
    sum = _mm_add_pi32 (sum, tmp);
  }
  sum = _mm_add_pi32 (sum, _mm_srli_si64 (sum, 32));
  const int retval = _mm_cvtsi64_si32 (sum);
  _mm_empty();
  return retval;
}

static void train (const short* const t, short* const w, int n, const int e) {
  if (e) {
    const __m64 one = _mm_set1_pi16 (1);
    const __m64 err = _mm_set1_pi16 (short(e));
    while ((n -= 8) >= 0) {
      __m64 tmp = _mm_adds_pi16 (*(__m64 *) &t[n], *(__m64 *) &t[n]);
      tmp = _mm_mulhi_pi16 (tmp, err);
      tmp = _mm_adds_pi16 (tmp, one);
      tmp = _mm_srai_pi16 (tmp, 1);
      tmp = _mm_adds_pi16 (tmp, *(__m64 *) &w[n]);
      *(__m64 *) &w[n] = tmp;

      tmp = _mm_adds_pi16 (*(__m64 *) &t[n + 4], *(__m64 *) &t[n + 4]);
      tmp = _mm_mulhi_pi16 (tmp, err);
      tmp = _mm_adds_pi16 (tmp, one);
      tmp = _mm_srai_pi16 (tmp, 1);
      tmp = _mm_adds_pi16 (tmp, *(__m64 *) &w[n + 4]);
      *(__m64 *) &w[n + 4] = tmp;
    }
    _mm_empty();
  }
}
#else
int dot_product(short *t, short *w, int n) {
  int sum=0;
  n=(n+15)&-16;
  for (int i=0; i<n; i+=2)
    sum+=(t[i]*w[i]+t[i+1]*w[i+1]) >> 8;
  return sum;
}

void train(short *t, short *w, int n, int err) {
  n=(n+15)&-16;
  for (int i=0; i<n; ++i) {
    int wt=w[i]+(((t[i]*err*2>>16)+1)>>1);
    if (wt<-32768) wt=-32768;
    if (wt>32767) wt=32767;
    w[i]=wt;
  }
}
#endif // slow!

std::vector<float> model_predictions(1167, 0.5);
unsigned int prediction_index = 0;
float conversion_factor = 1.0 / 4095;

class Mixer {
  const int N, M, S;
  Array<short, 16> tx;
  Array<short, 16> wx;
  Array<int> cxt;
  int ncxt;
  int base;
  int nx;
  Array<int> pr;
  Mixer* mp;
public:
  Mixer(int n, int m, int s=1, int w=0);

  // Adjust weights to minimize coding cost of last prediction
  void update() {
    for (int i=0; i<ncxt; ++i) {
      int err=((y<<12)-pr[i])*7;
      train(&tx[0], &wx[cxt[i]*N], nx, err);
    }
    nx=base=ncxt=0;
  }

  // Input x (call up to N times)
  void add(int x) {
    model_predictions[prediction_index] = squash(x) * conversion_factor;
    ++prediction_index;
    tx[nx++]=x;
  }

  // Set a context (call S times, sum of ranges <= M)
  void set(int cx, int range) {
    cxt[ncxt++]=base+cx;
    base+=range;
  }

  // predict next bit
  int p() {
    while (nx&7) tx[nx++]=0;  // pad
    if (mp) {  // combine outputs
      prediction_index = 0;
      mp->update();
      for (int i=0; i<ncxt; ++i) {
        pr[i]=squash((dot_product(&tx[0], &wx[cxt[i]*N], nx) * 9)>>9);
        mp->add(stretch(pr[i]));
      }
      mp->set(0, 1);
      return mp->p();
    }
    else {  // S=1 context
      int z = dot_product(&tx[0], &wx[0], nx);
      base = squash((z*16)>>13);
      return pr[0]=squash(z>>9);
    }
  }
  ~Mixer();
};

Mixer::~Mixer() {
  delete mp;
}


Mixer::Mixer(int n, int m, int s, int w):
    N((n+7)&-8), M(m), S(s), tx(N), wx(N*M),
    cxt(S), ncxt(0), base(0), nx(0), pr(S), mp(0) {
  for (int i=0; i<S; ++i)
    pr[i]=2048;
  for (int i=0; i<N*M; ++i)
    wx[i]=w;
  if (S>1) mp=new Mixer(S, 1, 1, 0x7fff);
}

//////////////////////////// APM //////////////////////////////

class APM {
  int index;     // last p, context
  const int N;   // number of contexts
  Array<U16> t;  // [N][33]:  p, context -> p
public:
  APM(int n);
  int p(int pr=2048, int cxt=0, int rate=7) {
    pr=stretch(pr);
    int g=(y<<16)+(y<<rate)-y-y;
    t[index] += (g-t[index]) >> rate;
    t[index+1] += (g-t[index+1]) >> rate;
    const int w=pr&127;  // interpolation weight (33 points)
    index=((pr+2048)>>7)+cxt*33;
    return (t[index]*(128-w)+t[index+1]*w) >> 11;
  }
};

// maps p, cxt -> p initially
APM::APM(int n): index(0), N(n), t(n*33) {
  for (int i=0; i<N; ++i)
    for (int j=0; j<33; ++j)
      t[i*33+j] = i==0 ? squash((j-16)*128)*16 : t[j];
}

//////////////////////////// StateMap //////////////////////////

class StateMap {
protected:
  int cxt;  // context
  Array<U16> t;
public:
  StateMap();
  int p(int cx) {
    t[cxt]+=((y<<16)-t[cxt]+128) >> 8;
    return t[cxt=cx] >> 4;
  }
};

StateMap::StateMap(): cxt(0), t(256) {
  for (int i=0; i<256; ++i) {
    int n0=nex(i,2);
    int n1=nex(i,3);
    if (n0==0) n1*=64;
    if (n1==0) n0*=64;
    t[i] = 65536*(n1+1)/(n0+n1+2);
  }
}

//////////////////////////// hash //////////////////////////////

inline U32 hash(U32 a, U32 b, U32 c=0xffffffff, U32 d=0xffffffff,
    U32 e=0xffffffff) {
  U32 h=a*200002979u+b*30005491u+c*50004239u+d*70004807u+e*110002499u;
  return h^h>>9^a>>2^b>>3^c>>4^d>>5^e>>6;
}

///////////////////////////// BH ////////////////////////////////

template <int B> class BH {
  enum {M=8};
  Array<U8, 64> t;
  U32 n;
public:
  BH(int i): t(i*B), n(i-1) {
  }
  U8* operator[](U32 i);
};

template <int B>
inline  U8* BH<B>::operator[](U32 i) {
  int chk=(i>>16^i)&0xffff;
  i=i*M&n;
  U8 *p;
  U16 *cp;
  int j;
  for (j=0; j<M; ++j) {
    p=&t[(i+j)*B];
    cp=(U16*)p;
    if (p[2]==0) *cp=chk;
    if (*cp==chk) break;
  }
  if (j==0) return p+1;
  static U8 tmp[B];
  if (j==M) {
    --j;
    memset(tmp, 0, B);
    tmp[0] = chk;
    tmp[1] = chk>>8;
    if (M>2 && t[(i+j)*B+2]>t[(i+j-1)*B+2]) --j;
  }
  else memcpy(tmp, cp, B);
  memmove(&t[(i+1)*B], &t[i*B], j*B);
  memcpy(&t[i*B], tmp, B);
  return &t[i*B+1];
}

/////////////////////////// ContextMap /////////////////////////

inline int mix2(Mixer& m, int s, StateMap& sm) {
  int p1=sm.p(s);
  int n0=-!nex(s,2);
  int n1=-!nex(s,3);
  int st=stretch(p1)>>2;
  m.add(st);
  p1>>=4;
  int p0=255-p1;
  m.add(p1-p0);
  m.add(st*(n1-n0));
  m.add((p1&n0)-(p0&n1));
  m.add((p1&n1)-(p0&n0));
  return s>0;
}

class RunContextMap {
  BH<4> t;
  U8* cp;
public:
  RunContextMap(int m): t(m/4) {cp=t[0]+1;}
  void set(U32 cx) {  // update count
    if (cp[0]==0 || cp[1]!=buf(1)) cp[0]=1, cp[1]=buf(1);
    else if (cp[0]<255) ++cp[0];
    cp=t[cx]+1;
  }
  int p() {  // predict next bit
    if ((cp[1]+256)>>(8-bpos)==c0)
      return (((cp[1]>>(7-bpos))&1)*2-1)*ilog(cp[0]+1)*8;
    else
      return 0;
  }
  int mix(Mixer& m) {  // return run length
    m.add(p());
    return cp[0]!=0;
  }
};

class SmallStationaryContextMap {
  Array<U16> t;
  int cxt;
  U16 *cp;
public:
  SmallStationaryContextMap(int m): t(m/2), cxt(0) {
    for (U32 i=0; i<t.size(); ++i)
      t[i]=32768;
    cp=&t[0];
  }
  void set(U32 cx) {
    cxt=(cx*256)&(t.size()-256);
  }
  void mix(Mixer& m, int rate=7) {
    *cp += ((y<<16)-*cp+(1<<(rate-1))) >> rate;
    cp=&t[cxt+c0];
    m.add(stretch(*cp>>4));
  }
};

class ContextMap {
  const int C;
  class E {
    U16 chk[7];
    U8 last;
  public:
    U8 bh[7][7];
    U8* get(U16 chk);
  };
  Array<E, 64> t;
  Array<U8*> cp;
  Array<U8*> cp0;
  Array<U32> cxt;
  Array<U8*> runp;
  StateMap *sm;
  int cn;
  void update(U32 cx, int c);
  int mix1(Mixer& m, int cc, int bp, int c1, int y1);
public:
  ContextMap(U32 m, int c=1);
  ~ContextMap();
  void set(U32 cx, int next=-1);
  int mix(Mixer& m) {return mix1(m, c0, bpos, buf(1), y);}
};

inline U8* ContextMap::E::get(U16 ch) {
  if (chk[last&15]==ch) return &bh[last&15][0];
  int b=0xffff, bi=0;
  for (int i=0; i<7; ++i) {
    if (chk[i]==ch) return last=last<<4|i, &bh[i][0];
    int pri=bh[i][0];
    if ((last&15)!=i && last>>4!=i && pri<b) b=pri, bi=i;
  }
  return last=0xf0|bi, chk[bi]=ch, (U8*)memset(&bh[bi][0], 0, 7);
}

ContextMap::ContextMap(U32 m, int c): C(c), t(m>>6), cp(c), cp0(c),
    cxt(c), runp(c), cn(0) {
  sm=new StateMap[C];
  for (int i=0; i<C; ++i) {
    cp0[i]=cp[i]=&t[0].bh[0][0];
    runp[i]=cp[i]+3;
  }
}

ContextMap::~ContextMap() {
  delete[] sm;
}

inline void ContextMap::set(U32 cx, int next) {
  int i=cn++;
  i&=next;
  cx=cx*987654323+i;
  cx=cx<<16|cx>>16;
  cxt[i]=cx*123456791+i;
}

int ContextMap::mix1(Mixer& m, int cc, int bp, int c1, int y1) {
  // Update model with y
  int result=0;
  for (int i=0; i<cn; ++i) {
    if (cp[i]) {
      int ns=nex(*cp[i], y1);
      if (ns>=204 && rnd() << ((452-ns)>>3)) ns-=4;
      *cp[i]=ns;
    }

    // Update context pointers
    if (bpos>1 && runp[i][0]==0)
      cp[i]=0;
    else if (bpos==1||bpos==3||bpos==6)
      cp[i]=cp0[i]+1+(cc&1);
    else if (bpos==4||bpos==7)
      cp[i]=cp0[i]+3+(cc&3);
    else {
      cp0[i]=cp[i]=t[(cxt[i]+cc)&(t.size()-1)].get(cxt[i]>>16);

      // Update pending bit histories for bits 2-7
      if (bpos==0) {
        if (cp0[i][3]==2) {
          const int c=cp0[i][4]+256;
          U8 *p=t[(cxt[i]+(c>>6))&(t.size()-1)].get(cxt[i]>>16);
          p[0]=1+((c>>5)&1);
          p[1+((c>>5)&1)]=1+((c>>4)&1);
          p[3+((c>>4)&3)]=1+((c>>3)&1);
          p=t[(cxt[i]+(c>>3))&(t.size()-1)].get(cxt[i]>>16);
          p[0]=1+((c>>2)&1);
          p[1+((c>>2)&1)]=1+((c>>1)&1);
          p[3+((c>>1)&3)]=1+(c&1);
          cp0[i][6]=0;
        }
        // Update run count of previous context
        if (runp[i][0]==0)
          runp[i][0]=2, runp[i][1]=c1;
        else if (runp[i][1]!=c1)
          runp[i][0]=1, runp[i][1]=c1;
        else if (runp[i][0]<254)
          runp[i][0]+=2;
        else if (runp[i][0]==255)
          runp[i][0]=128;
        runp[i]=cp0[i]+3;
      }
    }

    // predict from last byte in context
    int rc=runp[i][0];
    if ((runp[i][1]+256)>>(8-bp)==cc) {
      int b=((runp[i][1]>>(7-bp))&1)*2-1;
      int c=ilog(rc+1)<<(2+(~rc&1));
      m.add(b*c);
    }
    else
      m.add(0);

    // predict from bit context
    result+=mix2(m, cp[i] ? *cp[i] : 0, sm[i]);
  }
  if (bp==7) cn=0;
  return result;
}

//////////////////////////// matchModel ///////////////////////////

int matchModel(Mixer& m) {
  const int MAXLEN=65534;
  static Array<int> t(MEM);
  static int h=0;
  static int ptr=0;
  static int len=0;
  static int result=0;
  
  static SmallStationaryContextMap scm1(0x20000);

  if (!bpos) {
    h=(h*997*8+buf(1)+1)&(t.size()-1);  // update context hash
    if (len) ++len, ++ptr;
    else {  // find match
      ptr=t[h];
      if (ptr && pos-ptr<(int)buf.size())
        while (buf(len+1)==buf[ptr-len-1] && len<MAXLEN) ++len;
    }
    t[h]=pos;  // update hash table
    result=len;
    scm1.set(pos);
  }

  // predict
  if (len>MAXLEN) len=MAXLEN;
  int sgn;
  if (len && buf(1)==buf[ptr-1] && c0==(buf[ptr]+256)>>(8-bpos)) {
    if ((buf[ptr]>>(7-bpos))&1) sgn=1;
    else sgn=-1;
  }
  else sgn=len=0;
  m.add(sgn*4*ilog(len));
  m.add(sgn*64*min(len, 32));
  scm1.mix(m);
  return result;
}

//////////////////////////// picModel //////////////////////////

void picModel(Mixer& m) {
  static U32 r0, r1, r2, r3;
  static Array<U8> t(0x10200);
  const int N=3;
  static int cxt[N];
  static StateMap sm[N];

  // update the model
  for (int i=0; i<N; ++i)
    t[cxt[i]]=nex(t[cxt[i]],y);

  // update the contexts (pixels surrounding the predicted one)
  r0+=r0+y;
  r1+=r1+((buf(215)>>(7-bpos))&1);
  r2+=r2+((buf(431)>>(7-bpos))&1);
  r3+=r3+((buf(647)>>(7-bpos))&1);
  cxt[0]=(r0&0x7)|((r1>>4)&0x38)|((r2>>3)&0xc0);
  cxt[1]=0x100+((r0&1)|((r1>>4)&0x3e)|((r2>>2)&0x40)|((r3>>1)&0x80));
  cxt[2]=0x200+((r0&0x3f)^(r1&0x3ffe)^((r2<<2)&0x7f00)^((r3<<5)&0xf800));

  // predict
  for (int i=0; i<N; ++i)
    m.add(stretch(sm[i].p(t[cxt[i]])));
}

//////////////////////////// wordModel /////////////////////////
U32 b2=0,b3=0,w4=0;
U32 w5=0,f4=0,tt=0;
U32 WRT_mpw[16]= { 4, 4, 3, 2, 2, 2, 1, 1,  1, 1, 1, 1, 0, 0, 0, 0 };
U32 WRT_mtt[16]= { 0, 0, 1, 2, 3, 4, 5, 5,  6, 6, 6, 6, 7, 7, 7, 7 };
int col=0;

static U32 frstchar=0,spafdo=0,spaces=0,spacecount=0, words=0,wordcount=0,wordlen=0,wordlen1=0;
void wordModel(Mixer& m) {
    static U32 word0=0, word1=0, word2=0, word3=0, word4=0, word5=0;
    static U32 xword0=0,xword1=0,xword2=0,cword0=0,ccword=0;
    static U32 number0=0, number1=0;
    static U32 text0=0;
    static ContextMap cm(MEM*31, 44);
    static int nl1=-3, nl=-2;
    static U32 mask = 0;
    static Array<int> wpos(0x10000);
    static int w=0;
    // Update word hashes
    if (bpos==0) {
        int c=c4&255,f=0;
        if (spaces&0x80000000) --spacecount;
        if (words&0x80000000) --wordcount;
        spaces=spaces*2;
        words=words*2;

        if (c>='A' && c<='Z') c+='a'-'A';
        if ((c>='a' && c<='z') || c==1 || c==2 ||(c>=128 &&(b2!=3))) {
            ++words, ++wordcount;
            word0^=hash(word0, c,0);
            text0=text0*997*16+c;
            wordlen++;
            wordlen=min(wordlen,45);
            f=0;
             w=word0&(wpos.size()-1);
        }
        else {
            if (word0) {
                word5=word4;
                word4=word3;
                word3=word2;
                word2=word1;
                word1=word0;
                wordlen1=wordlen;
                 wpos[w]=blpos;
                if (c==':'|| c=='=') cword0=word0;
                if (c==']'&& (frstchar!=':' || frstchar!='*')) xword0=word0;
                ccword=0;
                word0=wordlen=0;
                if((c=='.'||c=='!'||c=='?' ||c=='}' ||c==')') && buf(2)!=10) f=1; 
                
            }
            if ((c4&0xFFFF)==0x3D3D) xword1=word1,xword2=word2; // ==
            if ((c4&0xFFFF)==0x2727) xword1=word1,xword2=word2; // ''
            if (c==32 || c==10 ) { ++spaces, ++spacecount; if (c==10 ) nl1=nl, nl=pos-1;}
            else if (c=='.' || c=='!' || c=='?' || c==',' || c==';' || c==':') spafdo=0,ccword=c;//*31; 
            else { ++spafdo; spafdo=min(63,spafdo); }
        }
        if (c>='0' && c<='9') {
            number0^=hash(number0, c,1);
        }
        else if (number0) {
            number1=number0;//*11;
            number0=0,ccword=0;
        }
   
        col=min(255, pos-nl);
        int above=buf[nl1+col]; // text column context
        if (col<=2) frstchar=(col==2?min(c,96):0);
        if (frstchar=='[' && c==32)    {if(buf(3)==']' || buf(4)==']' ) frstchar=96,xword0=0;}
        cm.set(hash(513,spafdo, spaces,ccword));
        cm.set(hash(514,frstchar, c));
        cm.set(hash(515,col, frstchar));
        cm.set(hash(516,spaces, (words&255)));
        
        cm.set(hash(256,number0, word2));
        cm.set(hash(257,number0, word1));
        cm.set(hash(258,number1, c,ccword));
        cm.set(hash(259,number0, number1));
        cm.set(hash(260,word0, number1));

       
        cm.set(hash(518,wordlen1,col));
        cm.set(hash(519,c,spacecount/2));
        U32 h=wordcount*64+spacecount;
        cm.set(hash(520,c,h,ccword));
         cm.set(hash(517,frstchar,h));
        cm.set(hash(521,h,spafdo));

        U32 d=c4&0xf0ff;
        cm.set(hash(522,d,frstchar,ccword));

        h=word0*271;
        h=h+buf(1);
        cm.set(hash(262,h, 0));
        cm.set(hash(263,word0, 0));
        cm.set(hash(264,h, word1));
        cm.set(hash(265,word0, word1));
        cm.set(hash(266,h, word1,word2));
        cm.set(hash(268,text0&0xfffff, 0));
          cm.set(hash(269,word0, xword0));
          cm.set(hash(270,word0, xword1));
          cm.set(hash(271,word0, xword2));
          cm.set(hash(272,frstchar, xword2));
        
        cm.set(hash(273,word0, cword0));
        cm.set(hash(274,number0, cword0));
        cm.set(hash(275,h, word2));
        cm.set(hash(276,h, word3));
        cm.set(hash(277,h, word4));
        cm.set(hash(278,h, word5));
        cm.set(hash(279,h, word1,word3));
        cm.set(hash(280,h, word2,word3));
        if (f) {
            word5=word4;//*29;
            word4=word3;//*31;
            word3=word2;//*37;
            word2=word1;//*41;
            word1='.';
        }
        cm.set(hash(523,col,buf(1),above));
        cm.set(hash(524,buf(1),above));
        cm.set(hash(525,col,buf(1)));
        cm.set(hash(526,col,c==32));
        cm.set(hash(281, w, llog(blpos-wpos[w])>>4));
        cm.set(hash(282,buf(1),llog(blpos-wpos[w])>>2));
   
        
   int fl = 0;
    if ((c4&0xff) != 0) {
      if (isalpha(c4&0xff)) fl = 1;
      else if (ispunct(c4&0xff)) fl = 2;
      else if (isspace(c4&0xff)) fl = 3;
      else if ((c4&0xff) == 0xff) fl = 4;
      else if ((c4&0xff) < 16) fl = 5;
      else if ((c4&0xff) < 64) fl = 6;
      else fl = 7;
    }
    mask = (mask<<3)|fl;
 
    cm.set(hash(528,mask,0));
    cm.set(hash(529,mask,buf(1)));
    cm.set(hash(530,mask&0xff,col));
    cm.set(hash(531,mask,buf(2),buf(3)));
    cm.set(hash(532,mask&0x1ff,f4&0x00fff0));
    }
    cm.mix(m);
}

//////////////////////////// nestModel ///////////////////////

void nestModel(Mixer& m)
{
  static int ic=0, bc=0, pc=0,vc=0, qc=0, lvc=0, wc=0;
  static ContextMap cm(MEM/2, 10);
 // static U32 mask = 0;
  if (bpos==0) {
    int c=c4&255, matched=1, vv;
    const int lc = (c >= 'A' && c <= 'Z'?c+'a'-'A':c);
    if (lc == 'a' || lc == 'e' || lc == 'i' || lc == 'o' || lc == 'u') vv = 1; else
    if (lc >= 'a' && lc <= 'z') vv = 2; else
    if (lc == ' ' || lc == '.' || lc == ',' || lc == '!' || lc == '?' || lc == '\n') vv = 3; else
    if (lc >= '0' && lc <= '9') vv = 4; else
    if (lc == 'y') vv = 5; else
    if (lc == '\'') vv = 6; else vv=(c&32)?7:0;
    vc = (vc << 3) | vv;
    if (vv != lvc) {
      wc = (wc << 3) | vv;
      lvc = vv;
    }
    switch(c) {
      case ' ': qc = 0; break;
      case '(': ic += 513; break;
      case ')': ic -= 513; break;
      case '[': ic += 17; break;
      case ']': ic -= 17; break;
      case '<': ic += 23; qc += 34; break;
      case '>': ic -= 23; qc /= 5; break;
      case ':': pc = 20; break;
      case '{': ic += 22; break;
      case '}': ic -= 22; break;
      case '|': pc += 223; break;
      case '"': pc += 0x40; break;
      case '\'': pc += 0x42; break;
      case '\n': pc = qc = 0; break;
      case '.': pc = 0; break;
      case '!': pc = 0; break;
      case '?': pc = 0; break;
      case '#': pc += 0x08; break;
      case '%': pc += 0x76; break;
      case '$': pc += 0x45; break;
      case '*': pc += 0x35; break;
      case '-': pc += 0x3; break;
      case '@': pc += 0x72; break;
      case '&': qc += 0x12; break;
      case ';': qc /= 3; break;
      case '\\': pc += 0x29; break;
      case '/': pc += 0x11;
                if (buf.size() > 1 && buf(1) == '<') qc += 74;
                break;
      case '=': pc += 87; break;
      default: matched = 0;
    }
    if (matched) bc = 0; else bc += 1;
    if (bc > 300) bc = ic = pc = qc = 0;

    cm.set((3*vc+77*pc+373*ic+qc)&0xffff);
    cm.set((31*vc+27*pc+281*qc)&0xffff);
    cm.set((13*vc+271*ic+qc+bc)&0xffff);
    cm.set((17*pc+7*ic)&0xffff);
    cm.set((13*vc+ic)&0xffff);
    cm.set((vc/3+pc)&0xffff);
    cm.set((7*wc+qc)&0xffff);
    cm.set((vc&0xffff)|((f4&0xf)<<16));
    cm.set(((3*pc)&0xffff)|((f4&0xf)<<16));
    cm.set((ic&0xffff)|((f4&0xf)<<16));
  }
  cm.mix(m);
}

//////////////////////////// recordModel ///////////////////////

void recordModel(Mixer& m) {
  static int cpos1[256] , cpos2[256], cpos3[256], cpos4[256];
  static int wpos1[0x10000];
  static int rlen=2, rlen1=3, rlen2=4;
  static int rcount1=0, rcount2=0;
  static ContextMap cm(32768, 3), cn(32768/2, 3), co(32768*2, 3), cp(MEM, 3);

  // Find record length
  if (!bpos) {
    int w=c4&0xffff, c=w&255, d=w>>8;
#if 1
    int r=pos-cpos1[c];
    if (r>1 && r==cpos1[c]-cpos2[c]
        && r==cpos2[c]-cpos3[c] && r==cpos3[c]-cpos4[c]
        && (r>15 || ((c==buf(r*5+1)) && c==buf(r*6+1)))) {
      if (r==rlen1) ++rcount1;
      else if (r==rlen2) ++rcount2;
      else if (rcount1>rcount2) rlen2=r, rcount2=1;
      else rlen1=r, rcount1=1;
    }
    if (rcount1>15 && rlen!=rlen1) rlen=rlen1, rcount1=rcount2=0;
    if (rcount2>15 && rlen!=rlen2) rlen=rlen2, rcount1=rcount2=0;

    // Set 2 dimensional contexts
#endif
    cm.set(c<<8| (min(255, pos-cpos1[c])/4) );
    cm.set(w<<9| llog(pos-wpos1[w])>>2);
    
    cm.set(rlen|buf(rlen)<<10|buf(rlen*2)<<18);
    cn.set(w|rlen<<8);
    cn.set(d|rlen<<16);
    cn.set(c|rlen<<8);

    co.set(buf(1)<<8|min(255, pos-cpos1[buf(1)]));
    co.set(buf(1)<<17|buf(2)<<9|llog(pos-wpos1[w])>>2);
    int col=pos%rlen;
    co.set(buf(1)<<8|buf(rlen));

    cp.set(rlen|buf(rlen)<<10|col<<18);
    cp.set(rlen|buf(1)<<10|col<<18);
    cp.set(col|rlen<<12);

    // update last context positions
    cpos4[c]=cpos3[c];
    cpos3[c]=cpos2[c];
    cpos2[c]=cpos1[c];
    cpos1[c]=pos;
    wpos1[w]=pos;
  }
  cm.mix(m);
  cn.mix(m);
  co.mix(m);
  cp.mix(m);
}

void recordModel1(Mixer& m) {
  static int cpos1[256];
  static int wpos1[0x10000];
  static ContextMap cm(32768, 2), cn(32768/2, 4+1), co(32768*4, 4),cp(32768*2, 3), cq(32768*2, 3);

  // Find record length
  if (!bpos) {
    int w=c4&0xffff, c=w&255, d=w&0xf0ff,e=c4&0xffffff;
   
    cm.set(c<<8| (min(255, pos-cpos1[c])/4));
    cm.set(w<<9| llog(pos-wpos1[w])>>2);
  
    cn.set(w);
    cn.set(d<<8);
    cn.set(c<<16);
    cn.set((f4&0xfffff)); 
    int col=pos&3;
    cn.set(col|2<<12);

    co.set(c    );
    co.set(w<<8 );
    co.set(w5&0x3ffff);
    co.set(e<<3);
    
    cp.set(d    );
    cp.set(c<<8 );
    cp.set(w<<16);

    cq.set(w<<3 );
    cq.set(c<<19);
    cq.set(e);
    // update last context positions

    cpos1[c]=pos;
    wpos1[w]=pos;
  }
  cm.mix(m);
  cn.mix(m);
  co.mix(m);
  cq.mix(m);
  cp.mix(m);
}

//////////////////////////// sparseModel ///////////////////////

void sparseModel(Mixer& m, int seenbefore, int howmany) {
  static ContextMap cm(MEM*2, 40+2);
  if (bpos==0) {
    cm.set(seenbefore);
    cm.set(howmany);
    cm.set(buf(1)|buf(5)<<8);
    cm.set(buf(1)|buf(6)<<8);
    cm.set(buf(3)|buf(6)<<8);
    cm.set(buf(4)|buf(8)<<8);
    cm.set(buf(1)|buf(3)<<8|buf(5)<<16);
    cm.set(buf(2)|buf(4)<<8|buf(6)<<16);
    cm.set(c4&0x00f0f0ff);
    cm.set(c4&0x00ff00ff);
    cm.set(c4&0xff0000ff);
    cm.set(c4&0x00f8f8f8);
    cm.set(c4&0xf8f8f8f8);
    cm.set(f4&0x00000fff);
    cm.set(f4);
    cm.set(c4&0x00e0e0e0);
    cm.set(c4&0xe0e0e0e0);
    cm.set(c4&0x810000c1);
    cm.set(c4&0xC3CCC38C);
    cm.set(c4&0x0081CC81);
    cm.set(c4&0x00c10081);
     for (int i=1; i<8; ++i) {
      cm.set(seenbefore|buf(i)<<8);
      cm.set((buf(i+2)<<8)|buf(i+1));
      cm.set((buf(i+3)<<8)|buf(i+1));
    }
  }
  cm.mix(m);
}

U32 x4=0;
void sparseModel1(Mixer& m, int seenbefore, int howmany) {
   static ContextMap cm(MEM*4, 31);
    static SmallStationaryContextMap scm1(0x10000), scm2(0x20000), scm3(0x2000),
     scm4(0x8000), scm5(0x2000),scm6(0x2000), scma(0x10000);
  if (bpos==0) {
    scm5.set(seenbefore);
    scm6.set(howmany);
  U32  h=x4<<6;
    cm.set(buf(1)+(h&0xffffff00));
    cm.set(buf(1)+(h&0x00ffff00));
    cm.set(buf(1)+(h&0x0000ff00));
      U32 d=c4&0xffff;
     h<<=6;
    cm.set(d+(h&0xffff0000));
    cm.set(d+(h&0x00ff0000));
     h<<=6, d=c4&0xffffff;
    cm.set(d+(h&0xff000000));

    for (int i=1; i<5; ++i) { 
      cm.set(seenbefore|buf(i)<<8);
      cm.set((buf(i+3)<<8)|buf(i+1));
    }
    cm.set(spaces&0x7fff);
    cm.set(spaces&0xff);
    cm.set(words&0x1ffff);
    cm.set(f4&0x000fffff);
    cm.set(tt&0x00000fff);
      h=w4<<6;
    cm.set(buf(1)+(h&0xffffff00));
    cm.set(buf(1)+(h&0x00ffff00));
    cm.set(buf(1)+(h&0x0000ff00));
      d=c4&0xffff;
     h<<=6;
    cm.set(d+(h&0xffff0000));
    cm.set(d+(h&0x00ff0000));
     h<<=6, d=c4&0xffffff;
    cm.set(d+(h&0xff000000));
    cm.set(w4&0xf0f0f0ff);
    
    //cm.set(f4);
    cm.set((w4&63)*128+(5<<17));
    cm.set((f4&0xffff)<<11|frstchar);
    cm.set(spafdo*8*((w4&3)==1));
	
      scm1.set(words&127);
      scm2.set((words&12)*16+(w4&12)*4+(buf(1)>>4));
      scm3.set(w4&15);
      scm4.set(spafdo*((w4&3)==1));
      scma.set(frstchar);
  }
  cm.mix(m);
  scm1.mix(m);
  scm2.mix(m);
  scm3.mix(m);
  scm4.mix(m);
  scm5.mix(m);
  scm6.mix(m);
  scma.mix(m);
}

//////////////////////////// distanceModel ///////////////////////

void distanceModel(Mixer& m) {
  static ContextMap cr(MEM, 3);
  if( bpos == 0 ){
    static int pos00=0,pos20=0,posnl=0;
    int c=c4&0xff;
    if(c==0x00)pos00=pos;
    if(c==0x20)pos20=pos;
    if(c==0xff||c=='\r'||c=='\n')posnl=pos;
    cr.set(min(pos-pos00,255)|(c<<8));
    cr.set(min(pos-pos20,255)|(c<<8));
    cr.set(min(pos-posnl,255)|((c<<8)+234567));
  }
  cr.mix(m);
}

//////////////////////////// bmpModel /////////////////////////////////

inline U32 i4(int i) {
  return buf(i)+256*buf(i-1)+65536*buf(i-2)+16777216*buf(i-3);
}

inline int i2(int i) {
  return buf(i)+256*buf(i-1);
}

inline int sqrbuf(int i) {
  return buf(i)*buf(i);
}

int bmpModel(Mixer& m) {
  static int w=0;
  static int eoi=0;
  static U32 tiff=0;
  const int SC=0x20000;
  static SmallStationaryContextMap scm1(SC), scm2(SC),
    scm3(SC), scm4(SC), scm5(SC), scm6(SC), scm7(SC), scm8(SC), scm9(SC*2),
    scm10(512);
  static ContextMap cm(MEM*4, 15);

  if (!bpos && buf(54)=='B' && buf(53)=='M'
      && i4(44)==54 && i4(40)==40 && i4(24)==0) {
    w=((i4(36)+3)&(-4))*3; 
    const int height=i4(32);
    eoi=pos;
    if (w<0x30000 && height<0x10000) {
      eoi=pos+w*height; 
    }
    else
      eoi=pos;
  }

  if (!bpos) {
    if (c4==0x49492a00) tiff=pos; 
    if (pos-tiff==4 && c4!=0x08000000) tiff=0;
    if (tiff && pos-tiff==200) { 
      int dirsize=i2(pos-tiff-4); 
      w=0;
      int bpp=0, compression=0, width=0, height=0;
      for (int i=tiff+6; i<pos-12 && --dirsize>0; i+=12) {
        int tag=i2(pos-i); 
        int tagfmt=i2(pos-i-2); 
        int taglen=i4(pos-i-4); 
        int tagval=i4(pos-i-8); 
        if ((tagfmt==3||tagfmt==4) && taglen==1) {
          if (tag==256) width=tagval;
          if (tag==257) height=tagval;
          if (tag==259) compression=tagval;
          if (tag==277) bpp=tagval; 
        }
      }
      if (width>0 && height>0 && width*height>50 && compression==1
          && (bpp==1||bpp==3))
        eoi=tiff+width*height*bpp, w=width*bpp;
      if (eoi>pos) {}
      else
        tiff=w=0;
    }
  }
  if (pos>eoi) return w=0;
  if (!bpos) {
    int color=pos%3;
    int mean=buf(3)+buf(w-3)+buf(w)+buf(w+3);
    const int var=(sqrbuf(3)+sqrbuf(w-3)+sqrbuf(w)+sqrbuf(w+3)-(mean*mean/4))>>2;
    mean>>=2;
    const int logvar=ilog(var);
    int i=color<<4;
    cm.set(hash(++i, buf(3)));
    cm.set(hash(++i, buf(3), buf(1)));
    cm.set(hash(++i, buf(3), buf(1), buf(2)));
    cm.set(hash(++i, buf(w)));
    cm.set(hash(++i, buf(w), buf(1)));
    cm.set(hash(++i, buf(w), buf(1), buf(2)));
    cm.set(hash(++i, (buf(3)+buf(w))>>3, buf(1)>>4, buf(2)>>4));
    cm.set(hash(++i, buf(1), buf(2)));
    cm.set(hash(++i, buf(3), buf(1)-buf(4)));
    cm.set(hash(++i, buf(3)+buf(1)-buf(4)));
    cm.set(hash(++i, buf(w), buf(1)-buf(w+1)));
    cm.set(hash(++i, buf(w)+buf(1)-buf(w+1)));
    cm.set(hash(++i, buf(w*3-3), buf(w*3-6)));
    cm.set(hash(++i, buf(w*3+3), buf(w*3+6)));
    cm.set(hash(++i, mean, logvar>>4));
    scm1.set(buf(3)+buf(w)-buf(w+3));
    scm2.set(buf(3)+buf(w-3)-buf(w));
    scm3.set(buf(3)*2-buf(6));
    scm4.set(buf(w)*2-buf(w*2));
    scm5.set(buf(w+3)*2-buf(w*2+6));
    scm6.set(buf(w-3)*2-buf(w*2-6));
    scm7.set(buf(w-3)+buf(1)-buf(w-2));
    scm8.set(buf(w)+buf(w-3)-buf(w*2-3));
    scm9.set(mean>>1|(logvar<<1&0x180));
  }
  scm1.mix(m);
  scm2.mix(m);
  scm3.mix(m);
  scm4.mix(m);
  scm5.mix(m);
  scm6.mix(m);
  scm7.mix(m);
  scm8.mix(m);
  scm9.mix(m);
  scm10.mix(m);
  cm.mix(m);
  return w;
}

//////////////////////////// jpegModel /////////////////////////

class IntBuf {
  Array<int> b;
public:
  IntBuf(int i=0): b(i) {}
  int& operator[](int i) {
    return b[i&(b.size()-1)];
  }
};

#define jassert(x) if (!(x)) { \
  jpeg=0; \
  return 0;}

struct HUF {U32 min, max; int val;};

int jpegModel(Mixer& m) {
  enum {SOF0=0xc0, SOF1, SOF2, SOF3, DHT, RST0=0xd0, SOI=0xd8, EOI, SOS, DQT,
    DNL, DRI, APP0=0xe0, COM=0xfe, FF}; 
  static int jpeg=0; 
  static int app; 
  static int sof=0, sos=0, data=0; 
  static Array<int> ht(8); 
  static int htsize=0; 
  static U32 huffcode=0; 
  static int huffbits=0; 
  static int huffsize=0; 
  static int rs=-1; 
  static int mcupos=0; 
  static Array<HUF> huf(128); 
  static int mcusize=0; 
  static int linesize=0;
  static int hufsel[2][10]; 
  static Array<U8> hbuf(2048); 
  static Array<int> color(10); 
  static Array<int> pred(4); 
  static int dc=0; 
  static int width=0; 
  static int row=0, column=0; 
  static Buf cbuf(0x20000);
  static int cpos=0; 
  static int rs1;
  static int ssum=0, ssum1=0, ssum2=0, ssum3=0;
  static IntBuf cbuf2(0x20000);
  static Array<int> adv_pred(7), sumu(8), sumv(8);
  static Array<int> ls(10); 
  static Array<int> lcp(4), zpos(64);
  static int dqt_state = -1, dqt_end = 0, qnum = 0;
  static Array<U8> qtab(256);
  static Array<int> qmap(10);

  const static U8 zzu[64]={ 
    0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,
    3,2,1,0,1,2,3,4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};
  const static U8 zzv[64]={
    0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,
    4,5,6,7,7,6,5,4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};
  if (!blpos && bpos==1) {
	jpeg=0;
  sof=0, sos=0, data=0; 
 	htsize=0,huffcode=0, huffbits=0,huffsize=0; 
 	rs=-1, mcupos=0, mcusize=0, linesize=0; 
 	dc=0, width=0,row=0, column=0;
 	cpos=0,ssum=0, ssum1=0, ssum2=0, ssum3=0;
 	dqt_state = -1, dqt_end = 0, qnum = 0;
  }

  if (!bpos && !blpos) jpeg=0;
  if (bpos && !jpeg) return 0;
  if (!bpos && app>=0) --app;
  if (app>0) return 0;
  if (!bpos) {
    if (!jpeg && buf(4)==FF && buf(3)==SOI && buf(2)==FF && buf(1)>>4==0xe) {
      jpeg=1;
      sos=sof=htsize=data=mcusize=linesize=0, app=2;
      huffcode=huffbits=huffsize=mcupos=cpos=0, rs=-1;
      memset(&huf[0], 0, huf.size()*sizeof(HUF));
      memset(&pred[0], 0, pred.size()*sizeof(int));
    }
    if (jpeg && data && buf(2)==FF && buf(1) && (buf(1)&0xf8)!=RST0) {
      jassert(buf(1)==EOI);
      jpeg=0;
    }
    if (!jpeg) return 0;
    if (!data && !app && buf(4)==FF && (buf(3)>>4==0xe || buf(3)==COM))
      app=buf(2)*256+buf(1)+2;
    if (buf(5)==FF && buf(4)==SOS) {
      int len=buf(3)*256+buf(2);
      if (len==6+2*buf(1) && buf(1) && buf(1)<=4) 
        sos=pos-5, data=sos+len+2, jpeg=2;
    }
    if (buf(4)==FF && buf(3)==DHT && htsize<8) ht[htsize++]=pos-4;
    if (buf(4)==FF && buf(3)==SOF0) sof=pos-4;
    if (buf(4)==FF && buf(3)==DQT)
      dqt_end=pos+buf(2)*256+buf(1)-1, dqt_state=0;
    else if (dqt_state>=0) {
      if (pos>=dqt_end)
        dqt_state = -1;
      else {
        if (dqt_state%65==0)
          qnum = buf(1);
        else {
          jassert(buf(1)>0);
          jassert(qnum>=0 && qnum<4);
          qtab[qnum*64+((dqt_state%65)-1)]=buf(1)-1;
        }
        dqt_state++;
      }
    }
    if (buf(2)==FF && (buf(1)&0xf8)==RST0) {
      huffcode=huffbits=huffsize=mcupos=0, rs=-1;
      memset(&pred[0], 0, pred.size()*sizeof(int));
    }
  }

  {
    if (pos==data && bpos==1) {
      jassert(htsize>0);
      int i;
      for (i=0; i<htsize; ++i) {
        int p=ht[i]+4; 
        int end=p+buf[p-2]*256+buf[p-1]-2; 
        int count=0; 
        while (p<end && end<pos && end<p+2100 && ++count<10) {
          int tc=buf[p]>>4, th=buf[p]&15;
          if (tc>=2 || th>=4) break;
          jassert(tc>=0 && tc<2 && th>=0 && th<4);
          HUF* h=&huf[tc*64+th*16];
          int val=p+17; 
          int hval=tc*1024+th*256; 
          int j;
          for (j=0; j<256; ++j)
            hbuf[hval+j]=buf[val+j];
          int code=0;
          for (j=0; j<16; ++j) {
            h[j].min=code;
            h[j].max=code+=buf[p+j+1];
            h[j].val=hval;
            val+=buf[p+j+1];
            hval+=buf[p+j+1];
            code*=2;
          }
          p=val;
          jassert(hval>=0 && hval<2048);
        }
        jassert(p==end);
      }
      huffcode=huffbits=huffsize=0, rs=-1;
      if (!sof && sos) return 0;
      int ns=buf[sos+4];
      int nf=buf[sof+9];
      jassert(ns<=4 && nf<=4);
      mcusize=0; 
      int hmax=0; 
      for (i=0; i<ns; ++i) {
        for (int j=0; j<nf; ++j) {
          if (buf[sos+2*i+5]==buf[sof+3*j+10]) {
            int hv=buf[sof+3*j+11]; 
            if (hv>>4>hmax) hmax=hv>>4;
            hv=(hv&15)*(hv>>4); 
            jassert(hv>=1 && hv+mcusize<=10);
            while (hv) {
              jassert(mcusize<10);
              hufsel[0][mcusize]=buf[sos+2*i+6]>>4&15;
              hufsel[1][mcusize]=buf[sos+2*i+6]&15;
              jassert (hufsel[0][mcusize]<4 && hufsel[1][mcusize]<4);
              color[mcusize]=i;
              int tq=buf[sof+3*j+12]; 
              jassert(tq>=0 && tq<4);
              qmap[mcusize]=tq;
              --hv;
              ++mcusize;
            }
          }
        }
      }
      jassert(hmax>=1 && hmax<=10);
      int j;
      for (j=0; j<mcusize; ++j) {
        ls[j]=0;
        for (int i=1; i<mcusize; ++i) if (color[(j+i)%mcusize]==color[j]) ls[j]=i;
        ls[j]=(mcusize-ls[j])<<6;
      }
      for (j=0; j<64; ++j) zpos[zzu[j]+8*zzv[j]]=j;
      width=buf[sof+7]*256+buf[sof+8]; 
      width=(width-1)/(hmax*8)+1; 
      jassert(width>0);
      mcusize*=64; 
      row=column=0;
    }
  }
  {
    if (mcusize && buf(1+(!bpos))!=FF) { 
      jassert(huffbits<=32);
      huffcode+=huffcode+y;
      ++huffbits;
      if (rs<0) {
        jassert(huffbits>=1 && huffbits<=16);
        const int ac=(mcupos&63)>0;
        jassert(mcupos>=0 && (mcupos>>6)<10);
        jassert(ac==0 || ac==1);
        const int sel=hufsel[ac][mcupos>>6];
        jassert(sel>=0 && sel<4);
        const int i=huffbits-1;
        jassert(i>=0 && i<16);
        const HUF *h=&huf[ac*64+sel*16];
        jassert(h[i].min<=h[i].max && h[i].val<2048 && huffbits>0);
        if (huffcode<h[i].max) {
          jassert(huffcode>=h[i].min);
          int k=h[i].val+huffcode-h[i].min;
          jassert(k>=0 && k<2048);
          rs=hbuf[k];
          huffsize=huffbits;
        }
      }
      if (rs>=0) {
        if (huffsize+(rs&15)==huffbits) {
          rs1=rs;
          int x=0; 
          if (mcupos&63) { 
            if (rs==0) {
              mcupos=(mcupos+63)&-64;
              jassert(mcupos>=0 && mcupos<=mcusize && mcupos<=640);
              while (cpos&63) {
                cbuf2[cpos]=0;
                cbuf[cpos++]=0;
              }
            }
            else { 
              jassert((rs&15)<=10);
              const int r=rs>>4;
              const int s=rs&15;
              jassert(mcupos>>6==(mcupos+r)>>6);
              mcupos+=r+1;
              x=huffcode&((1<<s)-1);
              if (s && !(x>>(s-1))) x-=(1<<s)-1;
              for (int i=r; i>=1; --i) {
                cbuf2[cpos]=0;
                cbuf[cpos++]=i<<4|s;
              }
              cbuf2[cpos]=x;
              cbuf[cpos++]=(s<<4)|(huffcode<<2>>s&3)|12;
              ssum+=s;
            }
          }
          else { 
            jassert(rs<12);
            ++mcupos;
            x=huffcode&((1<<rs)-1);
            if (rs && !(x>>(rs-1))) x-=(1<<rs)-1;
            jassert(mcupos>=0 && mcupos>>6<10);
            const int comp=color[mcupos>>6];
            jassert(comp>=0 && comp<4);
            dc=pred[comp]+=x;
            jassert((cpos&63)==0);
            cbuf2[cpos]=dc;
            cbuf[cpos++]=(dc+1023)>>3;
            if ((mcupos>>6)==0) {
              ssum1=0;
              ssum2=ssum3;
            } else {
              if (color[(mcupos>>6)-1]==color[0]) ssum1+=(ssum3=ssum);
              ssum2=ssum1;
            }
            ssum=rs;
          }
          jassert(mcupos>=0 && mcupos<=mcusize);
          if (mcupos>=mcusize) {
            mcupos=0;
            if (++column==width) column=0, ++row;
          }
          huffcode=huffsize=huffbits=0, rs=-1;
          {
            const int acomp=mcupos>>6, q=64*qmap[acomp];
            const int zz=mcupos&63, cpos_dc=cpos-zz;
            if (zz==0) {
              for (int i=0; i<8; ++i) sumu[i]=sumv[i]=0;
              int cpos_dc_ls_acomp = cpos_dc-ls[acomp];
              int cpos_dc_mcusize_width = cpos_dc-mcusize*width;
              for (int i=0; i<64; ++i) {
                sumu[zzu[i]]+=(zzv[i]&1?-1:1)*(zzv[i]?16*(16+zzv[i]):181)*(qtab[q+i]+1)*cbuf2[cpos_dc_mcusize_width+i];
                sumv[zzv[i]]+=(zzu[i]&1?-1:1)*(zzu[i]?16*(16+zzu[i]):181)*(qtab[q+i]+1)*cbuf2[cpos_dc_ls_acomp+i];
              }
            }
            else {
              sumu[zzu[zz-1]]-=(zzv[zz-1]?16*(16+zzv[zz-1]):181)*(qtab[q+zz-1]+1)*cbuf2[cpos-1];
              sumv[zzv[zz-1]]-=(zzu[zz-1]?16*(16+zzu[zz-1]):181)*(qtab[q+zz-1]+1)*cbuf2[cpos-1];
            }

            for (int i=0; i<3; ++i)
              for (int st=0; st<8; ++st) {
                const int zz2=min(zz+st, 63);
                int p=(sumu[zzu[zz2]]*i+sumv[zzv[zz2]]*(2-i))/2;
                p/=(qtab[q+zz2]+1)*181*(16+zzv[zz2])*(16+zzu[zz2])/256;
                if (zz2==0) p-=cbuf2[cpos_dc-ls[acomp]];
                p=(p<0?-1:+1)*ilog(10*abs(p)+1)/10;
                if (st==0) {
                  adv_pred[i]=p;
                  adv_pred[i+4]=p/4;
                }
                else if (abs(p)>abs(adv_pred[i])+1) {
                  adv_pred[i]+=(st*2+(p>0))<<6;
                  if (abs(p/4)>abs(adv_pred[i+4])+1) adv_pred[i+4]+=(st*2+(p>0))<<6;
                  break;
                }
              }
            x=2*sumu[zzu[zz]]+2*sumv[zzv[zz]];
            for (int i=0; i<8; ++i) x-=(zzu[zz]<i)*sumu[i]+(zzv[zz]<i)*sumv[i];
            x/=(qtab[q+zz]+1)*181;
            if (zz==0) x-=cbuf2[cpos_dc-ls[acomp]];
            adv_pred[3]=(x<0?-1:+1)*ilog(10*abs(x)+1)/10;

            for (int i=0; i<4; ++i) {
              const int a=(i&1?zzv[zz]:zzu[zz]), b=(i&2?2:1);
              if (a<b) x=255;
              else {
                const int zz2=zpos[zzu[zz]+8*zzv[zz]-(i&1?8:1)*b];
                x=(qtab[q+zz2]+1)*cbuf2[cpos_dc+zz2]/(qtab[q+zz]+1);
                x=(x<0?-1:+1)*ilog(10*abs(x)+1)/10;
              }
              lcp[i]=x;
            }
            if (column==0) adv_pred[1]=adv_pred[2], adv_pred[0]=1;
            if (row==0) adv_pred[1]=adv_pred[0], adv_pred[2]=1;
          }
        }
      }
    }
  }
  if (!jpeg || !data) return 0;
  if (buf(1+(!bpos))==FF) {
    m.add(128);
    m.set(1, 8);
    m.set(0, 257);
    m.set(buf(1), 256);
    return 1;
  }
  const int N=28;
  static BH<9> t(MEM); 
  static Array<U32> cxt(N); 
  static Array<U8*> cp(N); 
  static StateMap sm[N];
  static Mixer m1(32, 770, 3);
  static APM a1(0x8000), a2(0x10000);
  if (cp[N-1]) {
    for (int i=0; i<N; ++i)
      *cp[i]=nex(*cp[i],y);
  }
  m1.update();
  const int comp=color[mcupos>>6];
  const int coef=(mcupos&63)|comp<<6;
  const int hc=(huffcode*4+((mcupos&63)==0)*2+(comp==0))|1<<(huffbits+2);
  static int hbcount=2;
  if (++hbcount>2 || huffbits==0) hbcount=0;
  jassert(coef>=0 && coef<256);
  const int zu=zzu[mcupos&63], zv=zzv[mcupos&63];
  if (hbcount==0) {
    int n=0;
    cxt[0]=hash(++n, hc, coef, adv_pred[2], ssum2>>6);
    cxt[1]=hash(++n, hc, coef, adv_pred[0], ssum2>>6);
    cxt[2]=hash(++n, hc, coef, adv_pred[1], ssum2>>6);
    cxt[3]=hash(++n, hc, rs1, adv_pred[2]);
    cxt[4]=hash(++n, hc, rs1, adv_pred[0]);
    cxt[5]=hash(++n, hc, rs1, adv_pred[1]);
    cxt[6]=hash(++n, hc, adv_pred[2], adv_pred[0]);
    cxt[7]=hash(++n, hc, cbuf[cpos-width*mcusize], adv_pred[3]);
    cxt[8]=hash(++n, hc, cbuf[cpos-ls[mcupos>>6]], adv_pred[3]);
    cxt[9]=hash(++n, hc, lcp[0], lcp[1], adv_pred[1]);
    cxt[10]=hash(++n, hc, lcp[0], lcp[1], mcupos&63);
    cxt[11]=hash(++n, hc, zu, lcp[0], lcp[2]/3);
    cxt[12]=hash(++n, hc, zv, lcp[1], lcp[3]/3);
    cxt[13]=hash(++n, hc, mcupos>>1);
    cxt[14]=hash(++n, hc, mcupos&63, column>>1);
    cxt[15]=hash(++n, hc, column>>3, lcp[0]+256*(lcp[2]/4), lcp[1]+256*(lcp[3]/4));
    cxt[16]=hash(++n, hc, ssum>>3, mcupos&63);
    cxt[17]=hash(++n, hc, rs1, mcupos&63);
    cxt[18]=hash(++n, hc, mcupos>>3, ssum2>>5, adv_pred[3]);
    cxt[19]=hash(++n, hc, lcp[0]/4, lcp[1]/4, adv_pred[5]);
    cxt[20]=hash(++n, hc, cbuf[cpos-width*mcusize], adv_pred[6]);
    cxt[21]=hash(++n, hc, cbuf[cpos-ls[mcupos>>6]], adv_pred[4]);
    cxt[22]=hash(++n, hc, adv_pred[2]);
    cxt[23]=hash(n, hc, adv_pred[0]);
    cxt[24]=hash(n, hc, adv_pred[1]);
    cxt[25]=hash(++n, hc, zv, lcp[1], adv_pred[6]);
    cxt[26]=hash(++n, hc, zu, lcp[0], adv_pred[4]);
    cxt[27]=hash(++n, hc, lcp[0], lcp[1], adv_pred[3]);
  }
  m1.add(128);
 switch(hbcount)
  {
   case 0: for (int i=0; i<N; ++i) cp[i]=t[cxt[i]]+1, m1.add(stretch(sm[i].p(*cp[i]))); break;
   case 1: { int hc=1+(huffcode&1)*3; for (int i=0; i<N; ++i) cp[i]+=hc, m1.add(stretch(sm[i].p(*cp[i]))); } break;
   default: { int hc=1+(huffcode&1); for (int i=0; i<N; ++i) cp[i]+=hc, m1.add(stretch(sm[i].p(*cp[i]))); } break;
  }

  m1.set(column==0, 2);
  m1.set(coef, 256);
  m1.set(hc&511, 512);
  int pr=m1.p();
  m.add(stretch(pr));
  pr=a1.p(pr, (hc&511)|((adv_pred[1]==0?0:(abs(adv_pred[1])-4)&63)<<9), 1023);
  pr=a2.p(pr, (hc&255)|(coef<<8), 255);
  m.add(stretch(pr));
  m.set(1, 8);
  m.set(1+(hc&255), 257);
  m.set(buf(1), 256);
  return 1;
}

//////////////////////////// exeModel /////////////////////////

U32 execxt(int i, int x=0) {
  int prefix=(buf(i+2)==0x0f)+2*(buf(i+2)==0x66)+3*(buf(i+2)==0x67)
    +4*(buf(i+3)==0x0f)+8*(buf(i+3)==0x66)+12*(buf(i+3)==0x67);
  int opcode=buf(i+1);
  int modrm=i ? buf(i)&0xc7 : 0;
  return prefix|opcode<<4|modrm<<12|x<<20;
}

void exeModel(Mixer& m) {
  const int N=12;
  static ContextMap cm(MEM, N);
  if (!bpos) {
    for (int i=0; i<N; ++i)
      cm.set(execxt(i, buf(1)*(i>4)));
  }
  cm.mix(m);
}

//////////////////////////// indirectModel /////////////////////

void indirectModel(Mixer& m) {
  static ContextMap cm(MEM, 12);
  static U32 t1[256];
  static U16 t2[0x10000];
  static U16 t3[0x8000];
  static U16 t4[0x8000];

  if (!bpos) {
    U32 d=c4&0xffff, c=d&255, d2=(buf(1)&31)+32*(buf(2)&31)+1024*(buf(3)&31);
    U32 d3=(buf(1)>>3&31)+32*(buf(3)>>3&31)+1024*(buf(4)>>3&31);
    U32& r1=t1[d>>8];
    r1=r1<<8|c;
    U16& r2=t2[c4>>8&0xffff];
    r2=r2<<8|c;
    U16& r3=t3[(buf(2)&31)+32*(buf(3)&31)+1024*(buf(4)&31)];
    r3=r3<<8|c;
    U16& r4=t4[(buf(2)>>3&31)+32*(buf(4)>>3&31)+1024*(buf(5)>>3&31)];
    r4=r4<<8|c;
    const U32 t=c|t1[c]<<8;
    const U32 t0=d|t2[d]<<16;
    const U32 ta=d2|t3[d2]<<16;
    const U32 tc=d3|t4[d3]<<16;
    cm.set(t);
    cm.set(t0);
    cm.set(ta);
    cm.set(tc);
    cm.set(t&0xff00);
    cm.set(t0&0xff0000);
    cm.set(ta&0xff0000);
    cm.set(tc&0xff0000);
    cm.set(t&0xffff);
    cm.set(t0&0xffffff);
    cm.set(ta&0xffffff);
    cm.set(tc&0xffffff);
  }
  cm.mix(m);
}

//////////////////////////// dmcModel //////////////////////////

struct DMCNode {
  unsigned int nx[2];
  U8 state;
  unsigned int c0:12, c1:12;
};

void dmcModel(Mixer& m) {
  static int top=0, curr=0;
  static Array<DMCNode> t(MEM*2);
  static StateMap sm;
  static int threshold=256;

  // clone next state
  if (top>0 && top<(int)t.size()) {
    int next=t[curr].nx[y];
    int n=y?t[curr].c1:t[curr].c0;
    int nn=t[next].c0+t[next].c1;
    if (n>=threshold*2 && nn-n>=threshold*3) {
      int r=n*4096/nn;
      t[next].c0 -= t[top].c0 = t[next].c0*r>>12;
      t[next].c1 -= t[top].c1 = t[next].c1*r>>12;
      t[top].nx[0]=t[next].nx[0];
      t[top].nx[1]=t[next].nx[1];
      t[top].state=t[next].state;
      t[curr].nx[y]=top;
      ++top;
      if (top==((int)t.size()*4)/8) {
        threshold=512;
      } else if (top==((int)t.size()*6)/8) {
        threshold=768;
      }
    }
  }

  if (top==(int)t.size() && bpos==1) top=0;
  if (top==0) {
    for (int i=0; i<256; ++i) {
      for (int j=0; j<256; ++j) {
        if (i<127) {
          t[j*256+i].nx[0]=j*256+i*2+1;
          t[j*256+i].nx[1]=j*256+i*2+2;
        }
        else {
          t[j*256+i].nx[0]=(i-127)*256;
          t[j*256+i].nx[1]=(i+1)*256;
        }
        t[j*256+i].c0=128;
        t[j*256+i].c1=128;
      }
    }
    top=65536;
    curr=0;
    threshold=256;
  }

   if (y) {
    if (t[curr].c1<=3840) t[curr].c1+=256;
   } else  if (t[curr].c0<=3840)   t[curr].c0+=256;
  t[curr].state=nex(t[curr].state, y);
  curr=t[curr].nx[y];

  // predict
  const int pr1=sm.p(t[curr].state);
  const int n1=t[curr].c1;
  const int n0=t[curr].c0;
  const int pr2=(n1+5)*4096/(n0+n1+10);
  m.add(stretch(pr1));
  m.add(stretch(pr2));
}

//////////////////////////// contextModel //////////////////////

typedef enum {DEFAULT, JPEG, EXE, TEXT} Filetype;

U32 last_prediction = 2048;

int contextModel2() {
  static ContextMap cm(MEM*31, 9);
  static RunContextMap rcm7(MEM), rcm9(MEM), rcm10(MEM);
  static Mixer m(1200, 10800, 8, 32);
  static U32 cxt[16];  // order 0-11 contexts
  static Filetype filetype=EXE;

  m.update();
  m.add(64);

  // Test for special file types
  int isjpeg=jpegModel(m);  // 1 if JPEG is detected, else 0
  int ismatch=ilog(matchModel(m));  // Length of longest matching context
  int isbmp=bmpModel(m);  // Image width (bytes) if BMP or TIFF detected, or 0

  if (isjpeg) {
    return m.p();
  }
  else if (isbmp>0) {
    static int col=0;
    if (++col>=24) col=0;
    m.set(2, 8);
    m.set(col, 24);
    m.set((buf(isbmp)+buf(3))>>4, 32);
    m.set(c0, 256);
    return m.p();
  }

  // Normal model
  if (bpos==0) {
    for (int i=15; i>0; --i)  // update order 0-11 context hashes
      cxt[i]=cxt[i-1]*257+(c4&255)+1;
    for (int i=0; i<7; ++i)
      cm.set(cxt[i]);
    rcm7.set(cxt[7]);
    cm.set(cxt[8]);
    rcm9.set(cxt[10]);
    rcm10.set(cxt[12]);
    cm.set(cxt[14]);
  }
  int order=cm.mix(m);
  
  rcm7.mix(m);
  rcm9.mix(m);
  rcm10.mix(m);

  if (level>=4) {
    sparseModel(m,ismatch,order);
    sparseModel1(m,ismatch,order);
    distanceModel(m);
    picModel(m);
    recordModel(m);
    recordModel1(m);
    wordModel(m);
    nestModel(m);
    indirectModel(m);
    dmcModel(m);
    if (filetype==EXE) exeModel(m);
  }

  order = order-5;
  if(order<0) order=0;

  U32 d=c0<<(8-bpos),c=(d+(bpos==1?b3/2:0))&192;
  if(!bpos)c=words*16&192;

  U32 c1=buf(1);

  m.set(order*256+(w4&240)+(b2>>4),1536);
  m.set(order*256+(w4&3)*64+(words>>1&63),1536);
  m.set(bpos*256+c1,2048);
  m.set(min(bpos,5)*256+(tt&63)+c,1536);
  m.set(order*256+((d|c1>>bpos)&248)+bpos,1536);
  m.set(bpos*256+(((words<<bpos&255)>>bpos)|(d&255)),2048);
  U32 pr = last_prediction / 16;
  m.set(pr,256);
  m.set(c0,256);

  pr=m.p();
  return pr;
}


//////////////////////////// Predictor /////////////////////////

class Predictor {
  int pr;  // next prediction
public:
  Predictor();
  int p() const {return pr;}
  void update();
};

Predictor::Predictor(): pr(2048) {}

void Predictor::update() {
  static APM a(256), a1(0x10000), a2(0x10000), a3(0x10000),
                      a4(0x10000), a5(0x10000), a6(0x10000);
  static U32 x5=0;
  // Update global context: pos, bpos, c0, c4, buf
  c0+=c0+y;
  if (c0>=256) {
    ++blpos;
    	buf[pos++]=c0;
	c0-=256;
	c4=(c4<<8)+c0;
        int i=WRT_mpw[c0>>4];
		w4=w4*4+i;
		if (b2==3) i=2;
		w5=w5*4+i;
		b3=b2;
        b2=c0;   
	    x4=x4*256+c0,x5=(x5<<8)+c0;
	    if(c0=='.' || c0=='!' || c0=='?' || c0=='/'|| c0==')') {
			w5=(w5<<8)|0x3ff,f4=(f4&0xfffffff0)+2,x5=(x5<<8)+c0,x4=x4*256+c0;
            if(c0!='!') w4|=12, tt=(tt&0xfffffff8)+1,b3='.';
		}
		if (c0==32) --c0;
		tt=tt*8+WRT_mtt[c0>>4];
		f4=f4*16+(c0>>4);
	c0=1;
  }
  bpos=(bpos+1)&7;

  // Filter the context model with APMs
  int pr0=contextModel2();

  pr=a.p(pr0, c0);
  
  int pr1=a1.p(pr0, c0+256*buf(1));
  int pr2=a2.p(pr0, c0^(hash(buf(1), buf(2))&0xffff));
  int pr3=a3.p(pr0, c0^(hash(buf(1), buf(2), buf(3))&0xffff));
  pr0=(pr0+pr1+pr2+pr3+2)>>2;
  
      pr1=a4.p(pr, c0+256*buf(1));
      pr2=a5.p(pr, c0^(hash(buf(1), buf(2))&0xffff));
      pr3=a6.p(pr, c0^(hash(buf(1), buf(2), buf(3))&0xffff));
  pr=(pr+pr1+pr2+pr3+2)>>2;

  pr=(pr+pr0+1)>>1;
  last_prediction = pr;
}

Predictor paq8;

}  // namespace

PAQ8L::PAQ8L(int memory) {
  level = memory;
  buf.setsize(MEM*8);
}

float PAQ8L::Predict() {
  return paq8.p() * conversion_factor;
}

void PAQ8L::Perceive(int bit) {
  y = bit;
  paq8.update();
}

const std::vector<float>& PAQ8L::ModelPredictions() {
  return model_predictions;
}
