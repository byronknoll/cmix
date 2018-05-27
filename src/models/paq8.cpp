// This code is a hybrid of paq8l, paq8pxd (released by Kaido Orav) and paq8px (v101 and up, released by Jan Ondrus and Márcio Pais).

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

#include "paq8.h"
#include "../preprocess/preprocessor.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#define NDEBUG

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
typedef unsigned char U8;
typedef unsigned short U16;
typedef unsigned int U32;
typedef unsigned long long U64;

#ifndef min
inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}
#endif

using preprocessor::Filetype;

void quit(const char* message=0) {}

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
  const size_t sz=ALIGN+n*sizeof(T);
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

int level=DEFAULT_OPTION;
unsigned long long MEM() {
  return 0x10000UL<<level;
}

int y=0;

int c0=1;
U32 c4=0;
int bpos=0;
int blpos=0;
Buf buf;

int dt[1024];  // i -> 16K/(i+3)

struct ModelStats{
  U32 XML, x86_64, Record;
};

#define CacheSize (1<<5)

#if (CacheSize&(CacheSize-1)) || (CacheSize<8)
  #error Cache size must be a power of 2 bigger than 4
#endif

inline unsigned BitCount(unsigned v){
  v -= ((v>>1)&0x55555555);
  v = ((v>>2)&0x33333333) + (v&0x33333333);
  v = ((v>>4)+v)&0x0f0f0f0f;
  v = ((v>>8)+v)&0x00ff00ff;
  v = ((v>>16)+v)&0x0000ffff;
  return v;
}

inline unsigned ilog2(unsigned x) {
  x = x | (x >> 1);
  x = x | (x >> 2);
  x = x | (x >> 4);
  x = x | (x >> 8);
  x = x | (x >>16);
  return BitCount(x >> 1);
}

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

static const U8 State_table[256][4]={
  { 1, 2, 0, 0},{ 3, 5, 1, 0},{ 4, 6, 0, 1},{ 7, 10, 2, 0},
  { 8, 12, 1, 1},{ 9, 13, 1, 1},{ 11, 14, 0, 2},{ 15, 19, 3, 0},
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
  {140,252, 0,41}};

#define nex(state,sel) State_table[state][sel]

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
  for (int x=-2047; x<=2047; ++x) {
    int i=squash(x);
    for (int j=pi; j<=i; ++j)
      t[j]=x;
    pi=i+1;
  }
  t[4095]=2047;
}

#if !defined(__GNUC__)

#if (2 == _M_IX86_FP)
#define __SSE2__
#elif (1 == _M_IX86_FP)
#define __SSE__
#endif

#endif

static int dot_product (const short* const t, const short* const w, int n);

static void train (const short* const t, short* const w, int n, const int e);

#if defined(__SSE2__)
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

#elif defined(__SSE__)
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
#endif

#define NUM_INPUTS 1607
#define NUM_SETS 22

std::valarray<float> model_predictions(0.5, NUM_INPUTS + NUM_SETS + 11);
unsigned int prediction_index = 0;
float conversion_factor = 1.0 / 4095;

void AddPrediction(int x) {
  model_predictions[prediction_index++] = x * conversion_factor;
}

void ResetPredictions() {
  prediction_index = 0;
}

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

  void update() {
    for (int i=0; i<ncxt; ++i) {
      int err=((y<<12)-pr[i])*7;
      train(&tx[0], &wx[cxt[i]*N], nx, err);
    }
    nx=base=ncxt=0;
  }

  void add(int x) {
    AddPrediction(squash(x));
    tx[nx++]=x;
  }

  void set(int cx, int range) {
    cxt[ncxt++]=base+cx;
    base+=range;
  }

  int p() {
    while (nx&7) tx[nx++]=0;
    if (mp) {
      mp->update();
      for (int i=0; i<ncxt; ++i) {
        pr[i]=squash((dot_product(&tx[0], &wx[cxt[i]*N], nx) * 9)>>9);
        mp->add(stretch(pr[i]));
      }
      mp->set(0, 1);
      return mp->p();
    }
    else {
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

class APM {
  int index;
  const int N;
  Array<U16> t;
public:
  APM(int n);
  int p(int pr=2048, int cxt=0, int rate=7) {
    pr=stretch(pr);
    int g=(y<<16)+(y<<rate)-y-y;
    t[index] += (g-t[index]) >> rate;
    t[index+1] += (g-t[index+1]) >> rate;
    const int w=pr&127;
    index=((pr+2048)>>7)+cxt*33;
    return (t[index]*(128-w)+t[index+1]*w) >> 11;
  }
};

APM::APM(int n): index(0), N(n), t(n*33) {
  for (int i=0; i<N; ++i)
    for (int j=0; j<33; ++j)
      t[i*33+j] = i==0 ? squash((j-16)*128)*16 : t[j];
}

class StateMap {
protected:
  int cxt;
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

class StateMap32 {
protected:
  const int N;  // Number of contexts
  int cxt;      // Context of last prediction
  Array<U32> t;       // cxt -> prediction in high 22 bits, count in low 10 bits
  inline void update(int limit) {
    U32 *p=&t[cxt], p0=p[0];
    int n=p0&1023;  //count
    int pr=p0>>10;  //prediction
    if (n<limit) ++p0;
    else p0=(p0&0xfffffc00)|limit;
    int target=y<<22;
    int delta=((target-pr)>>3)*dt[n]; //the larger the count (n) the less it should adapt
    p0+=delta&0xfffffc00; 
    p[0]=p0;
  }

public:
  StateMap32(int n=256);
  void Reset(int Rate=0){
    for (int i=0; i<N; ++i)
      t[i]=(t[i]&0xfffffc00)|min(Rate, t[i]&0x3FF);
  }
  // update bit y (0..1), predict next bit in context cx
  int p(int cx, int limit=1023) {
    update(limit);
    return t[cxt=cx]>>20;
  }
};

StateMap32::StateMap32(int n): N(n), cxt(0), t(n) {
  for (int i=0; i<N; ++i)
    t[i]=(1u<<31)+0;  //initial p=0.5, initial count=0
}

inline U32 hash(U32 a, U32 b, U32 c=0xffffffff, U32 d=0xffffffff,
    U32 e=0xffffffff) {
  U32 h=a*200002979u+b*30005491u+c*50004239u+d*70004807u+e*110002499u;
  return h^h>>9^a>>2^b>>3^c>>4^d>>5^e>>6;
}

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
inline U8* BH<B>::operator[](U32 i) {
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

//////////////////////////// HashTable /////////////////////////

// A HashTable maps a 32-bit index to an array of B bytes.
// The first byte is a checksum using the upper 8 bits of the
// index.  The second byte is a priority (0 = empty) for hash
// replacement.  The index need not be a hash.

// HashTable<B> h(n) - create using n bytes  n and B must be
//     powers of 2 with n >= B*4, and B >= 2.
// h[i] returns array [1..B-1] of bytes indexed by i, creating and
//     replacing another element if needed.  Element 0 is the
//     checksum and should not be modified.

template <int B>
class HashTable {
  Array<U8,64> t;  // table: 1 element = B bytes: checksum priority data data
  const int N;  // size in bytes
public:
  HashTable(int n): t(n), N(n) {
  }
  U8* operator[](U32 i);
};

template <int B>
inline U8* HashTable<B>::operator[](U32 i) {
  i*=123456791;
  i=i<<16|i>>16;
  i*=234567891;
  int chk=i>>24;
  i=i*B&(N-B);
  U8 *p = &t[0];
  if (p[i]==chk) return p+i;
  if (p[i^B]==chk) return p+(i^B);
  if (p[i^B*2]==chk) return p+(i^B*2);
  if (p[i+1]>p[(i+1)^B] || p[i+1]>p[(i+1)^B*2]) i^=B;
  if (p[i+1]>p[(i+1)^B^B*2]) i^=B^B*2;
  memset(p+i, 0, B);
  p[i]=chk;
  return p+i;
}

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
  void set(U32 cx) {
    if (cp[0]==0 || cp[1]!=buf(1)) cp[0]=1, cp[1]=buf(1);
    else if (cp[0]<255) ++cp[0];
    cp=t[cx]+1;
  }
  int p() {
    if ((cp[1]+256)>>(8-bpos)==c0)
      return (((cp[1]>>(7-bpos))&1)*2-1)*ilog(cp[0]+1)*8;
    else
      return 0;
  }
  int mix(Mixer& m) {
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
    m.add(stretch(*cp>>4)/2);
  }
};

/*
  Map for modelling contexts of (nearly-)stationary data.
  The context is looked up directly. For each bit modelled, a 32bit element stores
  a 22 bit prediction and a 10 bit adaptation rate offset.

  - BitsOfContext: How many bits to use for each context. Higher bits are discarded.
  - BitsPerContext: How many bits [1..8] of input are to be modelled for each context.
    New contexts must be set at those intervals.
  - Rate: Initial adaptation rate offset [0..1023]. Lower offsets mean faster adaptation.
    Will be increased on every occurrence until the higher bound is reached.

  Uses 2^(BitsOfContext+BitsPerContext+2) bytes of memory.
*/

class StationaryMap {
  Array<U32> Data;
  int Context, Mask, bCount, B;
  U32 *cp;
public:
  StationaryMap(int BitsOfContext, int BitsPerContext, int Rate=0): Data(1<<(BitsOfContext+BitsPerContext)), Context(0), Mask(1<<BitsPerContext), bCount(1), B(1) {
    Reset(Rate);
    cp=&Data[0];
  }
  void set(U32 ctx) {
    Context = (ctx*Mask)&(Data.size()-Mask);
  }
  void Reset( int Rate = 0 ){
    for (unsigned int i=0; i<Data.size(); ++i)
      Data[i]=(0x7FF<<20)|min(0x3FF,Rate);
  }
  void mix(Mixer& m) {
    U32 Count = min(0x3FF, ((*cp)&0x3FF)+1);
    int Prediction = (*cp)>>10, Error = (y<<22)-Prediction;
    Error = ((Error/8)*dt[Count])/1024;
    Prediction = min(0x3FFFFF,max(0,Prediction+Error));
    *cp = (Prediction<<10)|Count;
    B+=(y && B>1);
    cp=&Data[Context+B];
    Prediction = (*cp)>>20;
    m.add(stretch(Prediction)/2);
    m.add((Prediction-2048)/4);
    bCount<<=1; B<<=1;
    if (bCount==Mask){
      bCount=1;
      B=1;
    }
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

  int result=0;
  for (int i=0; i<cn; ++i) {
    if (cp[i]) {
      int ns=nex(*cp[i], y1);
      if (ns>=204 && rnd() << ((452-ns)>>3)) ns-=4;
      *cp[i]=ns;
    }

    if (bpos>1 && runp[i][0]==0)
      cp[i]=0;
    else if (bpos==1||bpos==3||bpos==6)
      cp[i]=cp0[i]+1+(cc&1);
    else if (bpos==4||bpos==7)
      cp[i]=cp0[i]+3+(cc&3);
    else {
      cp0[i]=cp[i]=t[(cxt[i]+cc)&(t.size()-1)].get(cxt[i]>>16);

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

    int rc=runp[i][0];
    if ((runp[i][1]+256)>>(8-bp)==cc) {
      int b=((runp[i][1]>>(7-bp))&1)*2-1;
      int c=ilog(rc+1)<<(2+(~rc&1));
      m.add(b*c);
    }
    else
      m.add(0);

    result+=mix2(m, cp[i] ? *cp[i] : 0, sm[i]);
  }
  if (bp==7) cn=0;
  return result;
}

/*
Context map for large contexts (32bits).
Maps to a bit history state, a 3 MRU byte history, and 1 byte RunStats.

Bit and byte histories are stored in a hash table with 64 byte buckets. 
The buckets are indexed by a context ending after 0, 2 or 5 bits of the
current byte. Thus, each byte modeled results in 3 main memory accesses
per context, with all other accesses to cache.

On a byte boundary (bit 0), only 3 of the 7 bit history states are used.
Of the remaining 4 bytes, 3 are then used to store the last bytes seen
in this context, 7 bits to store the length of consecutive occurrences of
the previously seen byte, and 1 bit to signal if more than 1 byte as been
seen in this context. The byte history is then combined with the bit history
states to provide additional states that are then mapped to predictions.
*/

class ContextMap2 {
  const U32 C; // max number of contexts
  class Bucket { // hash bucket, 64 bytes
    U16 Checksums[7]; // byte context checksums
    U8 MRU; // last 2 accesses (0-6) in low, high nibble
  public:
    U8 BitState[7][7]; // byte context, 3-bit context -> bit history state
                       // BitState[][0] = 1st bit, BitState[][1,2] = 2nd bit, BitState[][3..6] = 3rd bit
                       // BitState[][0] is also a replacement priority, 0 = empty
    inline U8* Find(U16 Checksum) { // Find or create hash element matching checksum. 
                                    // If not found, insert or replace lowest priority (skipping 2 most recent).
      if (Checksums[MRU&15]==Checksum)
        return &BitState[MRU&15][0];
      int worst=0xFFFF, index=0;
      for (int i=0; i<7; ++i) {
        if (Checksums[i]==Checksum)
          return MRU=MRU<<4|i, (U8*)&BitState[i][0];
        if (BitState[i][0]<worst && (MRU&15)!=i && MRU>>4!=i) {
          worst = BitState[i][0];
          index=i;
        }
      }
      MRU = 0xF0|index;
      Checksums[index] = Checksum;
      return (U8*)memset(&BitState[index][0], 0, 7);
    }
  };
  Array<Bucket, 64> Table; // bit histories for bits 0-1, 2-4, 5-7
                           // For 0-1, also contains run stats in BitState[][3] and byte history in BitState[][4..6]
  Array<U8*> BitState; // C pointers to current bit history states
  Array<U8*> BitState0; // First element of 7 element array containing BitState[i]
  Array<U8*> ByteHistory; // C pointers to run stats plus byte history, 4 bytes, [RunStats,1..3]
  Array<U32> Contexts; // C whole byte contexts (hashes)
  Array<bool> HasHistory; // True if context has a full valid byte history (i.e., seen at least 3 times)
  StateMap32 **Maps6b, **Maps8b, **Maps12b;
  U32 index; // Next context to set by set()
  U32 bits;
  U8 lastByte, lastBit, bitPos;
  inline void Update() {
    U64 mask = Table.size()-1;
    for (U32 i=0; i<index; i++) {
      if (BitState[i])
        *BitState[i] = nex(*BitState[i], lastBit);

      if (bitPos>1 && ByteHistory[i][0]==0)
        BitState[i] = nullptr;
      else {
        U16 chksum = Contexts[i]>>16;
        switch (bitPos) {
          case 0: {
            BitState[i] = BitState0[i] = Table[(Contexts[i]+bits)&mask].Find(chksum);
            // Update pending bit histories for bits 2-7
            if (BitState0[i][3]==2) {
              const int c = BitState0[i][4]+256;
              U8 *p = Table[(Contexts[i]+(c>>6))&mask].Find(chksum);
              p[0] = 1+((c>>5)&1);
              p[1+((c>>5)&1)] = 1+((c>>4)&1);
              p[3+((c>>4)&3)] = 1+((c>>3)&1);
              p = Table[(Contexts[i]+(c>>3))&mask].Find(chksum);
              p[0] = 1+((c>>2)&1);
              p[1+((c>>2)&1)] = 1+((c>>1)&1);
              p[3+((c>>1)&3)] = 1+(c&1);
              BitState0[i][6] = 0;
            }
            // Update byte history of previous context
            ByteHistory[i][3] = ByteHistory[i][2];
            ByteHistory[i][2] = ByteHistory[i][1];
            if (ByteHistory[i][0]==0)  // new context
              ByteHistory[i][0]=2, ByteHistory[i][1]=lastByte;
            else if (ByteHistory[i][1]!=lastByte)  // different byte in context
              ByteHistory[i][0]=1, ByteHistory[i][1]=lastByte;
            else if (ByteHistory[i][0]<254)  // same byte in context
              ByteHistory[i][0]+=2;
            else if (ByteHistory[i][0]==255) // more than one byte seen, but long run of current byte, reset to single byte seen
              ByteHistory[i][0] = 128;

            ByteHistory[i] = BitState0[i]+3;
            HasHistory[i] = *BitState0[i]>15;
            break;
          }
          case 2: case 5: {
            BitState[i] = BitState0[i] = Table[(Contexts[i]+bits)&mask].Find(chksum);
            break;
          }
          case 1: case 3: case 6: BitState[i] = BitState0[i]+1+lastBit; break;
          case 4: case 7: BitState[i] = BitState0[i]+3+(bits&3); break;
        }
      }
    }
  }
public:
  // Construct using Size bytes of memory for Count contexts
  ContextMap2(const U64 Size, const U32 Count) : C(Count), Table(Size>>6), BitState(Count), BitState0(Count), ByteHistory(Count), Contexts(Count), HasHistory(Count){
    Maps6b = new StateMap32*[C];
    Maps8b = new StateMap32*[C];
    Maps12b = new StateMap32*[C];
    for (U32 i=0; i<C; i++) {
      Maps6b[i] = new StateMap32((1<<6)+8);
      Maps8b[i] = new StateMap32(1<<8);
      Maps12b[i] = new StateMap32((1<<12)+(1<<9));
      BitState[i] = BitState0[i] = &Table[i].BitState[0][0];
      ByteHistory[i] = BitState[i]+3;
    }
    index = 0;
    lastByte = lastBit = 0;
    bits = 1;  bitPos = 0;
  }
  ~ContextMap2() {
    for (U32 i=0; i<C; i++) {
      delete Maps6b[i];
      delete Maps8b[i];
      delete Maps12b[i];
    }
    delete[] Maps6b;
    delete[] Maps8b;
    delete[] Maps12b;
  }
  inline void set(U32 ctx) { // set next whole byte context to ctx
    ctx = ctx*987654323+index; // permute (don't hash) ctx to spread the distribution
    ctx = ctx<<16|ctx>>16;
    Contexts[index] = ctx*123456791+index;
    index++;
  }
  void Train(const U8 B){
    for (bitPos=0; bitPos<8; bitPos++){
      Update();
      lastBit = (B>>(7-bitPos))&1;
      bits += bits+lastBit;
    }
    index = 0;
    bits = 1; bitPos = 0;
    lastByte = B;
  }
  int mix(Mixer& m) {
    int result = 0;
    lastBit = y;
    bitPos = bpos;
    bits+=bits+lastBit;
    lastByte = bits&0xFF;
    if (bitPos==0)
      bits = 1;
    Update();

    for (U32 i=0; i<index; i++) {
      // predict from bit context
      int state = (BitState[i])?*BitState[i]:0;
      result+=(state>0);
      int p1 = Maps8b[i]->p(state);
      int n0=nex(state, 2), n1=nex(state, 3), k=-~n1;
      k = (k*64)/(k-~n0);
      n0=-!n0, n1=-!n1;
      // predict from last byte in context
      if ((U32)((ByteHistory[i][1]+256)>>(8-bitPos))==bits){
        int RunStats = ByteHistory[i][0]; // count*2, +1 if 2 different bytes seen
        int sign=(ByteHistory[i][1]>>(7-bitPos)&1)*2-1;  // predicted bit + for 1, - for 0
        int value = ilog(RunStats+1)<<(3-(RunStats&1));
        m.add(sign*value);
      }
      else if (bitPos>0 && (ByteHistory[i][0]&1)>0) {
        if ((U32)((ByteHistory[i][2]+256)>>(8-bitPos))==bits)
          m.add((((ByteHistory[i][2]>>(7-bitPos))&1)*2-1)*128);
        else if (HasHistory[i] && (U32)((ByteHistory[i][3]+256)>>(8-bitPos))==bits)
          m.add((((ByteHistory[i][3]>>(7-bitPos))&1)*2-1)*128);
        else
          m.add(0);
      }
      else
        m.add(0);

      if (HasHistory[i]) {
        state  = (ByteHistory[i][1]>>(7-bitPos))&1;
        state |= ((ByteHistory[i][2]>>(7-bitPos))&1)*2;
        state |= ((ByteHistory[i][3]>>(7-bitPos))&1)*4;
      }
      else
        state = 8;

      int st = stretch(p1)>>2;
      m.add(st);
      m.add((p1-2047)>>3);
      p1>>=4;
      int p0 = 255-p1;
      m.add((st*(n1-n0)));
      m.add((p1&n0)-(p0&n1));
      m.add((p1&n1)-(p0&n0));
      m.add(stretch(Maps12b[i]->p((state<<9)|(bitPos<<6)|k))>>2);
      m.add(stretch(Maps6b[i]->p((state<<3)|bitPos))>>2);
    }
    if (bitPos==7) index = 0;
    return result;
  }
};

//////////////////////////// Text modelling /////////////////////////

#define TAB 0x09
#define NEW_LINE 0x0A
#define CARRIAGE_RETURN 0x0D
#define SPACE 0x20

inline bool CharInArray(const char c, const char a[], const int len) {
  if (a==nullptr)
    return false;
  int i=0;
  for (; i<len && c!=a[i]; i++);
  return i<len;
}

#define MAX_WORD_SIZE 64

class Word {
public:
  U8 Letters[MAX_WORD_SIZE];
  U8 Start, End;
  U32 Hash[4], Type, Language;
  Word() : Start(0), End(0), Hash{0,0,0,0}, Type(0) {
    memset(&Letters[0], 0, sizeof(U8)*MAX_WORD_SIZE);
  }
  bool operator==(const char *s) const {
    size_t len=strlen(s);
    return ((size_t)(End-Start+(Letters[Start]!=0))==len && memcmp(&Letters[Start], s, len)==0);
  }
  bool operator!=(const char *s) const {
    return !operator==(s);
  }
  void operator+=(const char c) {
    if (End<MAX_WORD_SIZE-1) {
      End+=(Letters[End]>0);
      Letters[End]=tolower(c);
    }
  }
  U8 operator[](U8 i) const {
    return (End-Start>=i)?Letters[Start+i]:0;
  }
  U8 operator()(U8 i) const {
    return (End-Start>=i)?Letters[End-i]:0;
  }
  U32 Length() const {
    if (Letters[Start]!=0)
      return End-Start+1;
    return 0;
  }
  void GetHashes() {
    Hash[0] = 0xc01dflu, Hash[1] = ~Hash[0];
    for (int i=Start; i<=End; i++) {
      U8 l = Letters[i];
      Hash[0]^=hash(Hash[0], l, i);
      Hash[1]^=hash(Hash[1], 
        ((l&0x80)==0)?l&0x5F:
        ((l&0xC0)==0x80)?l&0x3F:
        ((l&0xE0)==0xC0)?l&0x1F:
        ((l&0xF0)==0xE0)?l&0xF:l&0x7
      );
    }
    Hash[2] = (~Hash[0])^Hash[1];
    Hash[3] = (~Hash[1])^Hash[0];
  }
  bool ChangeSuffix(const char *OldSuffix, const char *NewSuffix) {
    size_t len=strlen(OldSuffix);
    if (Length()>len && memcmp(&Letters[End-len+1], OldSuffix, len)==0) {
      size_t n=strlen(NewSuffix);
      if (n>0) {
        memcpy(&Letters[End-int(len)+1], NewSuffix, min(MAX_WORD_SIZE-1,End+int(n))-End);
        End=min(MAX_WORD_SIZE-1, End-int(len)+int(n));
      }
      else
        End-=U8(len);
      return true;
    }
    return false;
  }
  bool MatchesAny(const char* a[], const int count) {
    int i=0;
    size_t len = (size_t)Length();
    for (; i<count && (len!=strlen(a[i]) || memcmp(&Letters[Start], a[i], len)!=0); i++);
    return i<count;
  }
  bool EndsWith(const char *Suffix) const {
    size_t len=strlen(Suffix);
    return (Length()>len && memcmp(&Letters[End-len+1], Suffix, len)==0);
  }
  bool StartsWith(const char *Prefix) const {
    size_t len=strlen(Prefix);
    return (Length()>len && memcmp(&Letters[Start], Prefix, len)==0);
  }
};

class Segment {
public:
  Word FirstWord; // useful following questions
  U32 WordCount;
  U32 NumCount;
};

class Sentence : public Segment {
public:
  enum Types { // possible sentence types, excluding Imperative
    Declarative,
    Interrogative,
    Exclamative,
    Count
  };
  Types Type;
  U32 SegmentCount;
  U32 VerbIndex; // relative position of last detected verb
  U32 NounIndex; // relative position of last detected noun
  U32 CapitalIndex; // relative position of last capitalized word, excluding the initial word of this sentence
  Word lastVerb, lastNoun, lastCapital;
};

class Paragraph {
public:
  U32 SentenceCount, TypeCount[Sentence::Types::Count], TypeMask;
};

class Language {
public:
  enum Flags {
    Verb                   = (1<<0),
    Noun                   = (1<<1)
  };
  enum Ids {
    Unknown,
    English,
    French,
    German,
    Count
  };
  virtual ~Language() {};
  virtual bool IsAbbreviation(Word *W) = 0;
};

class English: public Language {
private:
  static const int NUM_ABBREV = 6;
  const char *Abbreviations[NUM_ABBREV]={ "mr","mrs","ms","dr","st","jr" };
public:
  enum Flags {
    Adjective              = (1<<2),
    Plural                 = (1<<3),
    Male                   = (1<<4),
    Female                 = (1<<5),
    Negation               = (1<<6),
    PastTense              = (1<<7)|Verb,
    PresentParticiple      = (1<<8)|Verb,
    AdjectiveSuperlative   = (1<<9)|Adjective,
    AdjectiveWithout       = (1<<10)|Adjective,
    AdjectiveFull          = (1<<11)|Adjective,
    AdverbOfManner         = (1<<12),
    SuffixNESS             = (1<<13),
    SuffixITY              = (1<<14)|Noun,
    SuffixCapable          = (1<<15),
    SuffixNCE              = (1<<16),
    SuffixNT               = (1<<17),
    SuffixION              = (1<<18),
    SuffixAL               = (1<<19)|Adjective,
    SuffixIC               = (1<<20)|Adjective,
    SuffixIVE              = (1<<21),
    SuffixOUS              = (1<<22)|Adjective,
    PrefixOver             = (1<<23),
    PrefixUnder            = (1<<24)
  };
  bool IsAbbreviation(Word *W) { return W->MatchesAny(Abbreviations, NUM_ABBREV); };
};

class French: public Language {
private:
  static const int NUM_ABBREV = 2;
  const char *Abbreviations[NUM_ABBREV]={ "m","mm" };
public:
  enum Flags {
    Adjective              = (1<<2),
    Plural                 = (1<<3)
  };
  bool IsAbbreviation(Word *W) { return W->MatchesAny(Abbreviations, NUM_ABBREV); };
};

class German : public Language {
private:
  static const int NUM_ABBREV = 3;
  const char *Abbreviations[NUM_ABBREV]={ "fr","hr","hrn" };
public:
  enum Flags {
    Adjective              = (1<<2),
    Plural                 = (1<<3),
    Female                 = (1<<4)
  };
  bool IsAbbreviation(Word *W) { return W->MatchesAny(Abbreviations, NUM_ABBREV); };
};

//////////////////////////// Stemming routines /////////////////////////

class Stemmer {
protected:
  U32 GetRegion(const Word *W, const U32 From) {
    bool hasVowel = false;
    for (int i=W->Start+From; i<=W->End; i++) {
      if (IsVowel(W->Letters[i])) {
        hasVowel = true;
        continue;
      }
      else if (hasVowel)
        return i-W->Start+1;
    }
    return W->Start+W->Length();
  }
  bool SuffixInRn(const Word *W, const U32 Rn, const char *Suffix) {
    return (W->Start!=W->End && Rn<=W->Length()-strlen(Suffix));
  }
public:
  virtual ~Stemmer() {};
  virtual bool IsVowel(const char c) = 0;
  virtual void Hash(Word *W) = 0;
  virtual bool Stem(Word *W) = 0;
};

/*
  English affix stemmer, based on the Porter2 stemmer.

  Changelog:
  (29/12/2017) v127: Initial release by Márcio Pais
  (02/01/2018) v128: Small changes to allow for compilation with MSVC
  Fix buffer overflow (thank you Jan Ondrus)
  (28/01/2018) v133: Refactoring, added processing of "non/non-" prefixes
  (04/02/2018) v135: Refactoring, added processing of gender-specific words and common words
*/

class EnglishStemmer: public Stemmer {
private:
  static const int NUM_VOWELS = 6;
  const char Vowels[NUM_VOWELS]={'a','e','i','o','u','y'};
  static const int NUM_DOUBLES = 9;
  const char Doubles[NUM_DOUBLES]={'b','d','f','g','m','n','p','r','t'};
  static const int NUM_LI_ENDINGS = 10;
  const char LiEndings[NUM_LI_ENDINGS]={'c','d','e','g','h','k','m','n','r','t'};
  static const int NUM_NON_SHORT_CONSONANTS = 3;
  const char NonShortConsonants[NUM_NON_SHORT_CONSONANTS]={'w','x','Y'};
  static const int NUM_MALE_WORDS = 9;
  const char *MaleWords[NUM_MALE_WORDS]={"he","him","his","himself","man","men","boy","husband","actor"};
  static const int NUM_FEMALE_WORDS = 8;
  const char *FemaleWords[NUM_FEMALE_WORDS]={"she","her","herself","woman","women","girl","wife","actress"};
  static const int NUM_COMMON_WORDS = 12;
  const char *CommonWords[NUM_COMMON_WORDS]={"the","be","to","of","and","in","that","you","have","with","from","but"};
  static const int NUM_SUFFIXES_STEP0 = 3;
  const char *SuffixesStep0[NUM_SUFFIXES_STEP0]={"'s'","'s","'"};
  static const int NUM_SUFFIXES_STEP1b = 6;
  const char *SuffixesStep1b[NUM_SUFFIXES_STEP1b]={"eedly","eed","ed","edly","ing","ingly"};
  const U32 TypesStep1b[NUM_SUFFIXES_STEP1b]={English::AdverbOfManner,0,English::PastTense,English::AdverbOfManner|English::PastTense,English::PresentParticiple,English::AdverbOfManner|English::PresentParticiple};
  static const int NUM_SUFFIXES_STEP2 = 22;
  const char *(SuffixesStep2[NUM_SUFFIXES_STEP2])[2]={
    {"ization", "ize"},
    {"ational", "ate"},
    {"ousness", "ous"},
    {"iveness", "ive"},
    {"fulness", "ful"},
    {"tional", "tion"},
    {"lessli", "less"},
    {"biliti", "ble"},
    {"entli", "ent"},
    {"ation", "ate"},
    {"alism", "al"},
    {"aliti", "al"},
    {"fulli", "ful"},
    {"ousli", "ous"},
    {"iviti", "ive"},
    {"enci", "ence"},
    {"anci", "ance"},
    {"abli", "able"},
    {"izer", "ize"},
    {"ator", "ate"},
    {"alli", "al"},
    {"bli", "ble"}
  };
  const U32 TypesStep2[NUM_SUFFIXES_STEP2]={
    English::SuffixION,
    English::SuffixION|English::SuffixAL,
    English::SuffixNESS,
    English::SuffixNESS,
    English::SuffixNESS,
    English::SuffixION|English::SuffixAL,
    English::AdverbOfManner,
    English::AdverbOfManner|English::SuffixITY,
    English::AdverbOfManner,
    English::SuffixION,
    0,
    English::SuffixITY,
    English::AdverbOfManner,
    English::AdverbOfManner,
    English::SuffixITY,
    0,
    0,
    English::AdverbOfManner,
    0,
    0,
    English::AdverbOfManner,
    English::AdverbOfManner
  };
  static const int NUM_SUFFIXES_STEP3 = 8;
  const char *(SuffixesStep3[NUM_SUFFIXES_STEP3])[2]={
    {"ational", "ate"},
    {"tional", "tion"},
    {"alize", "al"},
    {"icate", "ic"},
    {"iciti", "ic"},
    {"ical", "ic"},
    {"ful", ""},
    {"ness", ""}
  };
  const U32 TypesStep3[NUM_SUFFIXES_STEP3]={English::SuffixION|English::SuffixAL,English::SuffixION|English::SuffixAL,0,0,English::SuffixITY,English::SuffixAL,English::AdjectiveFull,English::SuffixNESS};
  static const int NUM_SUFFIXES_STEP4 = 20;
  const char *SuffixesStep4[NUM_SUFFIXES_STEP4]={"al","ance","ence","er","ic","able","ible","ant","ement","ment","ent","ou","ism","ate","iti","ous","ive","ize","sion","tion"};
  const U32 TypesStep4[NUM_SUFFIXES_STEP4]={
    English::SuffixAL,
    English::SuffixNCE,
    English::SuffixNCE,
    0,
    English::SuffixIC,
    English::SuffixCapable,
    English::SuffixCapable,
    English::SuffixNT,
    0,
    0,
    English::SuffixNT,
    0,
    0,
    0,
    English::SuffixITY,
    English::SuffixOUS,
    English::SuffixIVE,
    0,
    English::SuffixION,
    English::SuffixION
  };
  static const int NUM_EXCEPTION_REGION1 = 3;
  const char *ExceptionsRegion1[NUM_EXCEPTION_REGION1]={"gener","arsen","commun"};
  static const int NUM_EXCEPTIONS1 = 18;
  const char *(Exceptions1[NUM_EXCEPTIONS1])[2]={
    {"skis", "ski"},
    {"skies", "sky"},
    {"dying", "die"},
    {"lying", "lie"},
    {"tying", "tie"},
    {"idly", "idle"},
    {"gently", "gentle"},
    {"ugly", "ugli"},
    {"early", "earli"},
    {"only", "onli"},
    {"singly", "singl"},
    {"sky", "sky"},
    {"news", "news"},
    {"howe", "howe"},
    {"atlas", "atlas"},
    {"cosmos", "cosmos"},
    {"bias", "bias"},
    {"andes", "andes"}
  };
  const U32 TypesExceptions1[NUM_EXCEPTIONS1]={
    English::Noun|English::Plural,
    English::Plural,
    English::PresentParticiple,
    English::PresentParticiple,
    English::PresentParticiple,
    English::AdverbOfManner,
    English::AdverbOfManner,
    English::Adjective,
    English::Adjective|English::AdverbOfManner,
    0,
    English::AdverbOfManner,
    English::Noun,
    English::Noun,
    0,
    English::Noun,
    English::Noun,
    English::Noun,
    0
  };
  static const int NUM_EXCEPTIONS2 = 8;
  const char *Exceptions2[NUM_EXCEPTIONS2]={"inning","outing","canning","herring","earring","proceed","exceed","succeed"};
  const U32 TypesExceptions2[NUM_EXCEPTIONS2]={English::Noun,English::Noun,English::Noun,English::Noun,English::Noun,English::Verb,English::Verb,English::Verb}; 
  inline bool IsConsonant(const char c) {
    return !IsVowel(c);
  }
  inline bool IsShortConsonant(const char c) {
    return !CharInArray(c, NonShortConsonants, NUM_NON_SHORT_CONSONANTS);
  }
  inline bool IsDouble(const char c) {
    return CharInArray(c, Doubles, NUM_DOUBLES);
  }
  inline bool IsLiEnding(const char c) {
    return CharInArray(c, LiEndings, NUM_LI_ENDINGS);
  }
  U32 GetRegion1(const Word *W) {
    for (int i=0; i<NUM_EXCEPTION_REGION1; i++) {
      if (W->StartsWith(ExceptionsRegion1[i]))
        return U32(strlen(ExceptionsRegion1[i]));
    }
    return GetRegion(W, 0);
  }
  bool EndsInShortSyllable(const Word *W) {
    if (W->End==W->Start)
      return false;
    else if (W->End==W->Start+1)
      return IsVowel((*W)(1)) && IsConsonant((*W)(0));
    else
      return (IsConsonant((*W)(2)) && IsVowel((*W)(1)) && IsConsonant((*W)(0)) && IsShortConsonant((*W)(0)));
  }
  bool IsShortWord(const Word *W) {
    return (EndsInShortSyllable(W) && GetRegion1(W)==W->Length());
  }
  inline bool HasVowels(const Word *W) {
    for (int i=W->Start; i<=W->End; i++) {
      if (IsVowel(W->Letters[i]))
        return true;
    }
    return false;
  }
  bool TrimStartingApostrophe(Word *W) {
    bool r=(W->Start!=W->End && (*W)[0]=='\'');
    W->Start+=(U8)r;
    return r;
  }
  void MarkYsAsConsonants(Word *W) {
    if ((*W)[0]=='y')
      W->Letters[W->Start]='Y';
    for (int i=W->Start+1; i<=W->End; i++) {
      if (IsVowel(W->Letters[i-1]) && W->Letters[i]=='y')
        W->Letters[i]='Y';
    }
  }
  bool ProcessPrefixes(Word *W) {
    if (W->StartsWith("irr") && W->Length()>5 && ((*W)[3]=='a' || (*W)[3]=='e'))
      W->Start+=2, W->Type|=English::Negation;
    else if (W->StartsWith("over") && W->Length()>5)
      W->Start+=4, W->Type|=English::PrefixOver;
    else if (W->StartsWith("under") && W->Length()>6)
      W->Start+=5, W->Type|=English::PrefixUnder;
    else if (W->StartsWith("unn") && W->Length()>5)
      W->Start+=2, W->Type|=English::Negation;
    else if (W->StartsWith("non") && W->Length()>(U32)(5+((*W)[3]=='-')))
      W->Start+=2+((*W)[3]=='-'), W->Type|=English::Negation;
    else
      return false;
    return true;
  }
  bool ProcessSuperlatives(Word *W) {
    if (W->EndsWith("est") && W->Length()>4) {
      U8 i=W->End;
      W->End-=3;
      W->Type|=English::AdjectiveSuperlative;

      if ((*W)(0)==(*W)(1) && (*W)(0)!='r' && !(W->Length()>=4 && memcmp("sugg",&W->Letters[W->End-3],4)==0)) {
        W->End-= ( ((*W)(0)!='f' && (*W)(0)!='l' && (*W)(0)!='s') ||
                   (W->Length()>4 && (*W)(1)=='l' && ((*W)(2)=='u' || (*W)(3)=='u' || (*W)(3)=='v'))) &&
                   (!(W->Length()==3 && (*W)(1)=='d' && (*W)(2)=='o'));
        if (W->Length()==2 && ((*W)[0]!='i' || (*W)[1]!='n'))
          W->End = i, W->Type&=~English::AdjectiveSuperlative;
      }
      else {
        switch((*W)(0)) {
          case 'd': case 'k': case 'm': case 'y': break;
          case 'g': {
            if (!( W->Length()>3 && ((*W)(1)=='n' || (*W)(1)=='r') && memcmp("cong",&W->Letters[W->End-3],4)!=0 ))
              W->End = i, W->Type&=~English::AdjectiveSuperlative;
            else
              W->End+=((*W)(2)=='a');
            break;
          }
          case 'i': { W->Letters[W->End]='y'; break; }
          case 'l': {
            if (W->End==W->Start+1 || memcmp("mo",&W->Letters[W->End-2],2)==0)
              W->End = i, W->Type&=~English::AdjectiveSuperlative;
            else
              W->End+=IsConsonant((*W)(1));
            break;
          }
          case 'n': {
            if (W->Length()<3 || IsConsonant((*W)(1)) || IsConsonant((*W)(2)))
              W->End = i, W->Type&=~English::AdjectiveSuperlative;
            break;
          }
          case 'r': {
            if (W->Length()>3 && IsVowel((*W)(1)) && IsVowel((*W)(2)))
              W->End+=((*W)(2)=='u') && ((*W)(1)=='a' || (*W)(1)=='i');
            else
              W->End = i, W->Type&=~English::AdjectiveSuperlative;
            break;
          }
          case 's': { W->End++; break; }
          case 'w': {
            if (!(W->Length()>2 && IsVowel((*W)(1))))
              W->End = i, W->Type&=~English::AdjectiveSuperlative;
            break;
          }
          case 'h': {
            if (!(W->Length()>2 && IsConsonant((*W)(1))))
              W->End = i, W->Type&=~English::AdjectiveSuperlative;
            break;
          }
          default: {
            W->End+=3;
            W->Type&=~English::AdjectiveSuperlative;
          }
        }
      }
    }
    return (W->Type&English::AdjectiveSuperlative)>0;
  }
  bool Step0(Word *W) {
    for (int i=0; i<NUM_SUFFIXES_STEP0; i++) {
      if (W->EndsWith(SuffixesStep0[i])) {
        W->End-=U8(strlen(SuffixesStep0[i]));
        W->Type|=English::Plural;
        return true;
      }
    }
    return false;
  }
  bool Step1a(Word *W) {
    if (W->EndsWith("sses")) {
      W->End-=2;
      W->Type|=English::Plural;
      return true;
    }
    if (W->EndsWith("ied") || W->EndsWith("ies")) {
      W->Type|=((*W)(0)=='d')?English::PastTense:English::Plural;
      W->End-=1+(W->Length()>4);
      return true;
    }
    if (W->EndsWith("us") || W->EndsWith("ss"))
      return false;
    if ((*W)(0)=='s' && W->Length()>2) {
      for (int i=W->Start;i<=W->End-2;i++) {
        if (IsVowel(W->Letters[i])) {
          W->End--;
          W->Type|=English::Plural;
          return true;
        }
      }
    }
    if (W->EndsWith("n't") && W->Length()>4) {
      switch ((*W)(3)) {
        case 'a': {
          if ((*W)(4)=='c')
            W->End-=2;
          else
            W->ChangeSuffix("n't","ll");
          break;
        }
        case 'i': { W->ChangeSuffix("in't","m"); break; }
        case 'o': {
          if ((*W)(4)=='w')
            W->ChangeSuffix("on't","ill");
          else
            W->End-=3;
          break;
        }
        default: W->End-=3;
      }
      W->Type|=English::Negation;
      return true;
    }
    if (W->EndsWith("hood") && W->Length()>7) {
      W->End-=4;
      return true;
    }
    return false;
  }
  bool Step1b(Word *W, const U32 R1) {
    for (int i=0; i<NUM_SUFFIXES_STEP1b; i++) {
      if (W->EndsWith(SuffixesStep1b[i])) {
        switch(i) {
          case 0: case 1: {
            if (SuffixInRn(W, R1, SuffixesStep1b[i]))
              W->End-=1+i*2;
            break;
          }
          default: {
            U8 j=W->End;
            W->End-=U8(strlen(SuffixesStep1b[i]));
            if (HasVowels(W)) {
              if (W->EndsWith("at") || W->EndsWith("bl") || W->EndsWith("iz") || IsShortWord(W))
                (*W)+='e';
              else if (W->Length()>2) {
                if ((*W)(0)==(*W)(1) && IsDouble((*W)(0)))
                  W->End--;
                else if (i==2 || i==3) {
                  switch((*W)(0)) {
                    case 'c': case 's': case 'v': { W->End+=!(W->EndsWith("ss") || W->EndsWith("ias")); break; }
                    case 'd': {
                      static const char nAllowed[4] = {'a','e','i','o'};
                      W->End+=IsVowel((*W)(1)) && (!CharInArray((*W)(2), nAllowed, 4)); break;
                    }
                    case 'k': { W->End+=W->EndsWith("uak"); break; }
                    case 'l': {
                      static const char Allowed1[10] = {'b','c','d','f','g','k','p','t','y','z'};
                      static const char Allowed2[4] = {'a','i','o','u'};
                      W->End+= CharInArray((*W)(1), Allowed1, 10) ||
                                (CharInArray((*W)(1), Allowed2, 4) && IsConsonant((*W)(2)));
                      break;
                    }
                  }
                }
                else if (i>=4) {
                  switch((*W)(0)) {
                    case 'd': {
                      if (IsVowel((*W)(1)) && (*W)(2)!='a' && (*W)(2)!='e' && (*W)(2)!='o')
                        (*W)+='e';
                      break;
                    }
                    case 'g': {
                      static const char Allowed[7] = {'a','d','e','i','l','r','u'};
                      if (
                        CharInArray((*W)(1), Allowed, 7) || (
                         (*W)(1)=='n' && (
                          (*W)(2)=='e' ||
                          ((*W)(2)=='u' && (*W)(3)!='b' && (*W)(3)!='d') ||
                          ((*W)(2)=='a' && ((*W)(3)=='r' || ((*W)(3)=='h' && (*W)(4)=='c'))) ||
                          (W->EndsWith("ring") && ((*W)(4)=='c' || (*W)(4)=='f'))
                         )
                        ) 
                      )
                        (*W)+='e';
                      break;
                    }
                    case 'l': {
                      if (!((*W)(1)=='l' || (*W)(1)=='r' || (*W)(1)=='w' || (IsVowel((*W)(1)) && IsVowel((*W)(2)))))
                        (*W)+='e';
                      if (W->EndsWith("uell") && W->Length()>4 && (*W)(4)!='q')
                        W->End--;
                      break;
                    }
                    case 'r': {
                      if ((
                        ((*W)(1)=='i' && (*W)(2)!='a' && (*W)(2)!='e' && (*W)(2)!='o') ||
                        ((*W)(1)=='a' && (!((*W)(2)=='e' || (*W)(2)=='o' || ((*W)(2)=='l' && (*W)(3)=='l')))) ||
                        ((*W)(1)=='o' && (!((*W)(2)=='o' || ((*W)(2)=='t' && (*W)(3)!='s')))) ||
                        (*W)(1)=='c' || (*W)(1)=='t') && (!W->EndsWith("str"))
                      )
                        (*W)+='e';
                      break;
                    }
                    case 't': {
                      if ((*W)(1)=='o' && (*W)(2)!='g' && (*W)(2)!='l' && (*W)(2)!='i' && (*W)(2)!='o')
                        (*W)+='e';
                      break;
                    }
                    case 'u': {
                      if (!(W->Length()>3 && IsVowel((*W)(1)) && IsVowel((*W)(2))))
                        (*W)+='e';
                      break;
                    }
                    case 'z': {
                      if (W->EndsWith("izz") && W->Length()>3 && ((*W)(3)=='h' || (*W)(3)=='u'))
                        W->End--;
                      else if ((*W)(1)!='t' && (*W)(1)!='z')
                        (*W)+='e';
                      break;
                    }
                    case 'k': {
                      if (W->EndsWith("uak"))
                        (*W)+='e';
                      break;
                    }
                    case 'b': case 'c': case 's': case 'v': {
                      if (!(
                        ((*W)(0)=='b' && ((*W)(1)=='m' || (*W)(1)=='r')) ||
                        W->EndsWith("ss") || W->EndsWith("ias") || (*W)=="zinc"
                      ))
                        (*W)+='e';
                      break;
                    }
                  }
                }
              }
            }
            else {
              W->End=j;
              return false;
            }
          }
        }
        W->Type|=TypesStep1b[i];
        return true;
      }
    }
    return false;
  }
  bool Step1c(Word *W) {
    if (W->Length()>2 && tolower((*W)(0))=='y' && IsConsonant((*W)(1))) {
      W->Letters[W->End]='i';
      return true;
    }
    return false;
  }
  bool Step2(Word *W, const U32 R1) {
    for (int i=0; i<NUM_SUFFIXES_STEP2; i++) {
      if (W->EndsWith(SuffixesStep2[i][0]) && SuffixInRn(W, R1, SuffixesStep2[i][0])) {
        W->ChangeSuffix(SuffixesStep2[i][0], SuffixesStep2[i][1]);
        W->Type|=TypesStep2[i];
        return true;
      }
    }
    if (W->EndsWith("logi") && SuffixInRn(W, R1, "ogi")) {
      W->End--;
      return true;
    }
    else if (W->EndsWith("li")) {
      if (SuffixInRn(W, R1, "li") && IsLiEnding((*W)(2))) {
        W->End-=2;
        W->Type|=English::AdverbOfManner;
        return true;
      }
      else if (W->Length()>3) {
        switch((*W)(2)) {
            case 'b': {
              W->Letters[W->End]='e';
              W->Type|=English::AdverbOfManner;
              return true;              
            }
            case 'i': {
              if (W->Length()>4) {
                W->End-=2;
                W->Type|=English::AdverbOfManner;
                return true;
              }
              break;
            }
            case 'l': {
              if (W->Length()>5 && ((*W)(3)=='a' || (*W)(3)=='u')) {
                W->End-=2;
                W->Type|=English::AdverbOfManner;
                return true;
              }
              break;
            }
            case 's': {
              W->End-=2;
              W->Type|=English::AdverbOfManner;
              return true;
            }
            case 'e': case 'g': case 'm': case 'n': case 'r': case 'w': {
              if (W->Length()>(U32)(4+((*W)(2)=='r'))) {
                W->End-=2;
                W->Type|=English::AdverbOfManner;
                return true;
              }
            }
        }
      }
    }
    return false;
  }
  bool Step3(Word *W, const U32 R1, const U32 R2) {
    bool res=false;
    for (int i=0; i<NUM_SUFFIXES_STEP3; i++) {
      if (W->EndsWith(SuffixesStep3[i][0]) && SuffixInRn(W, R1, SuffixesStep3[i][0])) {
        W->ChangeSuffix(SuffixesStep3[i][0], SuffixesStep3[i][1]);
        W->Type|=TypesStep3[i];
        res=true;
        break;
      }
    }
    if (W->EndsWith("ative") && SuffixInRn(W, R2, "ative")) {
      W->End-=5;
      W->Type|=English::SuffixIVE;
      return true;
    }
    if (W->Length()>5 && W->EndsWith("less")) {
      W->End-=4;
      W->Type|=English::AdjectiveWithout;
      return true;
    }
    return res;
  }
  bool Step4(Word *W, const U32 R2) {
    bool res=false;
    for (int i=0; i<NUM_SUFFIXES_STEP4; i++) {
      if (W->EndsWith(SuffixesStep4[i]) && SuffixInRn(W, R2, SuffixesStep4[i])) {
        W->End-=U8(strlen(SuffixesStep4[i])-(i>17));
        if (i!=10 || (*W)(0)!='m')
          W->Type|=TypesStep4[i];
        if (i==0 && W->EndsWith("nti")) {
          W->End--;
          res=true;
          continue;
        }
        return true;
      }
    }
    return res;
  }
  bool Step5(Word *W, const U32 R1, const U32 R2) {
    if ((*W)(0)=='e' && (*W)!="here") {
      if (SuffixInRn(W, R2, "e"))
        W->End--;
      else if (SuffixInRn(W, R1, "e")) {
        W->End--;
        W->End+=EndsInShortSyllable(W);
      }
      else
        return false;
      return true;
    }
    else if (W->Length()>1 && (*W)(0)=='l' && SuffixInRn(W, R2, "l") && (*W)(1)=='l') {
      W->End--;
      return true;
    }
    return false;
  }
public:
  inline bool IsVowel(const char c) final {
    return CharInArray(c, Vowels, NUM_VOWELS);
  }
  inline void Hash(Word *W) final {
    W->Hash[2] = W->Hash[3] = 0xb0a710ad;
    for (int i=W->Start; i<=W->End; i++) {
      U8 l = W->Letters[i];
      W->Hash[2]=W->Hash[2]*263*32 + l;
      if (IsVowel(l))
        W->Hash[3]=W->Hash[3]*997*8 + (l/4-22);
      else if (l>='b' && l<='z')
        W->Hash[3]=W->Hash[3]*271*32 + (l-97);
      else
        W->Hash[3]=W->Hash[3]*11*32 + l;
    }
  }
  bool Stem(Word *W) {
    if (W->Length()<2) {
      Hash(W);
      return false;
    }
    bool res = TrimStartingApostrophe(W);
    res|=ProcessPrefixes(W);
    res|=ProcessSuperlatives(W);
    for (int i=0; i<NUM_EXCEPTIONS1; i++) {
      if ((*W)==Exceptions1[i][0]) {
        if (i<11) {
          size_t len=strlen(Exceptions1[i][1]);
          memcpy(&W->Letters[W->Start], Exceptions1[i][1], len);
          W->End=W->Start+U8(len-1);
        }
        Hash(W);
        W->Type|=TypesExceptions1[i];
        W->Language = Language::English;
        return (i<11);
      }
    }

    // Start of modified Porter2 Stemmer
    MarkYsAsConsonants(W);
    U32 R1=GetRegion1(W), R2=GetRegion(W,R1);
    res|=Step0(W);
    res|=Step1a(W);
    for (int i=0; i<NUM_EXCEPTIONS2; i++) {
      if ((*W)==Exceptions2[i]) {
        Hash(W);
        W->Type|=TypesExceptions2[i];
        W->Language = Language::English;
        return res;
      }
    }
    res|=Step1b(W, R1);
    res|=Step1c(W);
    res|=Step2(W, R1);
    res|=Step3(W, R1, R2);
    res|=Step4(W, R2);
    res|=Step5(W, R1, R2);

    for (U8 i=W->Start; i<=W->End; i++) {
      if (W->Letters[i]=='Y')
        W->Letters[i]='y';
    }
    if (!W->Type || W->Type==English::Plural) {
      if (W->MatchesAny(MaleWords, NUM_MALE_WORDS))
        res = true, W->Type|=English::Male;
      else if (W->MatchesAny(FemaleWords, NUM_FEMALE_WORDS))
        res = true, W->Type|=English::Female;
    }
    if (!res)
      res=W->MatchesAny(CommonWords, NUM_COMMON_WORDS);
    Hash(W);
    if (res)
      W->Language = Language::English;
    return res;
  }
};

/*
  French suffix stemmer, based on the Porter stemmer.

  Changelog:
  (28/01/2018) v133: Initial release by Márcio Pais
  (04/02/2018) v135: Added processing of common words
  (25/02/2018) v139: Added UTF8 conversion
*/

class FrenchStemmer: public Stemmer {
private:
  static const int NUM_VOWELS = 17;
  const char Vowels[NUM_VOWELS]={'a','e','i','o','u','y','\xE2','\xE0','\xEB','\xE9','\xEA','\xE8','\xEF','\xEE','\xF4','\xFB','\xF9'};
  static const int NUM_COMMON_WORDS = 10;
  const char *CommonWords[NUM_COMMON_WORDS]={"de","la","le","et","en","un","une","du","que","pas"};
  static const int NUM_EXCEPTIONS = 3;
  const char *(Exceptions[NUM_EXCEPTIONS])[2]={
    {"monument", "monument"},
    {"yeux", "oeil"},
    {"travaux", "travail"},
  };
  const U32 TypesExceptions[NUM_EXCEPTIONS]={
    French::Noun,
    French::Noun|French::Plural,
    French::Noun|French::Plural
  };
  static const int NUM_SUFFIXES_STEP1 = 39;
  const char *SuffixesStep1[NUM_SUFFIXES_STEP1]={
    "ance","iqUe","isme","able","iste","eux","ances","iqUes","ismes","ables","istes", //11
    "atrice","ateur","ation","atrices","ateurs","ations", //6
    "logie","logies", //2
    "usion","ution","usions","utions", //4
    "ence","ences", //2
    "issement","issements", //2
    "ement","ements", //2
    "it\xE9","it\xE9s", //2
    "if","ive","ifs","ives", //4
    "euse","euses", //2
    "ment","ments" //2
  };
  static const int NUM_SUFFIXES_STEP2a = 35;
  const char *SuffixesStep2a[NUM_SUFFIXES_STEP2a]={
    "issaIent", "issantes", "iraIent", "issante",
    "issants", "issions", "irions", "issais",
    "issait", "issant", "issent", "issiez", "issons",
    "irais", "irait", "irent", "iriez", "irons",
    "iront", "isses", "issez", "\xEEmes",
    "\xEEtes", "irai", "iras", "irez", "isse",
    "ies", "ira", "\xEEt", "ie", "ir", "is",
    "it", "i"
  };
  static const int NUM_SUFFIXES_STEP2b = 38;
  const char *SuffixesStep2b[NUM_SUFFIXES_STEP2b]={
    "eraIent", "assions", "erions", "assent",
    "assiez", "\xE8rent", "erais", "erait",
    "eriez", "erons", "eront", "aIent", "antes",
    "asses", "ions", "erai", "eras", "erez",
    "\xE2mes", "\xE2tes", "ante", "ants",
    "asse", "\xE9""es", "era", "iez", "ais",
    "ait", "ant", "\xE9""e", "\xE9s", "er",
    "ez", "\xE2t", "ai", "as", "\xE9", "a"
  };
  static const int NUM_SET_STEP4 = 6;
  const char SetStep4[NUM_SET_STEP4]={'a','i','o','u','\xE8','s'};
  static const int NUM_SUFFIXES_STEP4 = 7;
  const char *SuffixesStep4[NUM_SUFFIXES_STEP4]={"i\xE8re","I\xE8re","ion","ier","Ier","e","\xEB"};
  static const int NUM_SUFFIXES_STEP5 = 5;
  const char *SuffixesStep5[NUM_SUFFIXES_STEP5]={"enn","onn","ett","ell","eill"};
  inline bool IsConsonant(const char c) {
    return !IsVowel(c);
  }
  void ConvertUTF8(Word *W) {
    for (int i=W->Start; i<W->End; i++) {
      U8 c = W->Letters[i+1]+((W->Letters[i+1]<0xA0)?0x60:0x40);
      if (W->Letters[i]==0xC3 && (IsVowel(c) || (W->Letters[i+1]&0xDF)==0x87)) {
        W->Letters[i] = c;
        if (i+1<W->End)
          memcpy(&W->Letters[i+1], &W->Letters[i+2], W->End-i-1);
        W->End--;
      }
    }
  }
  void MarkVowelsAsConsonants(Word *W) {
    for (int i=W->Start; i<=W->End; i++) {
      switch (W->Letters[i]) {
        case 'i': case 'u': {
          if (i>W->Start && i<W->End && (IsVowel(W->Letters[i-1]) || (W->Letters[i-1]=='q' && W->Letters[i]=='u')) && IsVowel(W->Letters[i+1]))
            W->Letters[i] = toupper(W->Letters[i]);
          break;
        }
        case 'y': {
          if ((i>W->Start && IsVowel(W->Letters[i-1])) || (i<W->End && IsVowel(W->Letters[i+1])))
            W->Letters[i] = toupper(W->Letters[i]);
        }
      }
    }
  }
  U32 GetRV(Word *W) {
    U32 len = W->Length(), res = W->Start+len;
    if (len>=3 && ((IsVowel(W->Letters[W->Start]) && IsVowel(W->Letters[W->Start+1])) || W->StartsWith("par") || W->StartsWith("col") || W->StartsWith("tap") ))
      return W->Start+3;
    else {
      for (int i=W->Start+1;i<=W->End;i++) {
        if (IsVowel(W->Letters[i]))
          return i+1;
      }
    }
    return res;
  }
  bool Step1(Word *W, const U32 RV, const U32 R1, const U32 R2, bool *ForceStep2a) {
    int i = 0;
    for (; i<11; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R2, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        if (i==3 /*able*/)
          W->Type|=French::Adjective;
        return true;
      }
    }
    for (; i<17; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R2, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        if (W->EndsWith("ic"))
          W->ChangeSuffix("c", "qU");
        return true;
      }
    }
    for (; i<25;i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R2, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]))-1-(i<19)*2;
        if (i>22) {
          W->End+=2;
          W->Letters[W->End]='t';
        }
        return true;
      }
    }
    for (; i<27; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R1, SuffixesStep1[i]) && IsConsonant((*W)((U8)strlen(SuffixesStep1[i])))) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        return true;
      }
    }
    for (; i<29; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, RV, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        if (W->EndsWith("iv") && SuffixInRn(W, R2, "iv")) {
          W->End-=2;
          if (W->EndsWith("at") && SuffixInRn(W, R2, "at"))
            W->End-=2;
        }
        else if (W->EndsWith("eus")) {
          if (SuffixInRn(W, R2, "eus"))
            W->End-=3;
          else if (SuffixInRn(W, R1, "eus"))
            W->Letters[W->End]='x';
        }
        else if ((W->EndsWith("abl") && SuffixInRn(W, R2, "abl")) || (W->EndsWith("iqU") && SuffixInRn(W, R2, "iqU")))
          W->End-=3;
        else if ((W->EndsWith("i\xE8r") && SuffixInRn(W, RV, "i\xE8r")) || (W->EndsWith("I\xE8r") && SuffixInRn(W, RV, "I\xE8r"))) {
          W->End-=2;
          W->Letters[W->End]='i';
        }
        return true;
      }
    }
    for (; i<31; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R2, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        if (W->EndsWith("abil")) {
          if (SuffixInRn(W, R2, "abil"))
            W->End-=4;
          else
            W->End--, W->Letters[W->End]='l';
        }
        else if (W->EndsWith("ic")) {
          if (SuffixInRn(W, R2, "ic"))
            W->End-=2;
          else
            W->ChangeSuffix("c", "qU");
        }
        else if (W->EndsWith("iv") && SuffixInRn(W, R2, "iv"))
          W->End-=2;
        return true;
      }
    }
    for (; i<35; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R2, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        if (W->EndsWith("at") && SuffixInRn(W, R2, "at")) {
          W->End-=2;
          if (W->EndsWith("ic")) {
            if (SuffixInRn(W, R2, "ic"))
              W->End-=2;
            else
              W->ChangeSuffix("c", "qU");
          }
        }
        return true;
      }
    }
    for (; i<37; i++) {
      if (W->EndsWith(SuffixesStep1[i])) {
        if (SuffixInRn(W, R2, SuffixesStep1[i])) {
          W->End-=U8(strlen(SuffixesStep1[i]));
          return true;
        }
        else if (SuffixInRn(W, R1, SuffixesStep1[i])) {
          W->ChangeSuffix(SuffixesStep1[i], "eux");
          return true;
        }
      }
    }
    for (; i<NUM_SUFFIXES_STEP1; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, RV+1, SuffixesStep1[i]) && IsVowel((*W)((U8)strlen(SuffixesStep1[i])))) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        (*ForceStep2a) = true;
        return true;
      }
    }
    if (W->EndsWith("eaux") || (*W)=="eaux") {
      W->End--;
      W->Type|=French::Plural;
      return true;
    }
    else if (W->EndsWith("aux") && SuffixInRn(W, R1, "aux")) {
      W->End--, W->Letters[W->End] = 'l';
      W->Type|=French::Plural;
      return true;
    }
    else if (W->EndsWith("amment") && SuffixInRn(W, RV, "amment")) {
      W->ChangeSuffix("amment", "ant");
      (*ForceStep2a) = true;
      return true;
    }
    else if (W->EndsWith("emment") && SuffixInRn(W, RV, "emment")) {
      W->ChangeSuffix("emment", "ent");
      (*ForceStep2a) = true;
      return true;
    }
    return false;
  }
  bool Step2a(Word *W, const U32 RV) {
    for (int i=0; i<NUM_SUFFIXES_STEP2a; i++) {
      if (W->EndsWith(SuffixesStep2a[i]) && SuffixInRn(W, RV+1, SuffixesStep2a[i]) && IsConsonant((*W)((U8)strlen(SuffixesStep2a[i])))) {
        W->End-=U8(strlen(SuffixesStep2a[i]));
        if (i==31 /*ir*/)
          W->Type|=French::Verb;
        return true;
      }
    }
    return false;
  }
  bool Step2b(Word *W, const U32 RV, const U32 R2) {
    for (int i=0; i<NUM_SUFFIXES_STEP2b; i++) {
      if (W->EndsWith(SuffixesStep2b[i]) && SuffixInRn(W, RV, SuffixesStep2b[i])) {
        switch (SuffixesStep2b[i][0]) {
          case 'a': case '\xE2': {
            W->End-=U8(strlen(SuffixesStep2b[i]));
            if (W->EndsWith("e") && SuffixInRn(W, RV, "e"))
              W->End--;
            return true;
          }
          default: {
            if (i!=14 || SuffixInRn(W, R2, SuffixesStep2b[i])) {
              W->End-=U8(strlen(SuffixesStep2b[i]));
              return true;
            }
          }
        }        
      }
    }
    return false;
  }
  void Step3(Word *W) {
    char *final = (char *)&W->Letters[W->End];
    if ((*final)=='Y')
      (*final) = 'i';
    else if ((*final)=='\xE7')
      (*final) = 'c';
  }
  bool Step4(Word *W, const U32 RV, const U32 R2) {
    bool res = false;
    if (W->Length()>=2 && W->Letters[W->End]=='s' && !CharInArray((*W)(1), SetStep4, NUM_SET_STEP4)) {
      W->End--;
      res = true;
    }
    for (int i=0; i<NUM_SUFFIXES_STEP4; i++) {
      if (W->EndsWith(SuffixesStep4[i]) && SuffixInRn(W, RV, SuffixesStep4[i])) {
        switch (i) {
          case 2: { //ion
            char prec = (*W)(3);
            if (SuffixInRn(W, R2, SuffixesStep4[i]) && SuffixInRn(W, RV+1, SuffixesStep4[i]) && (prec=='s' || prec=='t')) {
              W->End-=3;
              return true;
            }
            break;
          }
          case 5: { //e
            W->End--;
            return true;
          }
          case 6: { //\xEB
            if (W->EndsWith("gu\xEB")) {
              W->End--;
              return true;
            }
            break;
          }
          default: {
            W->ChangeSuffix(SuffixesStep4[i], "i");
            return true;
          }
        }
      }
    }
    return res;
  }
  bool Step5(Word *W) {
    for (int i=0; i<NUM_SUFFIXES_STEP5; i++) {
      if (W->EndsWith(SuffixesStep5[i])) {
        W->End--;
        return true;
      }
    }
    return false;
  }
  bool Step6(Word *W) {
    for (int i=W->End; i>=W->Start; i--) {
      if (IsVowel(W->Letters[i])) {
        if (i<W->End && (W->Letters[i]&0xFE)==0xE8) {
          W->Letters[i] = 'e';
          return true;
        }
        return false;
      }
    }
    return false;
  }
public:
  inline bool IsVowel(const char c) final {
    return CharInArray(c, Vowels, NUM_VOWELS);
  }
  inline void Hash(Word *W) final {
    W->Hash[2] = W->Hash[3] = ~0xeff1cace;
    for (int i=W->Start; i<=W->End; i++) {
      U8 l = W->Letters[i];
      W->Hash[2]=W->Hash[2]*251*32 + l;
      if (IsVowel(l))
        W->Hash[3]=W->Hash[3]*997*16 + l;
      else if (l>='b' && l<='z')
        W->Hash[3]=W->Hash[3]*271*32 + (l-97);
      else
        W->Hash[3]=W->Hash[3]*11*32 + l;
    }
  }
  bool Stem(Word *W) {
    ConvertUTF8(W);
    if (W->Length()<2) {
      Hash(W);
      return false;
    }
    for (int i=0; i<NUM_EXCEPTIONS; i++) {
      if ((*W)==Exceptions[i][0]) {
        size_t len=strlen(Exceptions[i][1]);
        memcpy(&W->Letters[W->Start], Exceptions[i][1], len);
        W->End=W->Start+U8(len-1);
        Hash(W);
        W->Type|=TypesExceptions[i];
        W->Language = Language::French;
        return true;
      }
    }
    MarkVowelsAsConsonants(W);
    U32 RV=GetRV(W), R1=GetRegion(W, 0), R2=GetRegion(W, R1);
    bool DoNextStep=false, res=Step1(W, RV, R1, R2, &DoNextStep);
    DoNextStep|=!res;
    if (DoNextStep) {
      DoNextStep = !Step2a(W, RV);
      res|=!DoNextStep;
      if (DoNextStep)
        res|=Step2b(W, RV, R2);
    }
    if (res)
      Step3(W);
    else
      res|=Step4(W, RV, R2);
    res|=Step5(W);
    res|=Step6(W);
    for (int i=W->Start; i<=W->End; i++)
      W->Letters[i] = tolower(W->Letters[i]);
    if (!res)
      res=W->MatchesAny(CommonWords, NUM_COMMON_WORDS);
    Hash(W);
    if (res)
      W->Language = Language::French;
    return res;
  }
};

/*
  German suffix stemmer, based on the Porter stemmer.

  Changelog:
  (27/02/2018) v140: Initial release by Márcio Pais
*/

class GermanStemmer : public Stemmer {
private:
  static const int NUM_VOWELS = 9;
  const char Vowels[NUM_VOWELS]={'a','e','i','o','u','y','\xE4','\xF6','\xFC'};
  static const int NUM_COMMON_WORDS = 10;
  const char *CommonWords[NUM_COMMON_WORDS]={"der","die","das","und","sie","ich","mit","sich","auf","nicht"};
  static const int NUM_ENDINGS = 10;
  const char Endings[NUM_ENDINGS]={'b','d','f','g','h','k','l','m','n','t'}; //plus 'r' for words ending in 's'
  static const int NUM_SUFFIXES_STEP1 = 6;
  const char *SuffixesStep1[NUM_SUFFIXES_STEP1]={"em","ern","er","e","en","es"};
  static const int NUM_SUFFIXES_STEP2 = 3;
  const char *SuffixesStep2[NUM_SUFFIXES_STEP2]={"en","er","est"};
  static const int NUM_SUFFIXES_STEP3 = 7;
  const char *SuffixesStep3[NUM_SUFFIXES_STEP3]={"end","ung","ik","ig","isch","lich","heit"};
  void ConvertUTF8(Word *W) {
    for (int i=W->Start; i<W->End; i++) {
      U8 c = W->Letters[i+1]+((W->Letters[i+1]<0x9F)?0x60:0x40);
      if (W->Letters[i]==0xC3 && (IsVowel(c) || c==0xDF)) {
        W->Letters[i] = c;
        if (i+1<W->End)
          memcpy(&W->Letters[i+1], &W->Letters[i+2], W->End-i-1);
        W->End--;
      }
    }
  }
  void ReplaceSharpS(Word *W) {
    for (int i=W->Start; i<=W->End; i++) {
      if (W->Letters[i]==0xDF) {
        W->Letters[i]='s';
        if (i+1<MAX_WORD_SIZE) {
          memcpy(&W->Letters[i+2], &W->Letters[i+1], MAX_WORD_SIZE-i-2);
          W->Letters[i+1]='s';
          W->End+=(W->End<MAX_WORD_SIZE-1);
        }
      }
    }
  }    
  void MarkVowelsAsConsonants(Word *W) {
    for (int i=W->Start+1; i<W->End; i++) {
      U8 c = W->Letters[i];
      if ((c=='u' || c=='y') && IsVowel(W->Letters[i-1]) && IsVowel(W->Letters[i+1]))
        W->Letters[i] = toupper(c);
    }
  }
  inline bool IsValidEnding(const char c, const bool IncludeR = false) {
    return CharInArray(c, Endings, NUM_ENDINGS) || (IncludeR && c=='r');
  }
  bool Step1(Word *W, const U32 R1) {
    int i = 0;
    for (; i<3; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R1, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        return true;
      }
    }
    for (; i<NUM_SUFFIXES_STEP1; i++) {
      if (W->EndsWith(SuffixesStep1[i]) && SuffixInRn(W, R1, SuffixesStep1[i])) {
        W->End-=U8(strlen(SuffixesStep1[i]));
        W->End-=U8(W->EndsWith("niss"));
        return true;
      }
    }
    if (W->EndsWith("s") && SuffixInRn(W, R1, "s") && IsValidEnding((*W)(1), true)) {
      W->End--;
      return true;
    }
    return false;
  }
  bool Step2(Word *W, const U32 R1) {
    for (int i=0; i<NUM_SUFFIXES_STEP2; i++) {
      if (W->EndsWith(SuffixesStep2[i]) && SuffixInRn(W, R1, SuffixesStep2[i])) {
        W->End-=U8(strlen(SuffixesStep2[i]));
        return true;
      }
    }
    if (W->EndsWith("st") && SuffixInRn(W, R1, "st") && W->Length()>5 && IsValidEnding((*W)(2))) {
      W->End-=2;
      return true;
    }
    return false;
  }
  bool Step3(Word *W, const U32 R1, const U32 R2) {
    int i = 0;
    for (; i<2; i++) {
      if (W->EndsWith(SuffixesStep3[i]) && SuffixInRn(W, R2, SuffixesStep3[i])) {
        W->End-=U8(strlen(SuffixesStep3[i]));
        if (W->EndsWith("ig") && (*W)(2)!='e' && SuffixInRn(W, R2, "ig"))
          W->End-=2;
        if (i)
          W->Type|=German::Noun;
        return true;
      }
    }
    for (; i<5; i++) {
      if (W->EndsWith(SuffixesStep3[i]) && SuffixInRn(W, R2, SuffixesStep3[i]) && (*W)((U8)strlen(SuffixesStep3[i]))!='e') {
        W->End-=U8(strlen(SuffixesStep3[i]));
        if (i>2)
          W->Type|=German::Adjective;
        return true;
      }
    }
    for (; i<NUM_SUFFIXES_STEP3; i++) {
      if (W->EndsWith(SuffixesStep3[i]) && SuffixInRn(W, R2, SuffixesStep3[i])) {
        W->End-=U8(strlen(SuffixesStep3[i]));
        if ((W->EndsWith("er") || W->EndsWith("en")) && SuffixInRn(W, R1, "e?"))
          W->End-=2;
        if (i>5)
          W->Type|=German::Noun|German::Female;
        return true;
      }
    }
    if (W->EndsWith("keit") && SuffixInRn(W, R2, "keit")) {
      W->End-=4;
      if (W->EndsWith("lich") && SuffixInRn(W, R2, "lich"))
        W->End-=4;
      else if (W->EndsWith("ig") && SuffixInRn(W, R2, "ig"))
        W->End-=2;
      W->Type|=German::Noun|German::Female;
      return true;
    }
    return false;
  }
public:
  inline bool IsVowel(const char c) final {
    return CharInArray(c, Vowels, NUM_VOWELS);
  }
  inline void Hash(Word *W) final {
    W->Hash[2] = W->Hash[3] = ~0xbea7ab1e;
    for (int i=W->Start; i<=W->End; i++) {
      U8 l = W->Letters[i];
      W->Hash[2]=W->Hash[2]*263*32 + l;
      if (IsVowel(l))
        W->Hash[3]=W->Hash[3]*997*16 + l;
      else if (l>='b' && l<='z')
        W->Hash[3]=W->Hash[3]*251*32 + (l-97);
      else
        W->Hash[3]=W->Hash[3]*11*32 + l;
    }
  }
  bool Stem(Word *W) {
    ConvertUTF8(W);
    if (W->Length()<2) {
      Hash(W);
      return false;
    }
    ReplaceSharpS(W);
    MarkVowelsAsConsonants(W);
    U32 R1=GetRegion(W, 0), R2=GetRegion(W, R1);
    R1 = min(3, R1);
    bool res = Step1(W, R1);
    res|=Step2(W, R1);
    res|=Step3(W, R1, R2);
    for (int i=W->Start; i<=W->End; i++) {
      switch (W->Letters[i]) {
        case 0xE4: { W->Letters[i] = 'a'; break; }
        case 0xF6: case 0xFC: { W->Letters[i]-=0x87; break; }
        default: W->Letters[i] = tolower(W->Letters[i]);
      }
    }
    if (!res)
      res=W->MatchesAny(CommonWords, NUM_COMMON_WORDS);
    Hash(W);
    if (res)
      W->Language = Language::German;
    return res;
  }
};

//////////////////////////// Models //////////////////////////////

// All of the models below take a Mixer as a parameter and write
// predictions to it.

//////////////////////////// TextModel ///////////////////////////

template <class T, const U32 Size> class Cache {
  static_assert(Size>1 && (Size&(Size-1))==0, "Cache size must be a power of 2 bigger than 1");
private:
  Array<T> Data;
  U32 Index;
public:
  explicit Cache() : Data(Size) { Index=0; }
  T& operator()(U32 i) {
    return Data[(Index-i)&(Size-1)];
  }
  void operator++(int) {
    Index++;
  }
  void operator--(int) {
    Index--;
  }
  T& Next() {
    Index++;
    Data[Index&(Size-1)] = T();
    return Data[Index&(Size-1)];
  }
};

/*
  Text model

  Changelog:
  (04/02/2018) v135: Initial release by Márcio Pais
  (11/02/2018) v136: Uses 16 contexts, sets 3 mixer contexts
  (15/02/2018) v138: Uses 21 contexts, sets 4 mixer contexts
  (25/02/2018) v139: Uses 26 contexts
  (27/02/2018) v140: Sets 6 mixer contexts
  (12/05/2018) v142: Sets 7 mixer contexts
*/

class TextModel {
private:
  const U32 MIN_RECOGNIZED_WORDS = 4;
  ContextMap2 Map;
  Array<Stemmer*> Stemmers;
  Array<Language*> Languages;
  Cache<Word, 8> Words[Language::Count];
  Cache<Segment, 4> Segments;
  Cache<Sentence, 4> Sentences;
  Cache<Paragraph, 2> Paragraphs;
  Array<U32> WordPos;
  U32 BytePos[256];
  Word *cWord, *pWord; // current word, previous word
  Segment *cSegment; // current segment
  Sentence *cSentence; // current sentence
  Paragraph *cParagraph; // current paragraph
  enum Parse {
    Unknown,
    ReadingWord,
    PossibleHyphenation,
    WasAbbreviation,
    AfterComma,
    AfterQuote,
    AfterAbbreviation,
    ExpectDigit
  } State, pState;
  struct {
    U32 Count[Language::Count-1]; // number of recognized words of each language in the last 64 words seen
    U64 Mask[Language::Count-1];  // binary mask with the recognition status of the last 64 words for each language
    int Id;  // current language detected
    int pId; // detected language of the previous word
  } Lang;
  struct {
    U64 numbers[2];   // last 2 numbers seen
    U32 numHashes[2]; // hashes of the last 2 numbers seen
    U8  numLength[2]; // digit length of last 2 numbers seen
    U32 numMask;      // binary mask of the results of the arithmetic comparisons between the numbers seen
    U32 numDiff;      // log2 of the consecutive differences between the last 16 numbers seen, clipped to 2 bits per difference
    U32 lastUpper;    // distance to last uppercase letter
    U32 maskUpper;    // binary mask of uppercase letters seen (in the last 32 bytes)
    U32 lastLetter;   // distance to last letter
    U32 lastDigit;    // distance to last digit
    U32 lastPunct;    // distance to last punctuation character
    U32 lastNewLine;  // distance to last new line character
    U32 prevNewLine;  // distance to penultimate new line character
    U32 wordGap;      // distance between the last words
    U32 spaces;       // binary mask of whitespace characters seen (in the last 32 bytes)
    U32 spaceCount;   // count of whitespace characters seen (in the last 32 bytes)
    U32 commas;       // number of commas seen in this line (not this segment/sentence)
    U32 quoteLength;  // length (in words) of current quote
    U32 maskPunct;    // mask of relative position of last comma related to other punctuation
    U32 nestHash;     // hash representing current nesting state
    U32 lastNest;     // distance to last nesting character
    U32 masks[4],
        wordLength[2];
    int UTF8Remaining;// remaining bytes for current UTF8-encoded Unicode code point (-1 if invalid byte found)
    U8 firstLetter;   // first letter of current word
    U8 firstChar;     // first character of current line
    U8 expectedDigit; // next expected digit of detected numerical sequence
    U8 prevPunct;     // most recent punctuation character seen
    Word TopicDescriptor; // last word before ':'
  } Info;
  U32 ParseCtx;
  void Update(Buf& buffer, ModelStats *Stats = nullptr);
  void SetContexts(Buf& buffer, ModelStats *Stats = nullptr);
public:
  TextModel(const U32 Size) : Map(Size, 26), Stemmers(Language::Count-1), Languages(Language::Count-1), WordPos(0x10000), State(Parse::Unknown), pState(State), Lang{ 0, 0, Language::Unknown, Language::Unknown }, Info{ 0 }, ParseCtx(0) {
    Stemmers[Language::English-1] = new EnglishStemmer();
    Stemmers[Language::French-1] = new FrenchStemmer();
    Stemmers[Language::German-1] = new GermanStemmer();
    Languages[Language::English-1] = new English();
    Languages[Language::French-1] = new French();
    Languages[Language::German-1] = new German();
    cWord = &Words[Lang.Id](0);
    pWord = &Words[Lang.Id](1);
    cSegment = &Segments(0);
    cSentence = &Sentences(0);
    cParagraph = &Paragraphs(0);
    memset(&BytePos[0], 0, 256*sizeof(U32));
  }
  ~TextModel() {
    for (int i=0; i<Language::Count-1; i++) {
      delete Stemmers[i];
      delete Languages[i];
    }
  }
  void Predict(Mixer& mixer, Buf& buffer, ModelStats *Stats = nullptr) {
    if (bpos==0) {
      Update(buffer, Stats);
      SetContexts(buffer, Stats);
    }
    Map.mix(mixer);
    mixer.set(hash((Lang.Id!=Language::Unknown)?1+Stemmers[Lang.Id-1]->IsVowel(buffer(1)):0, Info.masks[1]&0xFF, c0)&0x7FF, 2048);
    mixer.set(hash(ilog2(Info.wordLength[0]+1), c0,
      (Info.lastDigit<Info.wordLength[0]+Info.wordGap)|
      ((Info.lastUpper<Info.lastLetter+Info.wordLength[1])<<1)|
      ((Info.lastPunct<Info.wordLength[0]+Info.wordGap)<<2)|
      ((Info.lastUpper<Info.wordLength[0])<<3)
    )&0x7FF, 2048);
    mixer.set(hash(Info.masks[1]&0x3FF, Info.lastUpper<Info.wordLength[0], Info.lastUpper<Info.lastLetter+Info.wordLength[1])&0x7FF, 2048);
    mixer.set(hash(Info.spaces&0x1FF,
      (Info.lastUpper<Info.wordLength[0])|
      ((Info.lastUpper<Info.lastLetter+Info.wordLength[1])<<1)|
      ((Info.lastPunct<Info.lastLetter)<<2)|
      ((Info.lastPunct<Info.wordLength[0]+Info.wordGap)<<3)|
      ((Info.lastPunct<Info.lastLetter+Info.wordLength[1]+Info.wordGap)<<4)
    )&0x7FF, 2048);
    mixer.set(hash(Info.firstLetter*(Info.wordLength[0]<4), min(6, Info.wordLength[0]), c0)&0x7FF, 2048);
    mixer.set(hash((*pWord)[0], (*pWord)(0), min(4, Info.wordLength[0]), Info.lastPunct<Info.lastLetter)&0x7FF, 2048);
    mixer.set(hash(min(4, Info.wordLength[0]), c0,
      Info.lastUpper<Info.wordLength[0],
      (Info.nestHash>0)?Info.nestHash&0xFF:0x100|(Info.firstLetter*(Info.wordLength[0]>0 && Info.wordLength[0]<4))
    )&0xFFF, 4096);
  }
};

void TextModel::Update(Buf& buffer, ModelStats *Stats) {
  Info.lastUpper  = min(0xFF, Info.lastUpper+1), Info.maskUpper<<=1;
  Info.lastLetter = min(0x1F, Info.lastLetter+1);
  Info.lastDigit  = min(0xFF, Info.lastDigit+1);
  Info.lastPunct  = min(0x3F, Info.lastPunct+1);
  Info.lastNewLine++, Info.prevNewLine++, Info.lastNest++;
  Info.spaceCount-=(Info.spaces>>31), Info.spaces<<=1;
  Info.masks[0]<<=2, Info.masks[1]<<=2, Info.masks[2]<<=4, Info.masks[3]<<=3;
  pState = State;  

  U8 c = buffer(1), pC=tolower(c);
  BytePos[c] = pos;
  if (c!=pC) {
    c = pC;
    Info.lastUpper = 0, Info.maskUpper|=1;
  }
  pC = buffer(2);
  ParseCtx = hash(State=Parse::Unknown, pWord->Hash[1], c, (ilog2(Info.lastNewLine)+1)*(Info.lastNewLine*3>Info.prevNewLine), Info.masks[1]&0xFC);

  if ((c>='a' && c<='z') || c=='\'' || c=='-' || c>0x7F) {    
    if (Info.wordLength[0]==0) {
      // check for hyphenation with "+"
      if (pC==NEW_LINE && ((Info.lastLetter==3 && buffer(3)=='+') || (Info.lastLetter==4 && buffer(3)==CARRIAGE_RETURN && buffer(4)=='+'))) {
        Info.wordLength[0] = Info.wordLength[1];
        for (int i=Language::Unknown; i<Language::Count; i++)
          Words[i]--;
        cWord = pWord, pWord = &Words[Lang.pId](1);
        *cWord = Word();
        for (U32 i=0; i<Info.wordLength[0]; i++)
          (*cWord)+=buffer(Info.wordLength[0]-i+Info.lastLetter);
        Info.wordLength[1] = (*pWord).Length();
        cSegment->WordCount--;
        cSentence->WordCount--;
      }
      else {
        Info.wordGap = Info.lastLetter;
        Info.firstLetter = c;
      }
    }
    Info.lastLetter = 0;
    Info.wordLength[0]++;
    Info.masks[0]+=(Lang.Id!=Language::Unknown)?1+Stemmers[Lang.Id-1]->IsVowel(c):1, Info.masks[1]++, Info.masks[3]+=Info.masks[0]&3;
    if (c=='\'') {
      Info.masks[2]+=12;
      if (Info.wordLength[0]==1) {
        if (Info.quoteLength==0 && pC==SPACE)
          Info.quoteLength = 1;
        else if (Info.quoteLength>0 && Info.lastPunct==1) {
          Info.quoteLength = 0;
          ParseCtx = hash(State=Parse::AfterQuote, pC);
        }
      }
    }
    (*cWord)+=c;
    cWord->GetHashes();
    ParseCtx = hash(State=Parse::ReadingWord, cWord->Hash[1]);
  }
  else {
    if (cWord->Length()>0) {
      if (Lang.Id!=Language::Unknown)
        memcpy(&Words[Language::Unknown](0), cWord, sizeof(Word));

      for (int i=Language::Count-1; i>Language::Unknown; i--) {
        Lang.Count[i-1]-=(Lang.Mask[i-1]>>63), Lang.Mask[i-1]<<=1;
        if (i!=Lang.Id)
          memcpy(&Words[i](0), cWord, sizeof(Word));
        if (Stemmers[i-1]->Stem(&Words[i](0)))
          Lang.Count[i-1]++, Lang.Mask[i-1]|=1;
      }      
      Lang.Id = Language::Unknown;
      U32 best = MIN_RECOGNIZED_WORDS;
      for (int i=Language::Count-1; i>Language::Unknown; i--) {
        if (Lang.Count[i-1]>=best) {
          best = Lang.Count[i-1] + (i==Lang.pId); //bias to prefer the previously detected language
          Lang.Id = i;
        }
        Words[i]++;
      }
      Words[Language::Unknown]++;
      Lang.pId = Lang.Id;
      pWord = &Words[Lang.Id](1), cWord = &Words[Lang.Id](0);
      *cWord = Word();
      WordPos[pWord->Hash[1]&(WordPos.size()-1)] = pos;
      if (cSegment->WordCount==0)
        memcpy(&cSegment->FirstWord, pWord, sizeof(Word));
      cSegment->WordCount++;
      if (cSentence->WordCount==0)
        memcpy(&cSentence->FirstWord, pWord, sizeof(Word));
      cSentence->WordCount++;
      Info.wordLength[1] = Info.wordLength[0], Info.wordLength[0] = 0;
      Info.quoteLength+=(Info.quoteLength>0);
      if (Info.quoteLength>0x1F)
        Info.quoteLength = 0;
      cSentence->VerbIndex++, cSentence->NounIndex++, cSentence->CapitalIndex++;
      if ((pWord->Type&Language::Verb)!=0) {
        cSentence->VerbIndex = 0;
        memcpy(&cSentence->lastVerb, pWord, sizeof(Word));
      }
      if ((pWord->Type&Language::Noun)!=0) {
        cSentence->NounIndex = 0;
        memcpy(&cSentence->lastNoun, pWord, sizeof(Word));
      }
      if (cSentence->WordCount>1 && Info.lastUpper<Info.wordLength[1]) {
        cSentence->CapitalIndex = 0;
        memcpy(&cSentence->lastCapital, pWord, sizeof(Word));
      }
    }
    bool skip = false;
    switch (c) {
      case '.': {
        if (Lang.Id!=Language::Unknown && Info.lastUpper==Info.wordLength[1] && Languages[Lang.Id-1]->IsAbbreviation(pWord)) {
          ParseCtx = hash(State=Parse::WasAbbreviation, pWord->Hash[1]);
          break;
        }
      }
      case '?': case '!': {
        cSentence->Type = (c=='.')?Sentence::Types::Declarative:(c=='?')?Sentence::Types::Interrogative:Sentence::Types::Exclamative;
        cSentence->SegmentCount++;
        cParagraph->SentenceCount++;
        cParagraph->TypeCount[cSentence->Type]++;
        cParagraph->TypeMask<<=2, cParagraph->TypeMask|=cSentence->Type;
        cSentence = &Sentences.Next();
        Info.masks[3]+=3;
        skip = true;
      }
      case ',': case ';': case ':': {
        if (c==',') {
          Info.commas++;
          ParseCtx = hash(State=Parse::AfterComma, ilog2(Info.quoteLength+1), ilog2(Info.lastNewLine), Info.lastUpper<Info.lastLetter+Info.wordLength[1]);
        }
        else if (c==':')
          memcpy(&Info.TopicDescriptor, pWord, sizeof(Word));
        if (!skip) {
          cSentence->SegmentCount++;
          Info.masks[3]+=4;
        }
        Info.lastPunct = 0, Info.prevPunct = c;
        Info.masks[0]+=3, Info.masks[1]+=2, Info.masks[2]+=15;
        cSegment = &Segments.Next();
        break;
      }
      case NEW_LINE: {
        Info.prevNewLine = Info.lastNewLine, Info.lastNewLine = 0;
        Info.commas = 0;
        if (Info.prevNewLine==1 || (Info.prevNewLine==2 && pC==CARRIAGE_RETURN))
          cParagraph = &Paragraphs.Next();
        else if ((Info.lastLetter==2 && pC=='+') || (Info.lastLetter==3 && pC==CARRIAGE_RETURN && buffer(3)=='+'))
          ParseCtx = hash(Parse::ReadingWord, pWord->Hash[1]), State=Parse::PossibleHyphenation;
      }
      case TAB: case CARRIAGE_RETURN: case SPACE: {
        Info.spaceCount++, Info.spaces|=1;
        Info.masks[1]+=3, Info.masks[3]+=5;
        if (c==SPACE && pState==Parse::WasAbbreviation) {
          ParseCtx = hash(State=Parse::AfterAbbreviation, pWord->Hash[1]);
        }
        break;
      }
      case '(' : Info.masks[2]+=1; Info.masks[3]+=6; Info.nestHash+=31; Info.lastNest=0; break;
      case '[' : Info.masks[2]+=2; Info.nestHash+=11; Info.lastNest=0; break;
      case '{' : Info.masks[2]+=3; Info.nestHash+=17; Info.lastNest=0; break;
      case '<' : Info.masks[2]+=4; Info.nestHash+=23; Info.lastNest=0; break;
      case 0xAB: Info.masks[2]+=5; break;
      case ')' : Info.masks[2]+=6; Info.nestHash-=31; Info.lastNest=0; break;
      case ']' : Info.masks[2]+=7; Info.nestHash-=11; Info.lastNest=0; break;
      case '}' : Info.masks[2]+=8; Info.nestHash-=17; Info.lastNest=0; break;
      case '>' : Info.masks[2]+=9; Info.nestHash-=23; Info.lastNest=0; break;
      case 0xBB: Info.masks[2]+=10; break;
      case '"': {
        Info.masks[2]+=11;
        // start/stop counting
        if (Info.quoteLength==0)
          Info.quoteLength = 1;
        else {
          Info.quoteLength = 0;
          ParseCtx = hash(State=Parse::AfterQuote, 0x100|pC);
        }
        break;
      }
      case '/' : case '-': case '+': case '*': case '=': case '%': Info.masks[2]+=13; break;
      case '\\': case '|': case '_': case '@': case '&': case '^': Info.masks[2]+=14; break;
    }
    if (c>='0' && c<='9') {
      Info.numbers[0] = Info.numbers[0]*10 + (c&0xF), Info.numLength[0] = min(19, Info.numLength[0]+1);
      Info.numHashes[0] = hash(Info.numHashes[0], c, Info.numLength[0]);
      Info.expectedDigit = -1;
      if (Info.numLength[0]<Info.numLength[1] && (pState==Parse::ExpectDigit || ((Info.numDiff&3)==0 && Info.numLength[0]<=1))) {
        U64 ExpectedNum = Info.numbers[1]+(Info.numMask&3)-2, PlaceDivisor=1;
        for (int i=0; i<Info.numLength[1]-Info.numLength[0]; i++, PlaceDivisor*=10);
        if (ExpectedNum/PlaceDivisor==Info.numbers[0]) {
          PlaceDivisor/=10;
          Info.expectedDigit = (ExpectedNum/PlaceDivisor)%10;
          State = Parse::ExpectDigit;
        }
      }
      else {
        U8 d = buffer(Info.numLength[0]+2);
        if (Info.numLength[0]<3 && buffer(Info.numLength[0]+1)==',' && d>='0' && d<='9')
          State = Parse::ExpectDigit;
      }
      Info.lastDigit = 0;
      Info.masks[3]+=7;
    }
    else if (Info.numbers[0]>0) {
      Info.numMask<<=2, Info.numMask|=1+(Info.numbers[0]>=Info.numbers[1])+(Info.numbers[0]>Info.numbers[1]);
      Info.numDiff<<=2, Info.numDiff|=min(3,ilog2(abs((int)(Info.numbers[0]-Info.numbers[1]))));
      Info.numbers[1] = Info.numbers[0], Info.numbers[0] = 0;
      Info.numHashes[1] = Info.numHashes[0], Info.numHashes[0] = 0;
      Info.numLength[1] = Info.numLength[0], Info.numLength[0] = 0;
      cSegment->NumCount++, cSentence->NumCount++;
    }
  }
  if (Info.lastNewLine==1)
    Info.firstChar = (Lang.Id!=Language::Unknown)?c:min(c,96);
  if (Info.lastNest>512)
    Info.nestHash = 0;
  int leadingBitsSet = 0;
  while (((c>>(7-leadingBitsSet))&1)!=0) leadingBitsSet++;

  if (Info.UTF8Remaining>0 && leadingBitsSet==1)
    Info.UTF8Remaining--;
  else
    Info.UTF8Remaining = (leadingBitsSet!=1)?(c!=0xC0 && c!=0xC1 && c<0xF5)?(leadingBitsSet-(leadingBitsSet>0)):-1:0;
  Info.maskPunct = (BytePos[',']>BytePos['.'])|((BytePos[',']>BytePos['!'])<<1)|((BytePos[',']>BytePos['?'])<<2)|((BytePos[',']>BytePos[':'])<<3)|((BytePos[',']>BytePos[';'])<<4);
}

void TextModel::SetContexts(Buf& buffer, ModelStats *Stats) {
  U8 c = buffer(1), lc = tolower(c), m2 = Info.masks[2]&0xF, column = min(0xFF, Info.lastNewLine);;
  U16 w = ((State==Parse::ReadingWord)?cWord->Hash[1]:pWord->Hash[1])&0xFFFF;
  U32 h = ((State==Parse::ReadingWord)?cWord->Hash[1]:pWord->Hash[2])*271+c;
  int i = State<<6;

  Map.set(ParseCtx);
  Map.set(hash(i++, cWord->Hash[0], pWord->Hash[0],
    (Info.lastUpper<Info.wordLength[0])|
    ((Info.lastDigit<Info.wordLength[0]+Info.wordGap)<<1)
  )); 
  Map.set(hash(i++, cWord->Hash[1], Words[Lang.pId](2).Hash[1], min(10,ilog2((U32)Info.numbers[0])),
    (Info.lastUpper<Info.lastLetter+Info.wordLength[1])|
    ((Info.lastLetter>3)<<1)|
    ((Info.lastLetter>0 && Info.wordLength[1]<3)<<2)
  ));
  Map.set(hash(i++, cWord->Hash[1]&0xFFF, Info.masks[1]&0x3FF, Words[Lang.pId](3).Hash[2],
    (Info.lastDigit<Info.wordLength[0]+Info.wordGap)|
    ((Info.lastUpper<Info.lastLetter+Info.wordLength[1])<<1)|
    ((Info.spaces&0x7F)<<2)
  ));
  Map.set(hash(i++, cWord->Hash[1], pWord->Hash[3], Words[Lang.pId](2).Hash[3]));
  Map.set(hash(i++, h&0x7FFF, Words[Lang.pId](2).Hash[1]&0xFFF, Words[Lang.pId](3).Hash[1]&0xFFF));
  Map.set(hash(i++, cWord->Hash[1], c, (cSentence->VerbIndex<cSentence->WordCount)?cSentence->lastVerb.Hash[1]:0));
  Map.set(hash(i++, pWord->Hash[2], Info.masks[1]&0xFC, lc, Info.wordGap));
  Map.set(hash(i++, (Info.lastLetter==0)?cWord->Hash[1]:pWord->Hash[1], c, cSegment->FirstWord.Hash[2], min(3,ilog2(cSegment->WordCount+1))));
  Map.set(hash(i++, cWord->Hash[1], c, Segments(1).FirstWord.Hash[3]));
  Map.set(hash(i++, max(31,lc), Info.masks[1]&0xFFC, (Info.spaces&0xFE)|(Info.lastPunct<Info.lastLetter), (Info.maskUpper&0xFF)|(((0x100|Info.firstLetter)*(Info.wordLength[0]>1))<<8)));
  Map.set(hash(i++, column, min(7,ilog2(Info.lastUpper+1)), ilog2(Info.lastPunct+1)));
  Map.set(
    (column&0xF8)|(Info.masks[1]&3)|((Info.prevNewLine-Info.lastNewLine>63)<<2)|
    (min(3, Info.lastLetter)<<8)|
    (Info.firstChar<<10)|
    ((Info.commas>4)<<18)|
    ((m2>=1 && m2<=5)<<19)|
    ((m2>=6 && m2<=10)<<20)|
    ((m2==11 || m2==12)<<21)|
    ((Info.lastUpper<column)<<22)|
    ((Info.lastDigit<column)<<23)|
    ((column<Info.prevNewLine-Info.lastNewLine)<<24)
  );
  Map.set(hash(
    (2*column)/3,
    min(13, Info.lastPunct)+(Info.lastPunct>16)+(Info.lastPunct>32)+Info.maskPunct*16,
    ilog2(Info.lastUpper+1),
    ilog2(Info.prevNewLine-Info.lastNewLine),
    ((Info.masks[1]&3)==0)|
    ((m2<6)<<1)|
    ((m2<11)<<2)
  ));
  Map.set(hash(i++, column>>1, Info.spaces&0xF));
  Map.set(hash(
    Info.masks[3]&0x3F,
    min((max(Info.wordLength[0],3)-2)*(Info.wordLength[0]<8),3),
    Info.firstLetter*(Info.wordLength[0]<5),
    w&0x3FF,
    (c==buffer(2))|
    ((Info.masks[2]>0)<<1)|
    ((Info.lastPunct<Info.wordLength[0]+Info.wordGap)<<2)|
    ((Info.lastUpper<Info.wordLength[0])<<3)|
    ((Info.lastDigit<Info.wordLength[0]+Info.wordGap)<<4)|
    ((Info.lastPunct<2+Info.wordLength[0]+Info.wordGap+Info.wordLength[1])<<5)
  ));
  Map.set(hash(i++, w, c, Info.numHashes[1]));
  Map.set(hash(i++, w, c, llog(pos-WordPos[w])>>1));
  Map.set(hash(i++, w, c, Info.TopicDescriptor.Hash[1]&0x7FFF));
  Map.set(hash(i++, Info.numLength[0], c, Info.TopicDescriptor.Hash[1]&0x7FFF));
  Map.set(hash(i++, (Info.lastLetter>0)?c:0x100, Info.masks[1]&0xFFC, Info.nestHash&0x7FF));
  Map.set(hash(i++, w*17+c, Info.masks[3]&0x1FF,
    ((cSentence->VerbIndex==0 && cSentence->lastVerb.Length()>0)<<6)|
    ((Info.wordLength[1]>3)<<5)|
    ((cSegment->WordCount==0)<<4)|
    ((cSentence->SegmentCount==0 && cSentence->WordCount<2)<<3)|
    ((Info.lastPunct>=Info.lastLetter+Info.wordLength[1]+Info.wordGap)<<2)|
    ((Info.lastUpper<Info.lastLetter+Info.wordLength[1])<<1)|
    (Info.lastUpper<Info.wordLength[0]+Info.wordGap+Info.wordLength[1])
  ));
  Map.set(hash(i++, c, pWord->Hash[2], Info.firstLetter*(Info.wordLength[0]<6),
    ((Info.lastPunct<Info.wordLength[0]+Info.wordGap)<<1)|
    (Info.lastPunct>=Info.lastLetter+Info.wordLength[1]+Info.wordGap)
  ));
  Map.set(hash(i++, w*23+c, Words[Lang.pId](1+(Info.wordLength[0]==0)).Letters[Words[Lang.pId](1+(Info.wordLength[0]==0)).Start], Info.firstLetter*(Info.wordLength[0]<7)));
  Map.set(hash(i++, column, Info.spaces&7, Info.nestHash&0x7FF));
  Map.set(hash(i++, cWord->Hash[1], (Info.lastUpper<column)|((Info.lastUpper<Info.wordLength[0])<<1), min(5, Info.wordLength[0])));
}

int matchModel(Mixer& m) {
  const int MAXLEN=65534;
  static Array<int> t(MEM());
  static unsigned int h=0;
  static int ptr=0;
  static int len=0;
  static int result=0;

  static SmallStationaryContextMap scm1(0x20000);

  if (!bpos) {
    h=(h*997*8+buf(1)+1)&(t.size()-1);
    if (len) ++len, ++ptr;
    else {
      ptr=t[h];
      if (ptr && pos-ptr<(int)buf.size())
        while (buf(len+1)==buf[ptr-len-1] && len<MAXLEN) ++len;
    }
    t[h]=pos;
    result=len;
    scm1.set(pos);
  }

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

void picModel(Mixer& m) {
  static U32 r0, r1, r2, r3;
  static Array<U8> t(0x10200);
  const int N=3;
  static int cxt[N];
  static StateMap sm[N];

  for (int i=0; i<N; ++i)
    t[cxt[i]]=nex(t[cxt[i]],y);

  r0+=r0+y;
  r1+=r1+((buf(215)>>(7-bpos))&1);
  r2+=r2+((buf(431)>>(7-bpos))&1);
  r3+=r3+((buf(647)>>(7-bpos))&1);
  cxt[0]=(r0&0x7)|((r1>>4)&0x38)|((r2>>3)&0xc0);
  cxt[1]=0x100+((r0&1)|((r1>>4)&0x3e)|((r2>>2)&0x40)|((r3>>1)&0x80));
  cxt[2]=0x200+((r0&0x3f)^(r1&0x3ffe)^((r2<<2)&0x7f00)^((r3<<5)&0xf800));

  for (int i=0; i<N; ++i)
    m.add(stretch(sm[i].p(t[cxt[i]])));
}

U32 b2=0,b3=0,w4=0;
U32 w5=0,f4=0,tt=0;
U32 WRT_mpw[16]= { 4, 4, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 };
U32 WRT_mtt[16]= { 0, 0, 1, 2, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7 };
U32 col=0;

static U32 frstchar=0,spafdo=0,spaces=0,spacecount=0, words=0,wordcount=0,wordlen=0,wordlen1=0;
void wordModel(Mixer& m) {
    static U32 word0=0, word1=0, word2=0, word3=0, word4=0, word5=0;
    static U32 wrdhsh=0;
    static U32 xword0=0,xword1=0,xword2=0,cword0=0,ccword=0;
    static U32 number0=0, number1=0;
    static U32 text0=0,data0=0,type0=0;
    static U32 lastLetter=0, firstLetter=0, lastUpper=0, lastDigit=0, wordGap=0;
    static ContextMap cm(MEM()*31, 61);
    static int nl1=-3, nl=-2;
    static U32 mask=0, mask2=0;
    static Array<int> wpos(0x10000);
    static int w=0;
    static Array<Word> StemWords(4);
    static Word *cWord=&StemWords[0], *pWord=&StemWords[3];
    static EnglishStemmer StemmerEN;
    static int StemIndex=0;

    if (bpos==0) {
        int c=c4&255,pC=(U8)(c4>>8),f=0;
        if (spaces&0x80000000) --spacecount;
        if (words&0x80000000) --wordcount;
        spaces=spaces*2;
        words=words*2;
        lastUpper=min(lastUpper+1,255);
        lastLetter=min(lastLetter+1,255);
        mask2<<=2;

        if (c>='A' && c<='Z') c+='a'-'A', lastUpper=0;
        if ((c>='a' && c<='z') || c=='\'' || c=='-')
           (*cWord)+=c;
        else if ((*cWord).Length()>0){
           StemmerEN.Stem(cWord);
           cWord->GetHashes();
           StemIndex=(StemIndex+1)&3;
           pWord=cWord;
           cWord=&StemWords[StemIndex];
           *cWord = Word();
        }
        if ((c>='a' && c<='z') ||  ((c>=128 &&(b3!=3)) || (c>0 && c<4 ))) { //4
            if (!wordlen){
                // model hyphenation with "+"  //book1 case +\n +\r\n
                if ((lastLetter=3 && (c4&0xFFFF00)==0x2B0A00 && buf(4)!=0x2B) || (lastLetter=4 && (c4&0xFFFFFF00)==0x2B0D0A00 && buf(5)!=0x2B) ||
                    (lastLetter=3 && (c4&0xFFFF00)==0x2D0A00 && buf(4)!=0x2D) || (lastLetter=4 && (c4&0xFFFFFF00)==0x2D0D0A00 && buf(5)!=0x2D)){
                    word0=word1;
                    word1=word2;
                    word2=word3;
                    word3=word4;
                    word4=word5;
                    word5=0;
                    wordlen = wordlen1;
                    if (c<128){
                       StemIndex=(StemIndex-1)&3;
                       cWord=pWord;
                       pWord=&StemWords[(StemIndex-1)&3];
                       *cWord = Word();
                       for (U32 i=0;i<=wordlen;i++)
                           (*cWord)+=tolower(buf(wordlen-i+1+2*(i!=wordlen)));
                    }
                }else{
                      wordGap = lastLetter;
                      firstLetter = c;
                      wrdhsh = 0;
                }
            }
            lastLetter=0;
            ++words, ++wordcount;
            if (c>4)word0^=hash(word0, c,0);
            text0=text0*997*16+c;
            wordlen++;
            wordlen=min(wordlen,45);
            f=0;
            w=word0&(wpos.size()-1);
            if ((c=='a' || c=='e' || c=='i' || c=='o' || c=='u') || (c=='y' && (wordlen>0 && pC!='a' && pC!='e' && pC!='i' && pC!='o' && pC!='u'))){
                mask2++;
                wrdhsh=wrdhsh*997*8+(c/4-22);
            }else if (c>='b' && c<='z'){
                mask2+=2;
                wrdhsh=wrdhsh*271*32+(c-97);
            }else
                wrdhsh=wrdhsh*11*32+c;
        } else {
            if (word0) {
                type0 = (type0<<2)|1;
                word5=word4;
                word4=word3;
                word3=word2;
                word2=word1;
                word1=word0;
                wordlen1=wordlen;
                 wpos[w]=blpos;
                if (c==':'|| c=='=') cword0=word0;
                if (c==']'&& (frstchar!=':')) xword0=word0;
                ccword=0;
                word0=wordlen=0;
                if((c=='.'||c=='!'||c=='?' ||c=='}' ||c==')') && buf(2)!=10) f=1;

            }
            if (c==SPACE || c==10 || c==5) { ++spaces, ++spacecount; if (c==10 || c==5) nl1=nl, nl=pos-1;}
            else if (c=='.' || c=='!' || c=='?' || c==',' || c==';' || c==':') spafdo=0,ccword=c,mask2+=3;
            else { ++spafdo; spafdo=min(63,spafdo); }
        }
        if ((c4&0xFFFF)==0x3D3D && frstchar==0x3d) xword1=word1;//,xword2=word2; // == wiki
            if ((c4&0xFFFF)==0x2727) xword2=word1;//,xword2=word2; // '' wiki
        //if ((c4&0xFFFF)==0x7D7D) xword3=word1;       //}} wiki
        lastDigit=min(0xFF,lastDigit+1);
        if (c>='0' && c<='9') {
            if (buf(3)>='0' && buf(3)<='9' && (buf(2)=='.')&& number0==0) {number0=number1; number1=0;} // 0.4645
            number0^=hash(number0, c,1);
            lastDigit = 0;
        }
        else if (number0) {
            type0 = (type0<<2)|2;
            number1=number0;
            number0=0,ccword=0;
        }
        if (!((c>='a' && c<='z') ||(c>='0' && c<='9') || (c>=128 ))){
            data0^=hash(data0, c,1);
        }else if (data0) {
            type0 = (type0<<2)|3;
            data0=0;
            }
        col=min(255, pos-nl);

        int above=buf[nl1+col]; // text column context
        if (col<=2) frstchar=(col==2?min(c,96):0);
        if (frstchar=='[' && c==32)    {if(buf(3)==']' || buf(4)==']' ) frstchar=96,xword0=0;}
          cm.set(hash(513,spafdo, spaces,ccword));
        cm.set(hash(514,frstchar, c));
        cm.set(hash(515,col, frstchar, (lastUpper<col)*4+(mask2&3)));//?
        cm.set(hash(516,spaces, (words&255)));

        cm.set(spaces&0x7fff);
        cm.set(spaces&0xff);

        cm.set(hash(257,number0, word1,wordGap));
        cm.set(hash(258,number1, c,ccword));
        cm.set(hash(259,number0, number1,wordGap));
        cm.set(hash(260,word0, number1, lastDigit<wordGap+wordlen));
        cm.set(hash(274,number0, cword0));
        cm.set(hash(518,wordlen1,col)); //?//?
        cm.set(hash(519,c,spacecount/2,wordGap));
        U32 h=wordcount*64+spacecount;
        cm.set(hash(520,c,h,ccword));  //?//?
        cm.set(hash(517,frstchar,h,lastLetter)); //+
        cm.set(hash(data0,word1, number1,type0&0xFFF));
        cm.set(hash(521,h,spafdo));
        U32 d=c4&0xf0ff;
        cm.set(hash(522,d,frstchar,ccword));

        h=word0*271;
        h=h+buf(1);

         cm.set(hash(262,h, 0));
         cm.set(hash( number0*271+buf(1), 0));
         cm.set(hash(263,word0, 0));
         if (wrdhsh) cm.set(hash(wrdhsh,buf(wpos[word1&(wpos.size()-1)]))); else cm.set(0);
         cm.set(hash(264,h, word1));
         cm.set(hash(265,word0, word1));
         cm.set(hash(266,h, word1,word2,lastUpper<wordlen));//?
         cm.set(hash(267,text0&0xffffff,0));
         //cm.set(text0&0xffffff);
         cm.set(text0&0xfffff);
         cm.set(hash(269,word0, xword0));
         cm.set(hash(270,h, xword1));
         cm.set(hash(271,h, xword2));
         cm.set(hash(272,frstchar, xword2));
         cm.set(hash(273,word0, cword0));
         cm.set(hash(275,h, word2));
         cm.set(hash(276,h, word3));
         cm.set(hash(277,h, word4));
         cm.set(hash(278,h, word5));
         cm.set(hash(279,h, word1,word3));
         cm.set(hash(280,h, word2,word3));
        cm.set(buf(1)|buf(3)<<8|buf(5)<<16);
        cm.set(buf(2)|buf(4)<<8|buf(6)<<16);
        cm.set(buf(1)|buf(4)<<8|buf(7)<<16);
        if (f) {
            word5=word4;//*29;
            word4=word3;//*31;
            word3=word2;//*37;
            word2=word1;//*41;
            word1='.';
        }
        if (col<((U32)255)){
           cm.set(hash(523,col,buf(1),above));
           cm.set(hash(524,buf(1),above));
           cm.set(hash(525,col,buf(1)));
           cm.set(hash(526,col,c==32));
        }
        else {
          cm.set(0); cm.set(0); cm.set(0); cm.set(0);
        }

       if (wordlen) cm.set(hash(281, word0, llog(blpos-wpos[word1&(wpos.size()-1)])>>4));
       else cm.set(0);

        cm.set(hash(282,buf(1),llog(blpos-wpos[word1&(wpos.size()-1)])>>2));
        cm.set(hash(283,buf(1),word0,llog(blpos-wpos[word2&(wpos.size()-1)])>>2));

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

  {
   cm.set(hash(528,mask,0));
     cm.set(hash(529,mask,buf(1)));
     cm.set(hash(530,mask&0xff,col));
     cm.set(hash(531,mask,buf(2),buf(3)));
     cm.set(hash(532,mask&0x1ff,f4&0x00fff0));
     cm.set(hash(h, llog(wordGap), mask&0x1FF,
      ((wordlen1 > 3)<<6)|
      ((wordlen > 0)<<5)|
      ((spafdo == wordlen + 2)<<4)|
      ((spafdo == wordlen + wordlen1 + 3)<<3)|
      ((spafdo >= lastLetter + wordlen1 + wordGap)<<2)|
      ((lastUpper < lastLetter + wordlen1)<<1)|
      (lastUpper < wordlen + wordlen1 + wordGap)
    ),type0&0xFFF);
    }
    if (wordlen1)    cm.set(hash(col,wordlen1,above&0x5F,c4&0x5F)); else cm.set(0); //wordlist
    if (wrdhsh)  cm.set(hash(mask2&0x3F, wrdhsh&0xFFF, (0x100|firstLetter)*(wordlen<6),(wordGap>4)*2+(wordlen1>5)) ); else cm.set(0);//?
    if ( lastLetter<16) cm.set(hash((*pWord).Hash[2], h)); else cm.set(0);
    }
    cm.mix(m);
}

void nestModel(Mixer& m)
{
  static int ic=0, bc=0, pc=0, qc=0, lvc=0, ac=0, ec=0, uc=0, sense1=0, sense2=0, w=0;
  static unsigned int vc=0, wc=0;
  static ContextMap cm(MEM()/2, 10+2);

  if (bpos==0) {
    int c=c4&255, matched=1, vv;
    w*=((vc&7)>0 && (vc&7)<3);
    if (c&0x80) w = w*11*32 + c;
    const int lc = (c >= 'A' && c <= 'Z'?c+'a'-'A':c);
    if (lc == 'a' || lc == 'e' || lc == 'i' || lc == 'o' || lc == 'u'){ vv = 1; w = w*997*8 + (lc/4-22); } else
    if (lc >= 'a' && lc <= 'z'){ vv = 2; w = w*271*32 + lc-97; } else
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
      case '(': ic += 31; break;
      case ')': ic -= 31; break;
      case '[': ic += 11; break;
      case ']': ic -= 11; break;
      case '<': ic += 23; qc += 34; break;
      case '>': ic -= 23; qc /= 5; break;
      case ':': pc = 20; break;
      case '{': ic += 17; break;
      case '}': ic -= 17; break;
      case '|': pc += 223; break;
      case '"': pc += 0x40; break;
      case '\'': pc += 0x42; if (c!=(U8)(c4>>8)) sense2^=1; else ac+=(2*sense2-1); break;
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
      case '=': pc += 87; if (c!=(U8)(c4>>8)) sense1^=1; else ec+=(2*sense1-1); break;
      default: matched = 0;
    }
    if (c4==0x266C743B) uc=min(7,uc+1);
    else if (c4==0x2667743B) uc-=(uc>0);
    if (matched) bc = 0; else bc += 1;
    if (bc > 300) bc = ic = pc = qc = uc = 0;

    cm.set(hash( (vv>0 && vv<3)?0:(lc|0x100), ic&0x3FF, ec&0x7, ac&0x7, uc ));
    cm.set(hash(ic, w, ilog2(bc+1)));
    cm.set((3*vc+77*pc+373*ic+qc)&0xffff);
    cm.set(((unsigned int)31*vc+27*pc+281*qc)&0xffff);
    cm.set(((unsigned int)13*vc+271*ic+qc+bc)&0xffff);
    cm.set((17*pc+7*ic)&0xffff);
    cm.set(((unsigned int)13*vc+ic)&0xffff);
    cm.set((vc/3+pc)&0xffff);
    cm.set(((unsigned int)7*wc+qc)&0xffff);
    cm.set((vc&0xffff)|((f4&0xf)<<16));
    cm.set(((3*pc)&0xffff)|((f4&0xf)<<16));
    cm.set((ic&0xffff)|((f4&0xf)<<16));
  }
  cm.mix(m);
}

inline U8 Clip(int Px){
  return min(0xFF,max(0,Px));
}
inline U8 Clamp4( int Px, U8 n1, U8 n2, U8 n3, U8 n4){
  return min( max(n1,max(n2,max(n3,n4))), max( min(n1,min(n2,min(n3,n4))), Px ));
}

inline U8 LogMeanDiffQt(U8 a, U8 b){
  return (a!=b)?((a>b)<<3)|ilog2((a+b)/max(2,abs(a-b)*2)+1):0;
}

struct dBASE {
  U8 Version;
  U32 nRecords;
  U16 RecordLength, HeaderLength;
  int Start, End;
};

void recordModel(Mixer& m, Filetype filetype, ModelStats *Stats = NULL) {
  static int cpos1[256] , cpos2[256], cpos3[256], cpos4[256];
  static int wpos1[0x10000];
  static int rlen[3] = {2,3,4}; // run length and 2 candidates
  static int rcount[2] = {0,0}; // candidate counts
  static U8 padding = 0; // detected padding byte
  static int prevTransition = 0, nTransition = 0; // position of the last padding transition
  static int col = 0, mxCtx = 0;
  static ContextMap cm(32768, 3), cn(32768/2, 3), co(32768*2, 3), cp(MEM(), 7);
  static StationaryMap Map0(10,8), Map1(10,8);
  static bool MayBeImg24b = false;
  static dBASE dbase;

  if (!bpos) {
    int w=c4&0xffff, c=w&255, d=w>>8;
#if 1
    if (Stats && (*Stats).Record &&
          ((*Stats).Record>>16) != (unsigned int)rlen[0]) {
      rlen[0] = (*Stats).Record>>16;
      rcount[0]=rcount[1]=0;
    }
    else{
      // detect dBASE tables
      if (blpos==0 || (dbase.Version>0 && blpos>=dbase.End))
        dbase.Version = 0;
      else if (dbase.Version==0 && (filetype==preprocessor::DEFAULT || filetype==preprocessor::TEXT) && blpos>=31){
        U8 b = buf(32);
        if ( ((b&7)==3 || (b&7)==4 || (b>>4)==3 || b==0xF5) &&
             ((b=buf(30))>0 && b<13) &&
             ((b=buf(29))>0 && b<32) &&
             ((dbase.nRecords = buf(28)|(buf(27)<<8)|(buf(26)<<16)|(buf(25)<<24)) > 0 && dbase.nRecords<0xFFFFF) &&
             ((dbase.HeaderLength = buf(24)|(buf(23)<<8)) > 32 && ( ((dbase.HeaderLength-32-1)%32)==0 || (dbase.HeaderLength>255+8 && (((dbase.HeaderLength-=255+8)-32-1)%32)==0) )) &&
             ((dbase.RecordLength = buf(22)|(buf(21)<<8)) > 8) &&
             (buf(20)==0 && buf(19)==0 && buf(17)<=1 && buf(16)<=1)
        ){
          dbase.Version = (((b=buf(32))>>4)==3)?3:b&7;
          dbase.Start = blpos - 32 + dbase.HeaderLength;
          dbase.End = dbase.Start + dbase.nRecords * dbase.RecordLength;
          if (dbase.Version==3){
            rlen[0] = 32;
            rcount[0]=rcount[1]=0;
          }
        }
      }
      else if (dbase.Version>0 && blpos==dbase.Start){
        rlen[0] = dbase.RecordLength;
        rcount[0]=rcount[1]=0;
      }

      int r=pos-cpos1[c];
      if (r>1 && r==cpos1[c]-cpos2[c]
          && r==cpos2[c]-cpos3[c] && (r>32 || r==cpos3[c]-cpos4[c])
          && (r>10 || ((c==buf(r*5+1)) && c==buf(r*6+1)))) {
        if (r==rlen[1]) ++rcount[0];
        else if (r==rlen[2]) ++rcount[1];
        else if (rcount[0]>rcount[1]) rlen[2]=r, rcount[1]=1;
        else rlen[1]=r, rcount[0]=1;
      }

      // check candidate lengths
      for (int i=0; i < 2; i++) {
        if (rcount[i] > max(0,12-(int)ilog2(rlen[i+1]))){
          if (rlen[0] != rlen[i+1]){
            if (MayBeImg24b && rlen[i+1]==3){
              rcount[0]>>=1;
              rcount[1]>>=1;
              continue;
            }
            else if ( (rlen[i+1] > rlen[0]) && (rlen[i+1] % rlen[0] == 0) ){
              // maybe we found a multiple of the real record size..?
              // in that case, it is probably an immediate multiple (2x).
              // that is probably more likely the bigger the length, so
              // check for small lengths too
              if ((rlen[0] > 32) && (rlen[i+1] == rlen[0]*2)){
                rcount[0]>>=1;
                rcount[1]>>=1;
                continue;
              }
            }
            rlen[0] = rlen[i+1];
            rcount[i] = 0;
            MayBeImg24b = (rlen[0]>30 && (rlen[0]%3)==0);
            nTransition = 0;
          }
          else
            // we found the same length again, that's positive reinforcement that
            // this really is the correct record size, so give it a little boost
            rcount[i]>>=2;

          // if the other candidate record length is orders of
          // magnitude larger, it will probably never have enough time
          // to increase its counter before it's reset again. and if
          // this length is not a multiple of the other, than it might
          // really be worthwhile to investigate it, so we won't set its
          // counter to 0
          if (rlen[i+1]<<4 > rlen[1+(i^1)])
            rcount[i^1] = 0;
        }
      }
    }
#endif
    cm.set(c<<8| (min(255, pos-cpos1[c])/4) );
    cm.set(w<<9| llog(pos-wpos1[w])>>2);
    cm.set(rlen[0]|buf(rlen[0])<<10|buf(rlen[0]*2)<<18);

    cn.set(w|rlen[0]<<8);
    cn.set(d|rlen[0]<<16);
    cn.set(c|rlen[0]<<8);

    co.set(buf(1)<<8|min(255, pos-cpos1[buf(1)]));
    co.set(buf(1)<<17|buf(2)<<9|llog(pos-wpos1[w])>>2);
    co.set(buf(1)<<8|buf(rlen[0]));

    col=pos%rlen[0];
    if (!col)
      nTransition = 0;
    cp.set(rlen[0]|buf(rlen[0])<<10|col<<18);
    cp.set(rlen[0]|buf(1)<<10|col<<18);
    cp.set(col|rlen[0]<<12);

    /*
    Consider record structures that include fixed-length strings.
    These usually contain the text followed by either spaces or 0's,
    depending on whether they're to be trimmed or they're null-terminated.
    That means we can guess the length of the string field by looking
    for small repetitions of one of these padding bytes followed by a
    different byte. By storing the last position where this transition
    ocurred, and what was the padding byte, we are able to model these
    runs of padding bytes.
    Special care is taken to skip record structures of less than 9 bytes,
    since those may be little-endian 64 bit integers. If they contain
    relatively low values (<2^40), we may consistently get runs of 3 or
    even more 0's at the end of each record, and so we could assume that
    to be the general case. But with integers, we can't be reasonably sure
    that a number won't have 3 or more 0's just before a final non-zero MSB.
    And with such simple structures, there's probably no need to be fancy
    anyway
    */

    if ((((c4>>8) == SPACE*0x010101) && (c != SPACE)) || (!(c4>>8) && c && ((padding != SPACE) || (pos-prevTransition > rlen[0])))){
      prevTransition = pos;
      nTransition+=(nTransition<31);
      padding = (U8)d;
    }
    if (rlen[0]>8){
      cp.set( hash( min(min(0xFF,rlen[0]),pos-prevTransition), min(0x3FF,col), (w&0xF0F0)|(w==((padding<<8)|padding)), nTransition ) );
      cp.set( hash( w, (buf(rlen[0]+1)==padding && buf(rlen[0])==padding), col/max(1,rlen[0]/32) ) );
    }
    else
      cp.set(0), cp.set(0);

    int last4 = (buf(rlen[0]*4)<<24)|(buf(rlen[0]*3)<<16)|(buf(rlen[0]*2)<<8)|buf(rlen[0]);
    cp.set( (last4&0xFF)|((last4&0xF000)>>4)|((last4&0xE00000)>>9)|((last4&0xE0000000)>>14)|((col/max(1,rlen[0]/16))<<18) );
    cp.set( (last4&0xF8F8)|(col<<16) );

    int i=0x300;
    if (MayBeImg24b)
      i = (col%3)<<8, Map0.set(Clip(((U8)(c4>>16))+c-(c4>>24))|i);
    else
      Map0.set(Clip(c*2-d)|i);
    Map1.set(Clip(c+buf(rlen[0])-buf(rlen[0]+1))|i);

    cpos4[c]=cpos3[c];
    cpos3[c]=cpos2[c];
    cpos2[c]=cpos1[c];
    cpos1[c]=pos;
    wpos1[w]=pos;

    mxCtx = (rlen[0]>128)?(min(0x7F,col/max(1,rlen[0]/128))):col;
  }
  cm.mix(m);
  cn.mix(m);
  co.mix(m);
  cp.mix(m);
  Map0.mix(m);
  Map1.mix(m);

  m.set( (rlen[0]>2)*( (bpos<<7)|mxCtx ), 1024 );
  m.set( ((buf(rlen[0])^((U8)(c0<<(8-bpos))))>>4)|(min(0x1F,col/max(1,rlen[0]/32))<<4), 512 );
  if (Stats)
    (*Stats).Record = (min(0xFFFF,rlen[0])<<16)|min(0xFFFF,col);
}

void recordModel1(Mixer& m) {
  static int cpos1[256];
  static int wpos1[0x10000];
  static ContextMap cm(32768, 2), cn(32768/2, 4+1), co(32768*4, 4),cp(32768*2, 3), cq(32768*2, 3);

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

    co.set(c );
    co.set(w<<8 );
    co.set(w5&0x3ffff);
    co.set(e<<3);

    cp.set(d );
    cp.set(c<<8 );
    cp.set(w<<16);

    cq.set(w<<3 );
    cq.set(c<<19);
    cq.set(e);

    cpos1[c]=pos;
    wpos1[w]=pos;
  }
  cm.mix(m);
  cn.mix(m);
  co.mix(m);
  cq.mix(m);
  cp.mix(m);
}

void sparseModel(Mixer& m, int seenbefore, int howmany) {
  static ContextMap cm(MEM()*2, 40+2);
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
   static ContextMap cm(MEM()*4, 31);
    static SmallStationaryContextMap scm1(0x10000), scm2(0x20000), scm3(0x2000),
     scm4(0x8000), scm5(0x2000),scm6(0x2000), scma(0x10000);
  if (bpos==0) {
    scm5.set(seenbefore);
    scm6.set(howmany);
  U32 h=x4<<6;
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

void distanceModel(Mixer& m) {
  static ContextMap cr(MEM(), 3);
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

inline U32 i4(int i) {
  return buf(i)+256*buf(i-1)+65536*buf(i-2)+16777216*buf(i-3);
}

inline int i2(int i) {
  return buf(i)+256*buf(i-1);
}

inline U32 m4(int i) {
  return buf(i-3)+256*buf(i-2)+65536*buf(i-1)+16777216*buf(i);
}

inline int m2(int i) {
  return buf(i)*256+buf(i-1);
}

inline int sqrbuf(int i) {
  return buf(i)*buf(i);
}

void im1bitModel(Mixer& m, int w) {
  static U32 r0, r1, r2, r3;  // last 4 rows, bit 8 is over current pixel
  static Array<U8> t(0x23000);  // model: cxt -> state
  const int N=11;  // number of contexts
  static int cxt[N];  // contexts
  static StateMap sm[N];

  // update the model
  int i;
  for (i=0; i<N; ++i)
    t[cxt[i]]=nex(t[cxt[i]],y);

  // update the contexts (pixels surrounding the predicted one)
  r0+=r0+y;
  r1+=r1+((buf(w-1)>>(7-bpos))&1);
  r2+=r2+((buf(w+w-1)>>(7-bpos))&1);
  r3+=r3+((buf(w+w+w-1)>>(7-bpos))&1);
  cxt[0]=(r0&0x7)|(r1>>4&0x38)|(r2>>3&0xc0);
  cxt[1]=0x100+((r0&1)|(r1>>4&0x3e)|(r2>>2&0x40)|(r3>>1&0x80));
  cxt[2]=0x200+((r0&1)|(r1>>4&0x1d)|(r2>>1&0x60)|(r3&0xC0));
  cxt[3]=0x300+(y|((r0<<1)&4)|((r1>>1)&0xF0)|((r2>>3)&0xA));
  cxt[4]=0x400+((r0>>4&0x2AC)|(r1&0xA4)|(r2&0x349)|(!(r3&0x14D)));
  cxt[5]=0x800+(y|((r1>>4)&0xE)|((r2>>1)&0x70)|((r3<<2)&0x380));
  cxt[6]=0xC00+(((r1&0x30)^(r3&0x0c0c))|(r0&3));
  cxt[7]=0x1000+((!(r0&0x444))|(r1&0xC0C)|(r2&0xAE3)|(r3&0x51C));
  cxt[8]=0x2000+((r0&7)|((r1>>1)&0x3F8)|((r2<<5)&0xC00));
  cxt[9]=0x3000+((r0&0x3f)^(r1&0x3ffe)^(r2<<2&0x7f00)^(r3<<5&0xf800));
  cxt[10]=0x13000+((r0&0x3e)^(r1&0x0c0c)^(r2&0xc800));

  // predict
  for (i=0; i<N; ++i) m.add(stretch(sm[i].p(t[cxt[i]])));

  m.set( (r0&7)|((r1&0x3E)>>2)|((r2&0x1C0)<<2), 2048);
  m.set( y|((r1&0x1C0)>>5)|((r2&0x1C0)>>2)|((r3&0x1C0)<<1), 1024);
  m.set( ((r1>>5)&0xFE)|y, 256);
  m.set( (r0&0x3)|((r1&0xF80)>>5), 128 );
}

//////////////////////////// im4bitModel /////////////////////////////////

// Model for 4-bit image data
void im4bitModel(Mixer& m, int w) {
  static HashTable<16> t(MEM()/2);
  const int S=11; // number of contexts
  static U8* cp[S];
  static StateMap sm[S];
  static U8 WW=0, W=0, NWW=0, NW=0, N=0, NE=0, NEE=0, NNWW = 0, NNW=0, NN=0, NNE=0, NNEE=0;
  static int col=0, line=0, run=0, prevColor=0, px=0;
  int i;
  if (!cp[0]){
    for (i=0;i<S;i++)
      cp[i]=t[263*i]+1;
  }
  for (i=0;i<S;i++)
    *cp[i]=nex(*cp[i],y);

  if (!bpos || bpos==4){
      WW=W, NWW=NW, NW=N, N=NE, NE=NEE, NNWW=NWW, NNW=NN, NN=NNE, NNE=NNEE;
      if (!bpos)
        W=c4&0xF, NEE=buf(w-1)>>4, NNEE=buf(w*2-1)>>4;
      else
        W=c0&0xF, NEE=buf(w-1)&0xF, NNEE=buf(w*2-1)&0xF;
      run=(W!=WW || !col)?(prevColor=WW,0):min(0xFFF,run+1);
      px=1, i=0;

      cp[i++]=t[hash(W,NW,N)]+1;
      cp[i++]=t[hash(N, min(0xFFF, col/8))]+1;
      cp[i++]=t[hash(W,NW,N,NN,NE)]+1;
      cp[i++]=t[hash(W, N, NE+NNE*16, NEE+NNEE*16)]+1;
      cp[i++]=t[hash(W, N, NW+NNW*16, NWW+NNWW*16)]+1;
      cp[i++]=t[hash(W, ilog2(run+1), prevColor, col/max(1,w/2) )]+1;
      cp[i++]=t[hash(NE, min(0x3FF, (col+line)/max(1,w*8)))]+1;
      cp[i++]=t[hash(NW, (col-line)/max(1,w*8))]+1;
      cp[i++]=t[hash(WW*16+W,NN*16+N,NNWW*16+NW)]+1;
      cp[i++]=t[N+NN*16]+1;
      cp[i++]=t[-1]+1;

      col*=(++col)<w*2;
      line+=(!col);
  }
  else{
    px+=px+y;
    int j=(y+1)<<(bpos&3);
    for (i=0;i<S;i++)
      cp[i]+=j;
  }

  // predict
  for (i=0; i<S; i++)
    m.add(stretch(sm[i].p(*cp[i])));

  m.set(W*16+px, 256);
  m.set(min(31,col/max(1,w/16))+N*32, 512);
  m.set((bpos&3)+4*W+64*min(7,ilog2(run+1)), 512);
  m.set(W+NE*16+(bpos&3)*256, 1024);
  m.set(px, 16);
  m.set(0,1);
}

void im8bitModel(Mixer& m, int w, int gray = 0) {
  const int nMaps = 57;
  static ContextMap cm(MEM()*4, 48);
  static StationaryMap Map[nMaps] = {
    {12,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8},
     {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8},
     {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {0,8}
  };
  static U8 WWW, WW, W, NWW, NW, N, NE, NEE, NNWW, NNW, NN, NNE, NNEE, NNN; //pixel neighborhood
  static int ctx, lastPos=0, col=0, x=0;
  static int columns[2] = {1,1}, column[2];
  // Select nearby pixels as context
  if (!bpos) {
    if (pos!=lastPos+1){
      x = 0;
      columns[0] = max(1,w/max(1,ilog2(w)*2));
      columns[1] = max(1,columns[0]/max(1,ilog2(columns[0])));
    }
    else
      x*=(++x)<w;
    lastPos = pos;
    column[0]=x/columns[0];
    column[1]=x/columns[1];
    int i=0;
    WWW=buf(3), WW=buf(2), W=buf(1), NWW=buf(w+2), NW=buf(w+1), N=buf(w), NE=buf(w-1), NEE=buf(w-2), NNWW=buf(w*2+2), NNW=buf(w*2+1), NN=buf(w*2), NNE=buf(w*2-1), NNEE=buf(w*2-2), NNN=buf(w*3);

    if (!gray){
      cm.set(hash(++i, W));
      cm.set(hash(++i, W, column[0]));
      cm.set(hash(++i, N));
      cm.set(hash(++i, N, column[0]));
      cm.set(hash(++i, NW));
      cm.set(hash(++i, NW, column[0]));
      cm.set(hash(++i, NE));
      cm.set(hash(++i, NE, column[0]));
      cm.set(hash(++i, NWW));
      cm.set(hash(++i, NEE));
      cm.set(hash(++i, WW));
      cm.set(hash(++i, NN));
      cm.set(hash(++i, W, N));
      cm.set(hash(++i, W, NW));
      cm.set(hash(++i, W, NE));
      cm.set(hash(++i, W, NEE));
      cm.set(hash(++i, W, NWW));
      cm.set(hash(++i, N, NW));
      cm.set(hash(++i, N, NE));
      cm.set(hash(++i, NW, NE));
      cm.set(hash(++i, W, WW));
      cm.set(hash(++i, N, NN));
      cm.set(hash(++i, NW, NNWW));
      cm.set(hash(++i, NE, NNEE));
      cm.set(hash(++i, NW, NWW));
      cm.set(hash(++i, NW, NNW));
      cm.set(hash(++i, NE, NEE));
      cm.set(hash(++i, NE, NNE));
      cm.set(hash(++i, N, NNW));
      cm.set(hash(++i, N, NNE));
      cm.set(hash(++i, N, NNN));
      cm.set(hash(++i, W, WWW));
      cm.set(hash(++i, WW, NEE));
      cm.set(hash(++i, WW, NN));
      cm.set(hash(++i, W, buf(w-3)));
      cm.set(hash(++i, W, buf(w-4)));
      cm.set(hash(++i, W, hash(N,NW)&0x7FF));
      cm.set(hash(++i, N, hash(NN,NNN)&0x7FF));
      cm.set(hash(++i, W, hash(NE,NEE)&0x7FF));
      cm.set(hash(++i, W, hash(NW,N,NE)&0x7FF));
      cm.set(hash(++i, N, hash(NE,NN,NNE)&0x7FF));
      cm.set(hash(++i, N, hash(NW,NNW,NN)&0x7FF));
      cm.set(hash(++i, W, hash(WW,NWW,NW)&0x7FF));
      cm.set(hash(++i, W, hash(NW,N)&0xFF, hash(WW,NWW)&0xFF));
      cm.set(hash(++i, column[0]));
      cm.set(hash(++i, N, column[1] ));
      cm.set(hash(++i, W, column[1] ));
      cm.set(++i);

      ctx = min(0x1F,(x-1)/min(0x20,columns[0]));
    }
    else{
      cm.set(hash(++i, N));
      cm.set(hash(++i, W));
      cm.set(hash(++i, NW));
      cm.set(hash(++i, NE));
      cm.set(hash(++i, N, NN));
      cm.set(hash(++i, W, WW));
      cm.set(hash(++i, NE, NNEE));
      cm.set(hash(++i, NW, NNWW));
      cm.set(hash(++i, W, NEE));
      cm.set(hash(++i, Clamp4(W+N-NW,W,NW,N,NE)/2, LogMeanDiffQt(Clip(N+NE-NNE), Clip(N+NW-NNW))));
      cm.set(hash(++i, W/4, NE/4, column[0]));
      cm.set(hash(++i, Clip(W*2-WW)/4, Clip(N*2-NN)/4));
      cm.set(hash(++i, Clamp4(N+NE-NNE,W,N,NE,NEE)/4, column[0]));
      cm.set(hash(++i, Clamp4(N+NW-NNW,W,NW,N,NE)/4, column[0]));
      cm.set(hash(++i, (W+NEE)/4, column[0]));
      cm.set(hash(++i, Clip(W+N-NW), column[0]));
      cm.set(hash(++i, Clamp4(N*3-NN*3+NNN,W,N,NN,NE), LogMeanDiffQt(W,Clip(NW*2-NNW))));
      cm.set(hash(++i, Clamp4(W*3-WW*3+WWW,W,N,NE,NEE), LogMeanDiffQt(N,Clip(NW*2-NWW))));
      cm.set(hash(++i, (W+Clamp4(NE*3-NNE*3+buf(w*3-1),W,N,NE,NEE))/2, LogMeanDiffQt(N,(NW+NE)/2)));
      cm.set(hash(++i, (N+NNN)/8, Clip(N*3-NN*3+NNN)/4));
      cm.set(hash(++i, (W+WWW)/8, Clip(W*3-WW*3+WWW)/4));
      cm.set(hash(++i, Clip((-buf(4)+5*WWW-10*WW+10*W+Clamp4(NE*4-NNE*6+buf(w*3-1)*4-buf(w*4-1),N,NE,buf(w-2),buf(w-3)))/5)));
      cm.set(hash(++i, Clip(N*2-NN), LogMeanDiffQt(N,Clip(NN*2-NNN))));
      cm.set(hash(++i, Clip(W*2-WW), LogMeanDiffQt(NE,Clip(N*2-NW))));
      cm.set(~0xde7ec7ed);

      Map[0].set( ((U8)Clip(W+N-NW))|(LogMeanDiffQt(Clip(N+NE-NNE),Clip(N+NW-NNW))<<8) );
      Map[1].set(Clamp4(W+N-NW,W,NW,N,NE));
      Map[2].set(Clip(W+N-NW));
      Map[3].set(Clamp4(W+NE-N,W,NW,N,NE));
      Map[4].set(Clip(W+NE-N));
      Map[5].set(Clamp4(N+NW-NNW,W,NW,N,NE));
      Map[6].set(Clip(N+NW-NNW));
      Map[7].set(Clamp4(N+NE-NNE,W,N,NE,NEE));
      Map[8].set(Clip(N+NE-NNE));
      Map[9].set((W+NEE)/2);
      Map[10].set(Clip(N*3-NN*3+NNN));
      Map[11].set(Clip(W*3-WW*3+WWW));
      Map[12].set((W+Clip(NE*3-NNE*3+buf(w*3-1)))/2);
      Map[13].set((W+Clip(NEE*3-buf(w*2-3)*3+buf(w*3-4)))/2);
      Map[14].set(Clip(NN+buf(w*4)-buf(w*6)));
      Map[15].set(Clip(WW+buf(4)-buf(6)));
      Map[16].set(Clip(N+NN-NNN));
      Map[17].set(Clip(W+WW-WWW));
      Map[18].set(Clip(W+NEE-NE));
      Map[19].set(Clip((buf(w*5)-6*buf(w*4)+15*NNN-20*NN+15*N+Clamp4(W*2-NWW,W,NW,N,NN))/6));
      Map[20].set(Clip((-3*WW+8*W+Clamp4(NEE*3-NNEE*3+buf(w*3-2),NE,NEE,buf(w-3),buf(w-4)))/6));
      Map[21].set(Clip(NN+NW-buf(w*3+1)));
      Map[22].set(Clip(NN+NE-buf(w*3-1)));
      Map[23].set(Clip((W*2+NW)-(WW+2*NWW)+buf(w+3)));
      Map[24].set(Clip(((NW+NWW)/2)*3-buf(w*2+3)*3+(buf(w*3+4)+buf(w*3+5))/2));
      Map[25].set(Clip(NEE+NE-buf(w*2-3)));
      Map[26].set(Clip(NWW+WW-buf(w+4)));
      Map[27].set(Clip(((W+NW)*3-NWW*6+buf(w+3)+buf(w*2+3))/2));
      Map[28].set(Clip((NE*2+NNE)-(NNEE+buf(w*3-2)*2)+buf(w*4-3)));
      Map[29].set(buf(w*6));
      Map[30].set((buf(w-4)+buf(w-6))/2);
      Map[31].set((buf(4)+buf(6))/2);
      Map[32].set((W+N+buf(w-5)+buf(w-7))/4);
      Map[33].set(Clip(buf(w-3)+W-NEE));
      Map[34].set(Clip(4*NNN-3*buf(w*4)));
      Map[35].set(Clip(WW+NEE-N));
      Map[36].set((Clip(W*2-NW)+Clip(W*2-NWW)+N+NE)/4);
      Map[37].set(Clamp4(N*2-NN,W,N,NE,NEE));
      Map[38].set((N+NNN)/2);
      Map[39].set(Clip(NN+W-NNW));
      Map[40].set(Clip(NWW+N-NNWW));
      Map[41].set(Clip((4*WWW-15*WW+20*W+NEE)/10));
      Map[42].set(Clip((buf(w*3-3)-4*NNEE+6*NE+Clip(W*3-NW*3+NNW))/4));
      Map[43].set(Clip((N*2+NE)-(NN+2*NNE)+buf(w*3-1)));
      Map[44].set(Clip((NW*2+NNW)-(NNWW+buf(w*3+2)*2)+buf(w*4+3)));
      Map[45].set(Clip(NNWW+W-buf(w*2+3)));
      Map[46].set(Clip((-buf(w*4)+5*NNN-10*NN+10*N+Clip(W*4-NWW*6+buf(w*2+3)*4-buf(w*3+4)))/5));
      Map[47].set(Clip(NEE+Clip(buf(w-3)*2-buf(w*2-4))-buf(w-4)));
      Map[48].set(Clip(NW+W-NWW));
      Map[49].set(Clip((N*2+NW)-(NN+2*NNW)+buf(w*3+1)));
      Map[50].set(Clip(NN+Clip(NEE*2-buf(w*2-3))-NNE));
      Map[51].set(Clip((-buf(4)+5*WWW-10*WW+10*W+Clip(NE*2-NNE))/5));
      Map[52].set(Clip((-buf(5)+4*buf(4)-5*WWW+5*W+Clip(NE*2-NNE))/4));
      Map[53].set(Clip((WWW-4*WW+6*W+Clip(NE*3-NNE*3+buf(w*3-1)))/4));
      Map[54].set(Clip((-NNEE+3*NE+Clip(W*4-NW*6+NNW*4-buf(w*3+1)))/3));
      Map[55].set(((W+N)*3-NW*2)/4);

      ctx = min(0x1F,x/max(1,w/min(32,columns[0])))|( ( ((abs(W-N)*16>W+N)<<1)|(abs(N-NW)>8) )<<5 )|((W+N)&0x180);
    }
  }

  cm.mix(m);
  if (gray){
    for (int i=0;i<nMaps;i++)
      Map[i].mix(m);
  }
  col=(col+1)&7;
  m.set(ctx, 2048);
  m.set(col, 8);
  m.set((N+W)>>4, 32);
  m.set(c0, 256);
  m.set( ((abs((int)(W-N))>4)<<9)|((abs((int)(N-NE))>4)<<8)|((abs((int)(W-NW))>4)<<7)|((W>N)<<6)|((N>NE)<<5)|((W>NW)<<4)|((W>WW)<<3)|((N>NN)<<2)|((NW>NNWW)<<1)|(NE>NNEE), 1024 );
  m.set(min(63,column[0]), 64);
  m.set(min(127,column[1]), 128);
}

void im24bitModel(Mixer& m, int w, int alpha=0) {
  const int SC=0x20000;
  const int nMaps = 38;
  static SmallStationaryContextMap scm1(SC), scm2(SC),
    scm3(SC), scm4(SC), scm5(SC), scm6(SC), scm7(SC), scm8(SC), scm9(SC*2), scm10(512);
  static ContextMap cm(MEM()*4, 15+32);
  static StationaryMap Map[nMaps] = { {12,8}, {12,8}, {12,8}, {12,8}, {10,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {8,8}, {0,8} };
  static U8 WWW, WW, W, NWW, NW, N, NE, NEE, NNWW, NNW, NN, NNE, NNEE, NNN; //pixel neighborhood
  static int color = -1, stride = 3;
  static int ctx[2], padding, lastPos, x = 0;
  static int columns[2] = {1,1}, column[2];

  // Select nearby pixels as context
  if (!bpos) {
    if ((color < 0) || (pos-lastPos != 1)){
      stride = 3+alpha;
      padding = w%stride;
      x = 0;
      columns[0] = max(1,w/max(1,ilog2(w)*3));
      columns[1] = max(1,columns[0]/max(1,ilog2(columns[0])));
    }
    lastPos = pos;
    x*=(++x)<w;
    if (x+padding<w)
      color*=(++color)<stride;
    else
      color=(padding>0)*(stride+1);

    int i=color<<5;
    column[0]=x/columns[0];
    column[1]=x/columns[1];

    WWW=buf(3*stride), WW=buf(2*stride), W=buf(stride), NWW=buf(w+2*stride), NW=buf(w+stride), N=buf(w), NE=buf(w-stride), NEE=buf(w-2*stride), NNWW=buf((w+stride)*2), NNW=buf(w*2+stride), NN=buf(w*2), NNE=buf(w*2-stride), NNEE=buf((w-stride)*2), NNN=buf(w*3);
    int mean=W+NW+N+NE;
    const int var=(W*W+NW*NW+N*N+NE*NE-mean*mean/4)>>2;
    mean>>=2;
    const int logvar=ilog(var);

    ctx[0] = (min(color,stride-1)<<9)|((abs(W-N)>3)<<8)|((W>N)<<7)|((W>NW)<<6)|((abs(N-NW)>3)<<5)|((N>NW)<<4)|((abs(N-NE)>3)<<3)|((N>NE)<<2)|((W>WW)<<1)|(N>NN);
    ctx[1] = ((LogMeanDiffQt(buf(1),Clip(buf(w+1)+buf(w-stride+1)-buf(w*2-stride+1)))>>1)<<5)|((LogMeanDiffQt(Clip(N+NE-NNE),Clip(N+NW-NNW))>>1)<<2)|min(color,stride-1);

    cm.set(hash((N+1)>>1, LogMeanDiffQt(N,Clip(NN*2-NNN))));
    cm.set(hash((W+1)>>1, LogMeanDiffQt(W,Clip(WW*2-WWW))));
    cm.set(hash(Clamp4(W+N-NW,W,NW,N,NE), LogMeanDiffQt(Clip(N+NE-NNE), Clip(N+NW-NNW))));
    cm.set(hash((NNN+N+4)/8, Clip(N*3-NN*3+NNN)>>1 ));
    cm.set(hash((WWW+W+4)/8, Clip(W*3-WW*3+WWW)>>1 ));
    cm.set(hash(++i, (W+Clip(NE*3-NNE*3+buf(w*3-stride)))/4, LogMeanDiffQt(N,(NW+NE)/2)));
    cm.set(hash(++i, Clip((-buf(4*stride)+5*WWW-10*WW+10*W+Clamp4(NE*4-NNE*6+buf(w*3-stride)*4-buf(w*4-stride),N,NE,buf(w-2*stride),buf(w-3*stride)))/5)/4));
    cm.set(hash(Clip(NEE+N-NNEE), LogMeanDiffQt(W,Clip(NW+NE-NNE))));
    cm.set(hash(Clip(NN+W-NNW), LogMeanDiffQt(W,Clip(NNW+WW-NNWW))));
    cm.set(hash(++i, buf(1)));
    cm.set(hash(++i, buf(2)));
    cm.set(hash(++i, Clip(W+N-NW)/2, Clip(W+buf(1)-buf(stride+1))/2));
    cm.set(hash(Clip(N*2-NN)/2, LogMeanDiffQt(N,Clip(NN*2-NNN))));
    cm.set(hash(Clip(W*2-WW)/2, LogMeanDiffQt(W,Clip(WW*2-WWW))));
    cm.set(Clamp4(N*3-NN*3+NNN,W,NW,N,NE)/2);
    cm.set(Clamp4(W*3-WW*3+WWW,W,N,NE,NEE)/2);
    cm.set(hash(++i, LogMeanDiffQt(W,buf(stride+1)), Clamp4((buf(1)*W)/max(1,buf(stride+1)),W,N,NE,NEE)));
    cm.set(hash(++i, Clamp4(N+buf(2)-buf(w+2),W,NW,N,NE)));
    cm.set(hash(++i, Clip(W+N-NW), column[0]));
    cm.set(hash(++i, Clip(N*2-NN), LogMeanDiffQt(W,Clip(NW*2-NNW))));
    cm.set(hash(++i, Clip(W*2-WW), LogMeanDiffQt(N,Clip(NW*2-NWW))));
    cm.set(hash( (W+NEE)/2, LogMeanDiffQt(W,(WW+NE)/2) ));
    cm.set(Clamp4(Clip(W*2-WW) + Clip(N*2-NN) - Clip(NW*2-NNWW), W, NW, N, NE));
    cm.set(hash(++i, W, buf(2) ));
    cm.set(hash( N, NN&0x1F, NNN&0x1F ));
    cm.set(hash( W, WW&0x1F, WWW&0x1F ));
    cm.set(hash(++i, N, column[0] ));
    cm.set(hash(++i, Clip(W+NEE-NE), LogMeanDiffQt(W,Clip(WW+NE-N))));
    cm.set(hash( NN, buf(w*4)&0x1F, buf(w*6)&0x1F, column[1] ));
    cm.set(hash( WW, buf(stride*4)&0x1F, buf(stride*6)&0x1F, column[1] ));
    cm.set(hash( NNN, buf(w*6)&0x1F, buf(w*9)&0x1F, column[1] ));
    cm.set(hash(++i, column[1]));

    cm.set(hash(++i, W, LogMeanDiffQt(W,WW)));
    cm.set(hash(++i, W, buf(1)));
    cm.set(hash(++i, W/4, LogMeanDiffQt(W,buf(1)), LogMeanDiffQt(W,buf(2)) ));
    cm.set(hash(++i, N, LogMeanDiffQt(N,NN)));
    cm.set(hash(++i, N, buf(1)));
    cm.set(hash(++i, N/4, LogMeanDiffQt(N,buf(1)), LogMeanDiffQt(N,buf(2)) ));
    cm.set(hash(++i, (W+N)>>3, buf(1)>>4, buf(2)>>4));
    cm.set(hash(++i, buf(1)/2, buf(2)/2));
    cm.set(hash(++i, W, buf(1)-buf(stride+1)));
    cm.set(hash(++i, W+buf(1)-buf(stride+1)));
    cm.set(hash(++i, N, buf(1)-buf(w+1)));
    cm.set(hash(++i, N+buf(1)-buf(w+1)));
    cm.set(hash(++i, buf(w*3-stride), buf(w*3-2*stride)));
    cm.set(hash(++i, buf(w*3+stride), buf(w*3+2*stride)));
    cm.set(hash(++i, mean, logvar>>4));
    scm1.set(W+N-NW);
    scm2.set(W+NE-N);
    scm3.set(W*2-WW);
    scm4.set(N*2-NN);
    scm5.set(NW*2-NNWW);
    scm6.set(NE*2-NNEE);
    scm7.set(NE+buf(1)-buf(w-stride+1));
    scm8.set(N+NE-NNE);
    scm9.set(mean>>1|(logvar<<1&0x180));

    Map[0].set( ((U8)Clip(W+N-NW))|(LogMeanDiffQt(Clip(N+NE-NNE),Clip(N+NW-NNW))<<8) );
    Map[1].set( ((U8)Clip(N*2-NN))|(LogMeanDiffQt(W,Clip(NW*2-NNW))<<8) );
    Map[2].set( ((U8)Clip(W*2-WW))|(LogMeanDiffQt(N,Clip(NW*2-NWW))<<8) );
    Map[3].set( ((U8)Clip(W+N-NW))|(LogMeanDiffQt(buf(1),Clip(buf(stride+1)+buf(w+1)-buf(w+stride+1)))<<8) );
    Map[4].set( (min(color,stride-1)<<8)|Clip(N+buf(1)-buf(w+1)) );
    Map[5].set( Clip((-buf(4*stride)+5*WWW-10*WW+10*W+Clamp4(NE*4-NNE*6+buf(w*3-stride)*4-buf(w*4-stride),N,NE,buf(w-2*stride),buf(w-3*stride)))/5));
    Map[6].set((W+Clamp4(NE*3-NNE*3+buf(w*3-stride),W,N,NE,NEE))/2);
    Map[7].set( Clip((buf(w*5)-6*buf(w*4)+15*NNN-20*NN+15*N+Clamp4(W*4-NWW*6+buf(w*2+3*stride)*4-buf(w*3+4*stride),W,NW,N,NN))/6) );
    Map[8].set( Clip((-3*WW+8*W+Clamp4(NEE*3-NNEE*3+buf(w*3-stride*2),NE,NEE,buf(w-3*stride),buf(w-4*stride)))/6) );
    Map[9].set( Clip((buf(w*3-3*stride)-4*NNEE+6*NE+Clip(W*4-NW*6+NNW*4-buf(w*3+stride)))/4) );
    Map[10].set( Clip(W+N-NW+buf(1)-Clip(buf(stride+1)+buf(w+1)-buf(w+stride+1))) );
    Map[11].set( Clip(W+N-NW+buf(2)-Clip(buf(stride+2)+buf(w+2)-buf(w+stride+2))) );
    Map[12].set( Clip(N*2-NN + buf(1) - Clip(buf(w+1)*2-buf(w*2+1))) );
    Map[13].set( Clip(W*2-WW + buf(1) - Clip(buf(stride+1)*2-buf(stride*2+1))) );
    Map[14].set( Clip(N+NN-NNN));
    Map[15].set( Clip(W+WW-WWW));
    Map[16].set(Clip(W+NEE-NE));
    Map[17].set(Clip(NN+W-NNW));
    Map[18].set(Clip(NN+NW-buf(w*3+stride)));
    Map[19].set(Clip(NN+NE-buf(w*3-stride)));
    Map[20].set(Clip(NN+buf(w*4)-buf(w*6)));
    Map[21].set(Clip(WW+buf(stride*4)-buf(stride*6)));
    Map[22].set(Clip(WW+NEE-N+buf(1)-Clip(buf(stride*2+1)+buf(w-stride*2+1)-buf(w+1))));
    Map[23].set(Clip(W+N-NW));
    Map[24].set(buf(1));
    Map[25].set((W+NEE)/2);
    Map[26].set(Clip(N*3-NN*3+NNN));
    Map[27].set(Clip(N+NW-NNW));
    Map[28].set(Clip(N+NE-NNE));
    Map[29].set(Clip(((N+3*NW)/4)*3-((NNW+NNWW)/2)*3+(buf(w*3+2*stride)*3+buf(w*3+3*stride))/4));
    Map[30].set(buf(w*6));
    Map[31].set((buf(w-4*stride)+buf(w-6*stride))/2);
    Map[32].set((Clip(W*2-NW)+Clip(W*2-NWW)+N+NE)/4);
    Map[33].set((buf(stride*6)+buf(stride*4))/2);
    Map[34].set(Clip((-4*WW+15*W+10*Clip(NE*3-NNE*3+buf(w*3-stride))-Clip(buf(w-3*stride)*3-buf(w*2-3*stride)*3+buf(w*3-3*stride)))/20));
    Map[35].set(((W+N)*3-NW*2)/4);
    Map[36].set(Clip((WWW-4*WW+6*W+Clip(NE*4-NNE*6+buf(w*3-stride)*4-buf(w*4-stride)))/4));
  }

  // Predict next bit
  scm1.mix(m,8);
  scm2.mix(m,8);
  scm3.mix(m,8);
  scm4.mix(m,8);
  scm5.mix(m);
  scm6.mix(m);
  scm7.mix(m);
  scm8.mix(m);
  scm9.mix(m,8);
  scm10.mix(m,8);
  cm.mix(m);
  for (int i=0;i<nMaps;i++)
    Map[i].mix(m);
  static int col=0;
  if (++col>=stride*8) col=0;
  m.set(0, 1);
  m.set(min(63,column[0])+((ctx[0]>>3)&0xC0), 256);
  m.set(min(127,column[1])+((ctx[0]>>2)&0x180), 512);
  m.set((ctx[0]&0x7FC)|(bpos>>1), 2048);
  m.set(col, stride*8);
  m.set(x%stride, stride);
  m.set(c0, 256);
  m.set((ctx[1]<<2)|(bpos>>1), 1024);
}

#define CheckIfGrayscale(x,a) {\
  if (!w && gray && (x%(3+a))==0){\
      for (int i=0; (i<(3+a)) && gray; i++){\
        U8 B=buf(4-i);\
        if (gray>>9){\
          gray=0x100|B;\
          pltorder=1-2*(B>0);\
          if (Stats)\
            (*Stats).Record = ((*Stats).Record&0xFFFF)|((3+a)<<16);\
          continue;\
        }\
        if (!i){\
          /* load first component of this entry */ \
          gray = (gray&(((int)B-(gray&0xFF)==pltorder)<<8));\
          gray|=(gray)?B:0;\
        }\
        else if (i==3)\
          gray&=((!B || (B==0xFF))*0x1FF); /* alpha/attribute component must be zero or 0xFF */\
        else\
          gray&=((B==(gray&0xFF))<<9)-1;\
      }\
  }\
}

struct TGAImage{
  U32 Header, IdLength, Bpp, ImgType, MapSize, Width, Height;
};
struct BMPImage{
  U32 Header, Offset, Bpp, Size, Palette, HdrLess, Width, Height, BitMask;
};

int imgModel(Mixer& m, ModelStats *Stats = NULL) {
  static int w=0, bpp=0;
  static int eoi=0;
  static BMPImage BMP;
  static TGAImage TGA;
  static int alpha=0, gray=0, pltorder=0;

  if (!bpos){
    // detect .BMP/DIB images
    if (pos>=(eoi+40) && !BMP.Header && (
       (buf(54)=='B' && buf(53)=='M' && ((BMP.Offset=i4(44))&0xFFFFFBF7)==0x36 && i4(40)==0x28) ||
       (BMP.HdrLess=(i4(40)==0x28))
    )){
      BMP.Width = i4(36);
      BMP.Height = abs((int)i4(32));
      BMP.Bpp = i2(26);
      BMP.Size = i4(20);
      BMP.Palette = i4(4);

      BMP.Header = (i4(24)==0) && (i2(28)==1) && (BMP.Bpp==1 || BMP.Bpp==4 || BMP.Bpp==8 || BMP.Bpp==24 || BMP.Bpp==32) && BMP.Width<30000 && BMP.Height<10000 && (!BMP.Palette || ((U32)(1<<BMP.Bpp))>=BMP.Palette);
      if (BMP.Header){
        BMP.Offset = (BMP.HdrLess)? ((BMP.Bpp<24)?((BMP.Palette)?BMP.Palette*4:4<<BMP.Bpp):0):BMP.Offset-54;
        gray = (BMP.Bpp==8)?0x300:0;
        // possible icon/cursor?
        // these have a 1bpp AND mask
        if (BMP.HdrLess && (BMP.Width*2==BMP.Height) && BMP.Bpp>1 &&
           (
            (BMP.Size>0 && BMP.Size==( (BMP.Width*BMP.Height*(BMP.Bpp+1))>>4 )) ||
            ((!BMP.Size || BMP.Size<((BMP.Width*BMP.Height*BMP.Bpp)>>3)) && (
              (BMP.Width==8)   || (BMP.Width==10) || (BMP.Width==14) || (BMP.Width==16) || (BMP.Width==20) ||
              (BMP.Width==22)  || (BMP.Width==24) || (BMP.Width==32) || (BMP.Width==40) || (BMP.Width==48) ||
              (BMP.Width==60)  || (BMP.Width==64) || (BMP.Width==72) || (BMP.Width==80) || (BMP.Width==96) ||
              (BMP.Width==128) || (BMP.Width==256)
            ))
          )
        )
          BMP.Height = BMP.BitMask = BMP.Width;
      }
    }
    else{
      BMP.Offset-=(BMP.Offset>0);
      CheckIfGrayscale(BMP.Offset,1);
    }

    if (!BMP.Offset && (BMP.Header>0 || BMP.BitMask>0) && pos>=eoi){
      if (!BMP.Header && BMP.BitMask){
        BMP.Header=BMP.Bpp=1;
        BMP.Width = BMP.BitMask;
        BMP.BitMask=0;
      }
      bpp = BMP.Bpp;
      w = (bpp>4)?(BMP.Width*(bpp>>3)+3)&(-4):(bpp==1)?(((BMP.Width-1)>>5)+1)*4:((BMP.Width*4+31)>>5)*4;
      alpha = (bpp==32);
      eoi = w*BMP.Height;
      eoi = (eoi>64)?eoi+pos:(BMP.Header=w=0);
    }

    // detect .TGA images
    if (pos>=(eoi+8) && !TGA.Header){
      if ((m4(8)&0xFFFFFF)==0x010100 && (m4(4)&0xFFFFFFC7)==0x00000100 && (buf(1)==16 || buf(1)==24 || buf(1)==32)){
        TGA.Header = pos;
        TGA.IdLength = buf(8);
        TGA.MapSize = buf(1)/8;
        TGA.Bpp = 8;
        TGA.ImgType = 1; // uncompressed color-mapped image
      }
      else if ((m4(8)&0xFFFEFF)==0x000200 && !m4(4)){
        TGA.Header = pos;
        TGA.IdLength = buf(8);
        TGA.ImgType = buf(6);
        TGA.Bpp = (TGA.ImgType==2)?24:8; // Type 2: uncompressed true-color image, no palette / Type 3: uncompressed grayscale image, no palette
      }
    }
    else if (!w && TGA.Header){
      U32 p = pos - TGA.Header;
      if (p==8){
        TGA.Width = i2(4);
        TGA.Height = i2(2);
        TGA.Header*=(!i4(8) && TGA.Width && TGA.Width<0x3FFF && TGA.Height && TGA.Height<0x3FFF);
      }
      else if (p==10){
        // 1st byte = BPP
        // 2nd byte = image spec:
        //      bits 0 - 3: 0 for 24bpp, 8 for 32bpp (not always)
        //      bit 5 indicates direction: top-bottom or bottom-top, ignored
        //      all other bits must be set to 0
        U16 i = m2(2);
        if ((i&0xFFF7)==(32<<8))
          TGA.Bpp = 32;
        if ((i&0xFFD7)!=(TGA.Bpp<<8))
          memset(&TGA, 0, sizeof(TGAImage));
      }
      if (TGA.Header && p==10 + TGA.IdLength + TGA.MapSize*256){
        w = (TGA.Width*TGA.Bpp)>>3;
        gray = (TGA.ImgType==3);
        bpp = TGA.Bpp;
        alpha = (bpp==32);
        eoi = w*TGA.Height;
        eoi = (eoi>64)?eoi+pos:(TGA.Header=w=0);
      }
    }
  }
  if (pos>eoi) return w=0;
  if (w){
    switch (bpp){
      case 1 : im1bitModel(m, w); break;
      case 4 : im4bitModel(m, w); break;
      case 8 : im8bitModel(m, w, gray); break;
      default: im24bitModel(m, w, alpha);
    }
  }
  if (bpos==7 && (pos+1)==eoi){
    memset(&TGA, 0, sizeof(TGAImage));
    BMP.Header=gray=alpha=0;
  }
  return w;
}

//////////////////////////// wavModel /////////////////////////////////

// Model a 16/8-bit stereo/mono uncompressed .wav file.
// Based on 'An asymptotically Optimal Predictor for Stereo Lossless Audio Compression'
// by Florin Ghido.

static int S,D;
static int wmode;

inline int s2(int i) { return int(short(buf(i)+256*buf(i-1))); }
inline int t2(int i) { return int(short(buf(i-1)+256*buf(i))); }

inline int X1(int i) {
  switch (wmode) {
    case 0: return buf(i)-128;
    case 1: return buf(i<<1)-128;
    case 2: return s2(i<<1);
    case 3: return s2(i<<2);
    case 4: return (buf(i)^128)-128;
    case 5: return (buf(i<<1)^128)-128;
    case 6: return t2(i<<1);
    case 7: return t2(i<<2);
    default: return 0;
  }
}

inline int X2(int i) {
  switch (wmode) {
    case 0: return buf(i+S)-128;
    case 1: return buf((i<<1)-1)-128;
    case 2: return s2((i+S)<<1);
    case 3: return s2((i<<2)-2);
    case 4: return (buf(i+S)^128)-128;
    case 5: return (buf((i<<1)-1)^128)-128;
    case 6: return t2((i+S)<<1);
    case 7: return t2((i<<2)-2);
    default: return 0;
  }
}

void wavModel(Mixer& m, int info, ModelStats *Stats = NULL) {
  static int pr[3][2], n[2], counter[2];
  static double F[49][49][2],L[49][49];
  static int rpos=0, lastPos=0;
  int j,k,l,i=0;
  long double sum;
  const double a=0.996,a2=1/a;
  const int SC=0x20000;
  static SmallStationaryContextMap scm1(SC), scm2(SC), scm3(SC), scm4(SC), scm5(SC), scm6(SC), scm7(SC);
  static ContextMap cm(MEM()*2, 10+1);
  static int bits, channels, w, ch, col=0;
  static int z1, z2, z3, z4, z5, z6, z7;

  if (!bpos){
    rpos=(pos==lastPos+1)?rpos+1:0;
    lastPos = pos;
  }

  if (!bpos && !rpos) {
    bits=((info%4)/2)*8+8;
    channels=info%2+1;
    col=0;
    w=channels*(bits>>3);
    wmode=info;
    if (channels==1) S=48,D=0; else S=36,D=12;
    for (int j=0; j<channels; j++) {
      for (k=0; k<=S+D; k++) for (l=0; l<=S+D; l++) F[k][l][j]=0, L[k][l]=0;
      F[1][0][j]=1;
      n[j]=counter[j]=pr[2][j]=pr[1][j]=pr[0][j]=0;
      z1=z2=z3=z4=z5=z6=z7=0;
    }
  }
  // Select previous samples and predicted sample as context
  if (!bpos && rpos>=w) {
    ch=rpos%w;
    const int msb=ch%(bits>>3);
    const int chn=ch/(bits>>3);
    if (!msb) {
      z1=X1(1), z2=X1(2), z3=X1(3), z4=X1(4), z5=X1(5);
      k=X1(1);
      for (l=0; l<=min(S,counter[chn]-1); l++) { F[0][l][chn]*=a; F[0][l][chn]+=X1(l+1)*k; }
      for (l=1; l<=min(D,counter[chn]); l++) { F[0][l+S][chn]*=a; F[0][l+S][chn]+=X2(l+1)*k; }
      if (channels==2) {
        k=X2(2);
        for (l=1; l<=min(D,counter[chn]); l++) { F[S+1][l+S][chn]*=a; F[S+1][l+S][chn]+=X2(l+1)*k; }
        for (l=1; l<=min(S,counter[chn]-1); l++) { F[l][S+1][chn]*=a; F[l][S+1][chn]+=X1(l+1)*k; }
        z6=X2(1)+X1(1)-X2(2), z7=X2(1);
      } else z6=2*X1(1)-X1(2), z7=X1(1);
      if (++n[chn]==1) {
        if (channels==1) for (k=1; k<=S+D; k++) for (l=k; l<=S+D; l++) F[k][l][chn]=(F[k-1][l-1][chn]-X1(k)*X1(l))*a2;
        else for (k=1; k<=S+D; k++) if (k!=S+1) for (l=k; l<=S+D; l++) if (l!=S+1) F[k][l][chn]=(F[k-1][l-1][chn]-(k-1<=S?X1(k):X2(k-S))*(l-1<=S?X1(l):X2(l-S)))*a2;
        for (i=1; i<=S+D; i++) {
           sum=F[i][i][chn];
           for (k=1; k<i; k++) sum-=L[i][k]*L[i][k];
           sum=floor(sum+0.5);
           sum=1/sum;
           if (sum>0) {
             L[i][i]=sqrt(sum);
             for (j=(i+1); j<=S+D; j++) {
               sum=F[i][j][chn];
               for (k=1; k<i; k++) sum-=L[j][k]*L[i][k];
               sum=floor(sum+0.5);
               L[j][i]=sum*L[i][i];
             }
           } else break;
        }
        if (i>S+D && counter[chn]>S+1) {
          for (k=1; k<=S+D; k++) {
            F[k][0][chn]=F[0][k][chn];
            for (j=1; j<k; j++) F[k][0][chn]-=L[k][j]*F[j][0][chn];
            F[k][0][chn]*=L[k][k];
          }
          for (k=S+D; k>0; k--) {
            for (j=k+1; j<=S+D; j++) F[k][0][chn]-=L[j][k]*F[j][0][chn];
            F[k][0][chn]*=L[k][k];
          }
        }
        n[chn]=0;
      }
      sum=0;
      for (l=1; l<=S+D; l++) sum+=F[l][0][chn]*(l<=S?X1(l):X2(l-S));
      pr[2][chn]=pr[1][chn];
      pr[1][chn]=pr[0][chn];
      pr[0][chn]=int(floor(sum));
      counter[chn]++;
    }
    const int y1=pr[0][chn], y2=pr[1][chn], y3=pr[2][chn];
    int x1=buf(1), x2=buf(2), x3=buf(3);
    if (wmode==4 || wmode==5) x1^=128, x2^=128;
    if (bits==8) x1-=128, x2-=128;
    const int t=((bits==8) || ((!msb)^(wmode<6)));
    i=ch<<4;
    if ((msb)^(wmode<6)) {
      cm.set(hash(++i, y1&0xff));
      cm.set(hash(++i, y1&0xff, ((z1-y2+z2-y3)>>1)&0xff));
      cm.set(hash(++i, x1, y1&0xff));
      cm.set(hash(++i, x1, x2>>3, x3));
      if (bits==8)
        cm.set(hash(++i, y1&0xFE, ilog2(abs((int)(z1-y2)))*2+(z1>y2) ));
      else
        cm.set(hash(++i, (y1+z1-y2)&0xff));
      cm.set(hash(++i, x1));
      cm.set(hash(++i, x1, x2));
      cm.set(hash(++i, z1&0xff));
      cm.set(hash(++i, (z1*2-z2)&0xff));
      cm.set(hash(++i, z6&0xff));
      cm.set(hash( ++i, y1&0xFF, ((z1-y2+z2-y3)/(bits>>3))&0xFF ));
    } else {
      cm.set(hash(++i, (y1-x1+z1-y2)>>8));
      cm.set(hash(++i, (y1-x1)>>8));
      cm.set(hash(++i, (y1-x1+z1*2-y2*2-z2+y3)>>8));
      cm.set(hash(++i, (y1-x1)>>8, (z1-y2+z2-y3)>>9));
      cm.set(hash(++i, z1>>12));
      cm.set(hash(++i, x1));
      cm.set(hash(++i, x1>>7, x2, x3>>7));
      cm.set(hash(++i, z1>>8));
      cm.set(hash(++i, (z1*2-z2)>>8));
      cm.set(hash(++i, y1>>8));
      cm.set(hash( ++i, (y1-x1)>>6 ));
    }
    scm1.set(t*ch);
    scm2.set(t*((z1-x1+y1)>>9)&0xff);
    scm3.set(t*((z1*2-z2-x1+y1)>>8)&0xff);
    scm4.set(t*((z1*3-z2*3+z3-x1)>>7)&0xff);
    scm5.set(t*((z1+z7-x1+y1*2)>>10)&0xff);
    scm6.set(t*((z1*4-z2*6+z3*4-z4-x1)>>7)&0xff);
    scm7.set(t*((z1*5-z2*10+z3*10-z4*5+z5-x1+y1)>>9)&0xff);
  }

  // Predict next bit
  scm1.mix(m);
  scm2.mix(m);
  scm3.mix(m);
  scm4.mix(m);
  scm5.mix(m);
  scm6.mix(m);
  scm7.mix(m);
  cm.mix(m);
  if (Stats)
    (*Stats).Record = (w<<16)|((*Stats).Record&0xFFFF);
  recordModel(m, preprocessor::AUDIO, Stats);
  if (++col>=w*8) col=0;
  m.set( ch+4*ilog2(col&(bits-1)), 4*8 );
  m.set(col%bits<8, 2);
  m.set(col%bits, bits);
  m.set(col, w*8);
  m.set(c0, 256);
}

struct WAVAudio{
  U32 Header, Size, Channels, BitsPerSample, Chunk, Data;
};

int audioModel(Mixer& m, ModelStats *Stats = NULL) {
  static U32 eoi=0, length=0, info=0;
  static WAVAudio WAV;

  if (!bpos){
    if (pos>=(int)(eoi+4) && !WAV.Header && m4(4)==0x52494646){
      WAV.Header = pos;
      WAV.Chunk = 0;
      length = 0;
    }
    else if (WAV.Header){
      int p = pos - WAV.Header;
      if (p==4){
          WAV.Size = i4(4); // file size -8 bytes, sanitize
          WAV.Header*=(WAV.Size <= 0x3FFFFFFF);
      }
      else if (p==8)
        WAV.Header*=(m4(4)==0x57415645); // "WAVE"
      else if (p==(int)(16+length) && (m4(8)!=0x666d7420 || i4(4)!=16)){ // "fmt ", chunk size=16. should be first chunk in the file, but sometimes it's not
        length = ((i4(4)+1)&(-2)) + 8; // word aligned
        WAV.Header*=!(m4(8)==0x666d7420 && i4(4)!=16); // was "fmt " chunk, but not 16bytes
      }
      else if (p==(int)(20+length)){ // check for uncompressed audio ( 0100 xx xx ) and channels ( xx xx 0? 00 )
        WAV.Channels = buf(2);
        WAV.Header*=((WAV.Channels==1 || WAV.Channels==2) && (m4(4)&0xFFFFFCFF)==0x01000000);
      }
      else if (p==(int)(32+length)){
        WAV.BitsPerSample = buf(2);
        WAV.Header*=((WAV.BitsPerSample==8 || WAV.BitsPerSample==16) && (m2(2)&0xE7FF)==0);
      }
      else if (p==(int)(40+length+WAV.Chunk) && m4(8)!=0x64617461){ // skip other chunks if not "data"
        WAV.Chunk+=((i4(4)+1)&(-2)) + 8; // size of this chunk (word aligned), plus 8 bytes for next chunk id and chunk size
        WAV.Header*=(WAV.Chunk<=0xFFFFF);
      }
      else if (p==(int)(40+length+WAV.Chunk)){
        WAV.Data = (i4(4)+1)&(-2);
        if (WAV.Data && (WAV.Data%(WAV.Channels*(WAV.BitsPerSample/8)))==0){
          info = (WAV.Channels + WAV.BitsPerSample/4-3) +1;
          eoi = pos + WAV.Data;
        }
      }
    }
  }

  if (pos>(int)eoi)
    return info=0;

  if (info)
    wavModel(m, info-1, Stats);

  if (bpos==7 && (pos+1)==(int)eoi)
    memset(&WAV, 0, sizeof(WAVAudio));

  return info;
}

class IntBuf {
  Array<int> b;
public:
  IntBuf(int i=0): b(i) {}
  int& operator[](int i) {
    return b[i&(b.size()-1)];
  }
};

#define finish(success){ \
  int length = pos - images[idx].offset; \
  memset(&images[idx], 0, sizeof(JPEGImage)); \
  mcusize=0,dqt_state=-1; \
  idx-=(idx>0); \
  images[idx].app-=length; \
  if (images[idx].app < 0) \
    images[idx].app = 0; \
}

#define jassert(x) if (!(x)) { \
  if (idx>0) \
    finish(false) \
  else \
    images[idx].jpeg=0; \
  return images[idx].next_jpeg;}

struct HUF {U32 min, max; int val;};

struct JPEGImage{
  int offset, // offset of SOI marker
  jpeg, // 1 if JPEG is header detected, 2 if image data
  next_jpeg, // updated with jpeg on next byte boundary
  app, // Bytes remaining to skip in this marker
  sof, sos, data, // pointers to buf
  htsize; // number of pointers in ht
  int ht[8]; // pointers to Huffman table headers
  U8 qtab[256]; // table
  int qmap[10]; // block -> table number
};

int jpegModel(Mixer& m) {

  // State of parser
  enum {SOF0=0xc0, SOF1, SOF2, SOF3, DHT, RST0=0xd0, SOI=0xd8, EOI, SOS, DQT,
    DNL, DRI, APP0=0xe0, COM=0xfe, FF};  // Second byte of 2 byte codes
  const static int MaxEmbeddedLevel = 3;
  static JPEGImage images[MaxEmbeddedLevel];
  static int idx=-1;
  static int lastPos=0;

  // Huffman decode state
  static U32 huffcode=0;  // Current Huffman code including extra bits
  static int huffbits=0;  // Number of valid bits in huffcode
  static int huffsize=0;  // Number of bits without extra bits
  static int rs=-1;  // Decoded huffcode without extra bits.  It represents
    // 2 packed 4-bit numbers, r=run of zeros, s=number of extra bits for
    // first nonzero code.  huffcode is complete when rs >= 0.
    // rs is -1 prior to decoding incomplete huffcode.

  static int mcupos=0;  // position in MCU (0-639).  The low 6 bits mark
    // the coefficient in zigzag scan order (0=DC, 1-63=AC).  The high
    // bits mark the block within the MCU, used to select Huffman tables.

  // Decoding tables
  static Array<HUF> huf(128);  // Tc*64+Th*16+m -> min, max, val
  static int mcusize=0;  // number of coefficients in an MCU
  static int hufsel[2][10];  // DC/AC, mcupos/64 -> huf decode table
  static Array<U8> hbuf(2048);  // Tc*1024+Th*256+hufcode -> RS

  // Image state
  static Array<int> color(10);  // block -> component (0-3)
  static Array<int> pred(4);  // component -> last DC value
  static int dc=0;  // DC value of the current block
  static int width=0;  // Image width in MCU
  static int row=0, column=0;  // in MCU (column 0 to width-1)
  static Buf cbuf(0x20000); // Rotating buffer of coefficients, coded as:
    // DC: level shifted absolute value, low 4 bits discarded, i.e.
    //   [-1023...1024] -> [0...255].
    // AC: as an RS code: a run of R (0-15) zeros followed by an S (0-15)
    //   bit number, or 00 for end of block (in zigzag order).
    //   However if R=0, then the format is ssss11xx where ssss is S,
    //   xx is the first 2 extra bits, and the last 2 bits are 1 (since
    //   this never occurs in a valid RS code).
  static int cpos=0;  // position in cbuf
  static int rs1;  // last RS code
  static int rstpos=0,rstlen=0; // reset position
  static int ssum=0, ssum1=0, ssum2=0, ssum3=0;
    // sum of S in RS codes in block and sum of S in first component

  static IntBuf cbuf2(0x20000);
  static Array<int> adv_pred(4), sumu(8), sumv(8), run_pred(6);
  static int prev_coef=0, prev_coef2=0, prev_coef_rs=0;
  static Array<int> ls(10);  // block -> distance to previous block
  static Array<int> blockW(10), blockN(10), SamplingFactors(4);
  static Array<int> lcp(7), zpos(64);

    //for parsing Quantization tables
  static int dqt_state = -1, dqt_end = 0, qnum = 0;

  const static U8 zzu[64]={  // zigzag coef -> u,v
    0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,
    3,2,1,0,1,2,3,4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};
  const static U8 zzv[64]={
    0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,
    4,5,6,7,7,6,5,4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};

  // Standard Huffman tables (cf. JPEG standard section K.3)
  // IMPORTANT: these are only valid for 8-bit data precision
  const static U8 bits_dc_luminance[16] = {
    0, 1, 5, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0
  };
  const static U8 values_dc_luminance[12] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  };

  const static U8 bits_dc_chrominance[16] = {
    0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0
  };
  const static U8 values_dc_chrominance[12] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  };

  const static U8 bits_ac_luminance[16] = {
    0, 2, 1, 3, 3, 2, 4, 3, 5, 5, 4, 4, 0, 0, 1, 0x7d
  };
  const static U8 values_ac_luminance[162] = {
    0x01, 0x02, 0x03, 0x00, 0x04, 0x11, 0x05, 0x12,
    0x21, 0x31, 0x41, 0x06, 0x13, 0x51, 0x61, 0x07,
    0x22, 0x71, 0x14, 0x32, 0x81, 0x91, 0xa1, 0x08,
    0x23, 0x42, 0xb1, 0xc1, 0x15, 0x52, 0xd1, 0xf0,
    0x24, 0x33, 0x62, 0x72, 0x82, 0x09, 0x0a, 0x16,
    0x17, 0x18, 0x19, 0x1a, 0x25, 0x26, 0x27, 0x28,
    0x29, 0x2a, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39,
    0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
    0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59,
    0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69,
    0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79,
    0x7a, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89,
    0x8a, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98,
    0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
    0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6,
    0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3, 0xc4, 0xc5,
    0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2, 0xd3, 0xd4,
    0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda, 0xe1, 0xe2,
    0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea,
    0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
    0xf9, 0xfa
  };

  const static U8 bits_ac_chrominance[16] = {
    0, 2, 1, 2, 4, 4, 3, 4, 7, 5, 4, 4, 0, 1, 2, 0x77
  };
  const static U8 values_ac_chrominance[162] = {
    0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21,
    0x31, 0x06, 0x12, 0x41, 0x51, 0x07, 0x61, 0x71,
    0x13, 0x22, 0x32, 0x81, 0x08, 0x14, 0x42, 0x91,
    0xa1, 0xb1, 0xc1, 0x09, 0x23, 0x33, 0x52, 0xf0,
    0x15, 0x62, 0x72, 0xd1, 0x0a, 0x16, 0x24, 0x34,
    0xe1, 0x25, 0xf1, 0x17, 0x18, 0x19, 0x1a, 0x26,
    0x27, 0x28, 0x29, 0x2a, 0x35, 0x36, 0x37, 0x38,
    0x39, 0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48,
    0x49, 0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58,
    0x59, 0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68,
    0x69, 0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
    0x79, 0x7a, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
    0x88, 0x89, 0x8a, 0x92, 0x93, 0x94, 0x95, 0x96,
    0x97, 0x98, 0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5,
    0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4,
    0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3,
    0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2,
    0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda,
    0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9,
    0xea, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
    0xf9, 0xfa
  };

  if (idx < 0){
    memset(&images[0], 0, sizeof(images));
    idx = 0;
    lastPos = pos;
  }

  // Be sure to quit on a byte boundary
  if (!bpos) images[idx].next_jpeg=images[idx].jpeg>1;
  if (bpos && !images[idx].jpeg) return images[idx].next_jpeg;
  if (!bpos && images[idx].app>0){
    --images[idx].app;
    if (idx<MaxEmbeddedLevel && buf(4)==FF && buf(3)==SOI && buf(2)==FF && ((buf(1)&0xFE)==0xC0 || buf(1)==0xC4 || (buf(1)>=0xDB && buf(1)<=0xFE)) )
      memset(&images[++idx], 0, sizeof(JPEGImage));
  }
  if (images[idx].app>0) return images[idx].next_jpeg;
  if (!bpos) {

    // Parse.  Baseline DCT-Huffman JPEG syntax is:
    // SOI APPx... misc... SOF0 DHT... SOS data EOI
    // SOI (= FF D8) start of image.
    // APPx (= FF Ex) len ... where len is always a 2 byte big-endian length
    //   including the length itself but not the 2 byte preceding code.
    //   Application data is ignored.  There may be more than one APPx.
    // misc codes are DQT, DNL, DRI, COM (ignored).
    // SOF0 (= FF C0) len 08 height width Nf [C HV Tq]...
    //   where len, height, width (in pixels) are 2 bytes, Nf is the repeat
    //   count (1 byte) of [C HV Tq], where C is a component identifier
    //   (color, 0-3), HV is the horizontal and vertical dimensions
    //   of the MCU (high, low bits, packed), and Tq is the quantization
    //   table ID (not used).  An MCU (minimum compression unit) consists
    //   of 64*H*V DCT coefficients for each color.
    // DHT (= FF C4) len [TcTh L1...L16 V1,1..V1,L1 ... V16,1..V16,L16]...
    //   defines Huffman table Th (1-4) for Tc (0=DC (first coefficient)
    //   1=AC (next 63 coefficients)).  L1..L16 are the number of codes
    //   of length 1-16 (in ascending order) and Vx,y are the 8-bit values.
    //   A V code of RS means a run of R (0-15) zeros followed by S (0-15)
    //   additional bits to specify the next nonzero value, negative if
    //   the first additional bit is 0 (e.g. code x63 followed by the
    //   3 bits 1,0,1 specify 7 coefficients: 0, 0, 0, 0, 0, 0, 5.
    //   Code 00 means end of block (remainder of 63 AC coefficients is 0).
    // SOS (= FF DA) len Ns [Cs TdTa]... 0 3F 00
    //   Start of scan.  TdTa specifies DC/AC Huffman tables (0-3, packed
    //   into one byte) for component Cs matching C in SOF0, repeated
    //   Ns (1-4) times.
    // EOI (= FF D9) is end of image.
    // Huffman coded data is between SOI and EOI.  Codes may be embedded:
    // RST0-RST7 (= FF D0 to FF D7) mark the start of an independently
    //   compressed region.
    // DNL (= FF DC) 04 00 height
    //   might appear at the end of the scan (ignored).
    // FF 00 is interpreted as FF (to distinguish from RSTx, DNL, EOI).

    // Detect JPEG (SOI followed by a valid marker)
    if (!images[idx].jpeg && buf(4)==FF && buf(3)==SOI && buf(2)==FF && ((buf(1)&0xFE)==0xC0 || buf(1)==0xC4 || (buf(1)>=0xDB && buf(1)<=0xFE)) ){
      images[idx].jpeg=1;
      images[idx].offset = pos-4;
      images[idx].sos=images[idx].sof=images[idx].htsize=images[idx].data=0, images[idx].app=(buf(1)>>4==0xE)*2;
      mcusize=huffcode=huffbits=huffsize=mcupos=cpos=0, rs=-1;
      memset(&huf[0], 0, sizeof(huf));
      memset(&pred[0], 0, pred.size()*sizeof(int));
      rstpos=rstlen=0;
    }

    // Detect end of JPEG when data contains a marker other than RSTx
    // or byte stuff (00), or if we jumped in position since the last byte seen
    if (images[idx].jpeg && images[idx].data && ((buf(2)==FF && buf(1) && (buf(1)&0xf8)!=RST0) || (pos-lastPos>1)) ) {
      jassert((buf(1)==EOI) || (pos-lastPos>1));
      finish(true);
    }
    lastPos = pos;
    if (!images[idx].jpeg) return images[idx].next_jpeg;

    // Detect APPx, COM or other markers, so we can skip them
    if (!images[idx].data && !images[idx].app && buf(4)==FF && (((buf(3)>0xC1) && (buf(3)<=0xCF) && (buf(3)!=DHT)) || ((buf(3)>=0xDC) && (buf(3)<=0xFE)))){
      images[idx].app=buf(2)*256+buf(1)+2;
      if (idx>0)
        jassert( pos + images[idx].app < images[idx].offset + images[idx-1].app );
    }

    // Save pointers to sof, ht, sos, data,
    if (buf(5)==FF && buf(4)==SOS) {
      int len=buf(3)*256+buf(2);
      if (len==6+2*buf(1) && buf(1) && buf(1)<=4)  // buf(1) is Ns
        images[idx].sos=pos-5, images[idx].data=images[idx].sos+len+2, images[idx].jpeg=2;
    }
    if (buf(4)==FF && buf(3)==DHT && images[idx].htsize<8) images[idx].ht[images[idx].htsize++]=pos-4;
    if (buf(4)==FF && (buf(3)&0xFE)==SOF0) images[idx].sof=pos-4;

    // Parse Quantizazion tables
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
          images[idx].qtab[qnum*64+((dqt_state%65)-1)]=buf(1)-1;
        }
        dqt_state++;
      }
    }

    // Restart
    if (buf(2)==FF && (buf(1)&0xf8)==RST0) {
      huffcode=huffbits=huffsize=mcupos=0, rs=-1;
      memset(&pred[0], 0, pred.size()*sizeof(int));
      rstlen=column+row*width-rstpos;
      rstpos=column+row*width;
    }
  }

  {
    // Build Huffman tables
    // huf[Tc][Th][m] = min, max+1 codes of length m, pointer to byte values
    if (pos==images[idx].data && bpos==1) {
      int i;
      for (i=0; i<images[idx].htsize; ++i) {
        int p=images[idx].ht[i]+4;  // pointer to current table after length field
        int end=p+buf[p-2]*256+buf[p-1]-2;  // end of Huffman table
        int count=0;  // sanity check
        while (p<end && end<pos && end<p+2100 && ++count<10) {
          int tc=buf[p]>>4, th=buf[p]&15;
          if (tc>=2 || th>=4) break;
          jassert(tc>=0 && tc<2 && th>=0 && th<4);
          HUF* h=&huf[tc*64+th*16]; // [tc][th][0];
          int val=p+17;  // pointer to values
          int hval=tc*1024+th*256;  // pointer to RS values in hbuf
          int j;
          for (j=0; j<256; ++j) // copy RS codes
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

      // load default tables
      if (!images[idx].htsize){
        for (int tc = 0; tc < 2; tc++) {
          for (int th = 0; th < 2; th++) {
            HUF* h = &huf[tc*64+th*16];
            int hval = tc*1024 + th*256;
            int code = 0, c = 0, x = 0;

            for (int i = 0; i < 16; i++) {
              switch (tc*2+th) {
                case 0: x = bits_dc_luminance[i]; break;
                case 1: x = bits_dc_chrominance[i]; break;
                case 2: x = bits_ac_luminance[i]; break;
                case 3: x = bits_ac_chrominance[i];
              }

              h[i].min = code;
              h[i].max = (code+=x);
              h[i].val = hval;
              hval+=x;
              code+=code;
              c+=x;
            }

            hval = tc*1024 + th*256;
            c--;

            while (c >= 0){
              switch (tc*2+th) {
                case 0: x = values_dc_luminance[c]; break;
                case 1: x = values_dc_chrominance[c]; break;
                case 2: x = values_ac_luminance[c]; break;
                case 3: x = values_ac_chrominance[c];
              }

              hbuf[hval+c] = x;
              c--;
            }
          }
        }
        images[idx].htsize = 4;
      }

      // Build Huffman table selection table (indexed by mcupos).
      // Get image width.
      if (!images[idx].sof && images[idx].sos) return images[idx].next_jpeg;
      int ns=buf[images[idx].sos+4];
      int nf=buf[images[idx].sof+9];
      jassert(ns<=4 && nf<=4);
      mcusize=0;  // blocks per MCU
      int hmax=0;  // MCU horizontal dimension
      for (i=0; i<ns; ++i) {
        for (int j=0; j<nf; ++j) {
          if (buf[images[idx].sos+2*i+5]==buf[images[idx].sof+3*j+10]) { // Cs == C ?
            int hv=buf[images[idx].sof+3*j+11];  // packed dimensions H x V
            SamplingFactors[j] = hv;
            if (hv>>4>hmax) hmax=hv>>4;
            hv=(hv&15)*(hv>>4);  // number of blocks in component C
            jassert(hv>=1 && hv+mcusize<=10);
            while (hv) {
              jassert(mcusize<10);
              hufsel[0][mcusize]=buf[images[idx].sos+2*i+6]>>4&15;
              hufsel[1][mcusize]=buf[images[idx].sos+2*i+6]&15;
              jassert (hufsel[0][mcusize]<4 && hufsel[1][mcusize]<4);
              color[mcusize]=i;
              int tq=buf[images[idx].sof+3*j+12];  // quantization table index (0..3)
              jassert(tq>=0 && tq<4);
              images[idx].qmap[mcusize]=tq; // quantizazion table mapping
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
      width=buf[images[idx].sof+7]*256+buf[images[idx].sof+8];  // in pixels
      width=(width-1)/(hmax*8)+1;  // in MCU
      jassert(width>0);
      mcusize*=64;  // coefficients per MCU
      row=column=0;

      // we can have more blocks than components then we have subsampling
      int x=0, y=0;
      for (j = 0; j<(mcusize>>6); j++) {
        int i = color[j];
        int w = SamplingFactors[i]>>4, h = SamplingFactors[i]&0xf;
        blockW[j] = x==0?mcusize-64*(w-1):64;
        blockN[j] = y==0?mcusize*width-64*w*(h-1):w*64;
        x++;
        if (x>=w) { x=0; y++; }
        if (y>=h) { x=0; y=0; }
      }
    }
  }

  // Decode Huffman
  {
    if (mcusize && buf(1+(!bpos))!=FF) {  // skip stuffed byte
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
        const HUF *h=&huf[ac*64+sel*16]; // [ac][sel];
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
        if (huffsize+(rs&15)==huffbits) { // done decoding
          rs1=rs;
          int x=0;  // decoded extra bits
          if (mcupos&63) {  // AC
            if (rs==0) { // EOB
              mcupos=(mcupos+63)&-64;
              jassert(mcupos>=0 && mcupos<=mcusize && mcupos<=640);
              while (cpos&63) {
                cbuf2[cpos]=0;
                cbuf[cpos]=(!rs)?0:(63-(cpos&63))<<4; cpos++; rs++;
              }
            }
            else {  // rs = r zeros + s extra bits for the next nonzero value
                    // If first extra bit is 0 then value is negative.
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
          else {  // DC: rs = 0S, s<12
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

          // UPDATE_ADV_PRED !!!!
          {
            const int acomp=mcupos>>6, q=64*images[idx].qmap[acomp];
            const int zz=mcupos&63, cpos_dc=cpos-zz;
            const bool norst=rstpos!=column+row*width;
            if (zz==0) {
              for (int i=0; i<8; ++i) sumu[i]=sumv[i]=0;
              // position in the buffer of first (DC) coefficient of the block
              // of this same component that is to the west of this one, not
              // necessarily in this MCU
              int offset_DC_W = cpos_dc - blockW[acomp];
              // position in the buffer of first (DC) coefficient of the block
              // of this same component that is to the north of this one, not
              // necessarily in this MCU
              int offset_DC_N = cpos_dc - blockN[acomp];
              for (int i=0; i<64; ++i) {
                sumu[zzu[i]]+=(zzv[i]&1?-1:1)*(zzv[i]?16*(16+zzv[i]):185)*(images[idx].qtab[q+i]+1)*cbuf2[offset_DC_N+i];
                sumv[zzv[i]]+=(zzu[i]&1?-1:1)*(zzu[i]?16*(16+zzu[i]):185)*(images[idx].qtab[q+i]+1)*cbuf2[offset_DC_W+i];
              }
            }
            else {
              sumu[zzu[zz-1]]-=(zzv[zz-1]?16*(16+zzv[zz-1]):185)*(images[idx].qtab[q+zz-1]+1)*cbuf2[cpos-1];
              sumv[zzv[zz-1]]-=(zzu[zz-1]?16*(16+zzu[zz-1]):185)*(images[idx].qtab[q+zz-1]+1)*cbuf2[cpos-1];
            }

            for (int i=0; i<3; ++i)
            {
              run_pred[i]=run_pred[i+3]=0;
              for (int st=0; st<10 && zz+st<64; ++st) {
                const int zz2=zz+st;
                int p=sumu[zzu[zz2]]*i+sumv[zzv[zz2]]*(2-i);
                p/=(images[idx].qtab[q+zz2]+1)*185*(16+zzv[zz2])*(16+zzu[zz2])/128;
                if (zz2==0 && (norst || ls[acomp]==64)) p-=cbuf2[cpos_dc-ls[acomp]];
                p=(p<0?-1:+1)*ilog(abs(p)+1);
                if (st==0) {
                  adv_pred[i]=p;
                }
                else if (abs(p)>abs(adv_pred[i])+2 && abs(adv_pred[i]) < 210) {
                  if (run_pred[i]==0) run_pred[i]=st*2+(p>0);
                  if (abs(p)>abs(adv_pred[i])+21 && run_pred[i+3]==0) run_pred[i+3]=st*2+(p>0);
                }
              }
            }
            x=0;
            for (int i=0; i<8; ++i) x+=(zzu[zz]<i)*sumu[i]+(zzv[zz]<i)*sumv[i];
            x=(sumu[zzu[zz]]*(2+zzu[zz])+sumv[zzv[zz]]*(2+zzv[zz])-x*2)*4/(zzu[zz]+zzv[zz]+16);
            x/=(images[idx].qtab[q+zz]+1)*185;
            if (zz==0 && (norst || ls[acomp]==64)) x-=cbuf2[cpos_dc-ls[acomp]];
            adv_pred[3]=(x<0?-1:+1)*ilog(abs(x)+1);

            for (int i=0; i<4; ++i) {
              const int a=(i&1?zzv[zz]:zzu[zz]), b=(i&2?2:1);
              if (a<b) x=65535;
              else {
                const int zz2=zpos[zzu[zz]+8*zzv[zz]-(i&1?8:1)*b];
                x=(images[idx].qtab[q+zz2]+1)*cbuf2[cpos_dc+zz2]/(images[idx].qtab[q+zz]+1);
                x=(x<0?-1:+1)*(ilog(abs(x)+1)+(x!=0?17:0));
              }
              lcp[i]=x;
            }
            if ((zzu[zz]*zzv[zz])!=0){
              const int zz2=zpos[zzu[zz]+8*zzv[zz]-9];
              x=(images[idx].qtab[q+zz2]+1)*cbuf2[cpos_dc+zz2]/(images[idx].qtab[q+zz]+1);
              lcp[4]=(x<0?-1:+1)*(ilog(abs(x)+1)+(x!=0?17:0));

              x=(images[idx].qtab[q+zpos[8*zzv[zz]]]+1)*cbuf2[cpos_dc+zpos[8*zzv[zz]]]/(images[idx].qtab[q+zz]+1);
              lcp[5]=(x<0?-1:+1)*(ilog(abs(x)+1)+(x!=0?17:0));

              x=(images[idx].qtab[q+zpos[zzu[zz]]]+1)*cbuf2[cpos_dc+zpos[zzu[zz]]]/(images[idx].qtab[q+zz]+1);
              lcp[6]=(x<0?-1:+1)*(ilog(abs(x)+1)+(x!=0?17:0));
            }
            else
              lcp[4]=lcp[5]=lcp[6]=65535;

            int prev1=0,prev2=0,cnt1=0,cnt2=0,r=0,s=0;
            prev_coef_rs = cbuf[cpos-64];
            for (int i=0; i<acomp; i++) {
              x=0;
              x+=cbuf2[cpos-(acomp-i)*64];
              if (zz==0 && (norst || ls[i]==64)) x-=cbuf2[cpos_dc-(acomp-i)*64-ls[i]];
              if (color[i]==color[acomp]-1) { prev1+=x; cnt1++; r+=cbuf[cpos-(acomp-i)*64]>>4; s+=cbuf[cpos-(acomp-i)*64]&0xF; }
              if (color[acomp]>1 && color[i]==color[0]) { prev2+=x; cnt2++; }
            }
            if (cnt1>0) prev1/=cnt1, r/=cnt1, s/=cnt1, prev_coef_rs=(r<<4)|s;
            if (cnt2>0) prev2/=cnt2;
            prev_coef=(prev1<0?-1:+1)*ilog(11*abs(prev1)+1)+(cnt1<<20);
            prev_coef2=(prev2<0?-1:+1)*ilog(11*abs(prev2)+1);

            if (column==0 && blockW[acomp]>64*acomp) run_pred[1]=run_pred[2], run_pred[0]=0, adv_pred[1]=adv_pred[2], adv_pred[0]=0;
            if (row==0 && blockN[acomp]>64*acomp) run_pred[1]=run_pred[0], run_pred[2]=0, adv_pred[1]=adv_pred[0], adv_pred[2]=0;
          } // !!!!

        }
      }
    }
  }

  // Estimate next bit probability
  if (!images[idx].jpeg || !images[idx].data) return images[idx].next_jpeg;
  if (buf(1+(!bpos))==FF) {
    m.add(128);
    m.set(0, 9);
    m.set(0, 1025);
    m.set(buf(1), 1024);
    return 1;
  }
  if (rstlen>0 && rstlen==column+row*width-rstpos && mcupos==0 && (int)huffcode==(1<<huffbits)-1) {
    m.add(4095);
    m.set(0, 9);
    m.set(0, 1025);
    m.set(buf(1), 1024);
    return 1;
  }

  // Context model
  const int N=32; // size of t, number of contexts
  static BH<9> t(MEM());  // context hash -> bit history
    // As a cache optimization, the context does not include the last 1-2
    // bits of huffcode if the length (huffbits) is not a multiple of 3.
    // The 7 mapped values are for context+{"", 0, 00, 01, 1, 10, 11}.
  static Array<U32> cxt(N);  // context hashes
  static Array<U8*> cp(N);  // context pointers
  static StateMap sm[N];
  static Mixer m1(N+1, 2050, 3);
  static APM a1(0x8000), a2(0x20000);

  // Update model
  if (cp[N-1]) {
    for (int i=0; i<N; ++i)
      *cp[i]=nex(*cp[i],y);
  }
  m1.update();

  // Update context
  const int comp=color[mcupos>>6];
  const int coef=(mcupos&63)|comp<<6;
  const int hc=(huffcode*4+((mcupos&63)==0)*2+(comp==0))|1<<(huffbits+2);
  const bool firstcol=column==0 && blockW[mcupos>>6]>mcupos;
  static int hbcount=2;
  if (++hbcount>2 || huffbits==0) hbcount=0;
  jassert(coef>=0 && coef<256);
  const int zu=zzu[mcupos&63], zv=zzv[mcupos&63];
  if (hbcount==0) {
    int n=hc*32;
    cxt[0]=hash(++n, coef, adv_pred[2]/12+(run_pred[2]<<8), ssum2>>6, prev_coef/72);
    cxt[1]=hash(++n, coef, adv_pred[0]/12+(run_pred[0]<<8), ssum2>>6, prev_coef/72);
    cxt[2]=hash(++n, coef, adv_pred[1]/11+(run_pred[1]<<8), ssum2>>6);
    cxt[3]=hash(++n, rs1, adv_pred[2]/7, run_pred[5]/2, prev_coef/10);
    cxt[4]=hash(++n, rs1, adv_pred[0]/7, run_pred[3]/2, prev_coef/10);
    cxt[5]=hash(++n, rs1, adv_pred[1]/11, run_pred[4]);
    cxt[6]=hash(++n, adv_pred[2]/14, run_pred[2], adv_pred[0]/14, run_pred[0]);
    cxt[7]=hash(++n, cbuf[cpos-blockN[mcupos>>6]]>>4, adv_pred[3]/17, run_pred[1], run_pred[5]);
    cxt[8]=hash(++n, cbuf[cpos-blockW[mcupos>>6]]>>4, adv_pred[3]/17, run_pred[1], run_pred[3]);
    cxt[9]=hash(++n, lcp[0]/22, lcp[1]/22, adv_pred[1]/7, run_pred[1]);
    cxt[10]=hash(++n, lcp[0]/22, lcp[1]/22, mcupos&63, lcp[4]/30);
    cxt[11]=hash(++n, zu/2, lcp[0]/13, lcp[2]/30, prev_coef/40+((prev_coef2/28)<<20));
    cxt[12]=hash(++n, zv/2, lcp[1]/13, lcp[3]/30, prev_coef/40+((prev_coef2/28)<<20));
    cxt[13]=hash(++n, rs1, prev_coef/42, prev_coef2/34, hash(lcp[0]/60,lcp[2]/14,lcp[1]/60,lcp[3]/14));
    cxt[14]=hash(++n, mcupos&63, column>>1);
    cxt[15]=hash(++n, column>>3, min(5+2*(!comp),zu+zv), hash(lcp[0]/10,lcp[2]/40,lcp[1]/10,lcp[3]/40));
    cxt[16]=hash(++n, ssum>>3, mcupos&63);
    cxt[17]=hash(++n, rs1, mcupos&63, run_pred[1]);
    cxt[18]=hash(++n, coef, ssum2>>5, adv_pred[3]/30, (comp)?hash(prev_coef/22,prev_coef2/50):ssum/((mcupos&0x3F)+1));
    cxt[19]=hash(++n, lcp[0]/40, lcp[1]/40, adv_pred[1]/28, hash( (comp)?prev_coef/40+((prev_coef2/40)<<20):lcp[4]/22, min(7,zu+zv), ssum/(2*(zu+zv)+1) ) );
    cxt[20]=hash(++n, zv, cbuf[cpos-blockN[mcupos>>6]], adv_pred[2]/28, run_pred[2]);
    cxt[21]=hash(++n, zu, cbuf[cpos-blockW[mcupos>>6]], adv_pred[0]/28, run_pred[0]);
    cxt[22]=hash(++n, adv_pred[2]/7, run_pred[2]);
    cxt[23]=hash(n, adv_pred[0]/7, run_pred[0]);
    cxt[24]=hash(n, adv_pred[1]/7, run_pred[1]);
    cxt[25]=hash(++n, zv, lcp[1]/14, adv_pred[2]/16, run_pred[5]);
    cxt[26]=hash(++n, zu, lcp[0]/14, adv_pred[0]/16, run_pred[3]);
    cxt[27]=hash(++n, lcp[0]/14, lcp[1]/14, adv_pred[3]/16);
    cxt[28]=hash(++n, coef, prev_coef/10, prev_coef2/20);
    cxt[29]=hash(++n, coef, ssum>>2, prev_coef_rs);
    cxt[30]=hash(++n, coef, adv_pred[1]/17, hash(lcp[(zu<zv)]/24,lcp[2]/20,lcp[3]/24));
    cxt[31]=hash(++n, coef, adv_pred[3]/11, hash(lcp[(zu<zv)]/50,lcp[2+3*(zu*zv>1)]/50,lcp[3+3*(zu*zv>1)]/50));
  }

  // Predict next bit
  m1.add(128);
  jassert(hbcount<=2);
  int p;
 switch(hbcount)
  {
   case 0: for (int i=0; i<N; ++i){ cp[i]=t[cxt[i]]+1, m1.add(p=stretch(sm[i].p(*cp[i]))); m.add(p);} break;
   case 1: { int hc=1+(huffcode&1)*3; for (int i=0; i<N; ++i){ cp[i]+=hc, m1.add(p=stretch(sm[i].p(*cp[i]))); m.add(p); }} break;
   default: { int hc=1+(huffcode&1); for (int i=0; i<N; ++i){ cp[i]+=hc, m1.add(p=stretch(sm[i].p(*cp[i]))); m.add(p); }} break;
  }

  m1.set(firstcol, 2);
  m1.set(coef+256*min(3,huffbits), 1024);
  m1.set((hc&0x1FE)*2+min(3,ilog2(zu+zv)), 1024);
  int pr=m1.p();
  m.add(stretch(pr));
  m.add(pr-2048);
  pr=a1.p(pr, (hc&511)|(((adv_pred[1]/16)&63)<<9), 1023);
  m.add(stretch(pr));
  m.add(pr-2048);
  pr=a2.p(pr, (hc&511)|(coef<<9), 1023);
  m.add(stretch(pr));
  m.add(pr-2048);
  m.set(1 + (zu+zv<5)+(huffbits>8)*2+firstcol*4, 9);
  m.set(1 + (hc&0xFF) + 256*min(3,(zu+zv)/3), 1025);
  m.set(coef+256*min(3,huffbits/2), 1024);
  return 1;
}

//////////////////////////// exeModel /////////////////////////
/*
  Model for x86/x64 code.
  Based on the previous paq* exe models and on DisFilter (http://www.farbrausch.de/~fg/code/disfilter/) by Fabian Giesen.

  Attemps to parse the input stream as x86/x64 instructions, and quantizes them into 32-bit representations that are then
  used as context to predict the next bits, and also extracts possible sparse contexts at several previous positions that
  are relevant to parsing (prefixes, opcode, Mod and R/M fields of the ModRM byte, Scale field of the SIB byte)

  Changelog:
  (18/08/2017) v98: Initial release by Márcio Pais
  (19/08/2017) v99: Bug fixes, tables for instruction categorization, other small improvements
  (29/08/2017) v102: Added variable context dependent on parsing break point
  (22/10/2017) v116: Fixed a bug in function ProcessMode (missing break, thank you Mauro Vezzosi)
*/

// formats
enum InstructionFormat {
  // encoding mode
  fNM = 0x0,      // no ModRM
  fAM = 0x1,      // no ModRM, "address mode" (jumps or direct addresses)
  fMR = 0x2,      // ModRM present
  fMEXTRA = 0x3,  // ModRM present, includes extra bits for opcode
  fMODE = 0x3,    // bitmask for mode

  // no ModRM: size of immediate operand
  fNI = 0x0,      // no immediate
  fBI = 0x4,      // byte immediate
  fWI = 0x8,      // word immediate
  fDI = 0xc,      // dword immediate
  fTYPE = 0xc,    // type mask

  // address mode: type of address operand
  fAD = 0x0,      // absolute address
  fDA = 0x4,      // dword absolute jump target
  fBR = 0x8,      // byte relative jump target
  fDR = 0xc,      // dword relative jump target

  // others
  fERR = 0xf,     // denotes invalid opcodes
};

// 1 byte opcodes
const static U8 Table1[256] = {
  // 0       1       2       3       4       5       6       7       8       9       a       b       c       d       e       f
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI, // 0
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI, // 1
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI, // 2
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI, // 3

  fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI, // 4
  fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI, // 5
  fNM|fNI,fNM|fNI,fMR|fNI,fMR|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fDI,fMR|fDI,fNM|fBI,fMR|fBI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI, // 6
  fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR, // 7

  fMR|fBI,fMR|fDI,fMR|fBI,fMR|fBI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // 8
  fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fAM|fDA,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI, // 9
  fAM|fAD,fAM|fAD,fAM|fAD,fAM|fAD,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fBI,fNM|fDI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI, // a
  fNM|fBI,fNM|fBI,fNM|fBI,fNM|fBI,fNM|fBI,fNM|fBI,fNM|fBI,fNM|fBI,fNM|fDI,fNM|fDI,fNM|fDI,fNM|fDI,fNM|fDI,fNM|fDI,fNM|fDI,fNM|fDI, // b

  fMR|fBI,fMR|fBI,fNM|fWI,fNM|fNI,fMR|fNI,fMR|fNI,fMR|fBI,fMR|fDI,fNM|fBI,fNM|fNI,fNM|fWI,fNM|fNI,fNM|fNI,fNM|fBI,fERR   ,fNM|fNI, // c
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fBI,fNM|fBI,fNM|fNI,fNM|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // d
  fAM|fBR,fAM|fBR,fAM|fBR,fAM|fBR,fNM|fBI,fNM|fBI,fNM|fBI,fNM|fBI,fAM|fDR,fAM|fDR,fAM|fAD,fAM|fBR,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI, // e
  fNM|fNI,fERR   ,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fMEXTRA,fMEXTRA,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fMEXTRA,fMEXTRA, // f
};

// 2 byte opcodes
const static U8 Table2[256] = {
  // 0       1       2       3       4       5       6       7       8       9       a       b       c       d       e       f
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fNM|fNI,fERR   ,fNM|fNI,fNM|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 0
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 1
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // 2
  fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fERR   ,fNM|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 3

  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // 4
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // 5
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // 6
  fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fMR|fNI,fMR|fNI, // 7

  fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR,fAM|fDR, // 8
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // 9
  fNM|fNI,fNM|fNI,fNM|fNI,fMR|fNI,fMR|fBI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fMR|fNI,fMR|fBI,fMR|fNI,fERR   ,fMR|fNI, // a
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // b

  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI,fNM|fNI, // c
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // d
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // e
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   , // f
};

// 3 byte opcodes 0F38XX
const static U8 Table3_38[256] = {
  // 0       1       2       3       4       5       6       7       8       9       a       b       c       d       e       f
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   , // 0
  fMR|fNI,fERR   ,fERR   ,fERR   ,fMR|fNI,fMR|fNI,fERR   ,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fMR|fNI,fMR|fNI,fMR|fNI,fERR   , // 1
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fERR   ,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   , // 2
  fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // 3
  fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 4
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 5
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 6
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 7
  fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 8
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 9
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // a
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // b
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // c
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // d
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // e
  fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // f
};

// 3 byte opcodes 0F3AXX
const static U8 Table3_3A[256] = {
  // 0       1       2       3       4       5       6       7       8       9       a       b       c       d       e       f
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI, // 0
  fERR   ,fERR   ,fERR   ,fERR   ,fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 1
  fMR|fBI,fMR|fBI,fMR|fBI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 2
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 3
  fMR|fBI,fMR|fBI,fMR|fBI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 4
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 5
  fMR|fBI,fMR|fBI,fMR|fBI,fMR|fBI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 6
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 7
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 8
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // 9
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // a
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // b
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // c
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // d
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // e
  fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // f
};

// escape opcodes using ModRM byte to get more variants
const static U8 TableX[32] = {
  // 0       1       2       3       4       5       6       7
  fMR|fBI,fERR   ,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // escapes for 0xf6
  fMR|fDI,fERR   ,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI,fMR|fNI, // escapes for 0xf7
  fMR|fNI,fMR|fNI,fERR   ,fERR   ,fERR   ,fERR   ,fERR   ,fERR   , // escapes for 0xfe
  fMR|fNI,fMR|fNI,fMR|fNI,fERR   ,fMR|fNI,fERR   ,fMR|fNI,fERR   , // escapes for 0xff
};

const static U8 InvalidX64Ops[19] = {0x06, 0x07, 0x16, 0x17, 0x1E, 0x1F, 0x27, 0x2F, 0x37, 0x3F, 0x60, 0x61, 0x62, 0x82, 0x9A, 0xD4, 0xD5, 0xD6, 0xEA,};
const static U8 X64Prefixes[8] = {0x26, 0x2E, 0x36, 0x3E, 0x9B, 0xF0, 0xF2, 0xF3,};

enum InstructionCategory {
  OP_INVALID              =  0,
  OP_PREFIX_SEGREG        =  1,
  OP_PREFIX               =  2,
  OP_PREFIX_X87FPU        =  3,
  OP_GEN_DATAMOV          =  4,
  OP_GEN_STACK            =  5,
  OP_GEN_CONVERSION       =  6,
  OP_GEN_ARITH_DECIMAL    =  7,
  OP_GEN_ARITH_BINARY     =  8,
  OP_GEN_LOGICAL          =  9,
  OP_GEN_SHF_ROT          = 10,
  OP_GEN_BIT              = 11,
  OP_GEN_BRANCH           = 12,
  OP_GEN_BRANCH_COND      = 13,
  OP_GEN_BREAK            = 14,
  OP_GEN_STRING           = 15,
  OP_GEN_INOUT            = 16,
  OP_GEN_FLAG_CONTROL     = 17,
  OP_GEN_SEGREG           = 18,
  OP_GEN_CONTROL          = 19,
  OP_SYSTEM               = 20,
  OP_X87_DATAMOV          = 21,
  OP_X87_ARITH            = 22,
  OP_X87_COMPARISON       = 23,
  OP_X87_TRANSCENDENTAL   = 24,
  OP_X87_LOAD_CONSTANT    = 25,
  OP_X87_CONTROL          = 26,
  OP_X87_CONVERSION       = 27,
  OP_STATE_MANAGEMENT     = 28,
  OP_MMX                  = 29,
  OP_SSE                  = 30,
  OP_SSE_DATAMOV          = 31,
};

const static U8 TypeOp1[256] = {
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //03
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_STACK         , OP_GEN_STACK         , //07
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , //0B
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_STACK         , OP_PREFIX            , //0F
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //13
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_STACK         , OP_GEN_STACK         , //17
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //1B
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_STACK         , OP_GEN_STACK         , //1F
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , //23
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_PREFIX_SEGREG     , OP_GEN_ARITH_DECIMAL , //27
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //2B
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_PREFIX_SEGREG     , OP_GEN_ARITH_DECIMAL , //2F
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , //33
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_PREFIX_SEGREG     , OP_GEN_ARITH_DECIMAL , //37
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //3B
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_PREFIX_SEGREG     , OP_GEN_ARITH_DECIMAL , //3F
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //43
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //47
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //4B
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //4F
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , //53
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , //57
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , //5B
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_STACK         , //5F
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_BREAK         , OP_GEN_CONVERSION    , //63
  OP_PREFIX_SEGREG     , OP_PREFIX_SEGREG     , OP_PREFIX            , OP_PREFIX            , //67
  OP_GEN_STACK         , OP_GEN_ARITH_BINARY  , OP_GEN_STACK         , OP_GEN_ARITH_BINARY  , //6B
  OP_GEN_INOUT         , OP_GEN_INOUT         , OP_GEN_INOUT         , OP_GEN_INOUT         , //6F
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //73
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //77
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //7B
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //7F
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //83
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //87
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //8B
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_STACK         , //8F
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //93
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //97
  OP_GEN_CONVERSION    , OP_GEN_CONVERSION    , OP_GEN_BRANCH        , OP_PREFIX_X87FPU     , //9B
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //9F
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //A3
  OP_GEN_STRING        , OP_GEN_STRING        , OP_GEN_STRING        , OP_GEN_STRING        , //A7
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_STRING        , OP_GEN_STRING        , //AB
  OP_GEN_STRING        , OP_GEN_STRING        , OP_GEN_STRING        , OP_GEN_STRING        , //AF
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //B3
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //B7
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //BB
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //BF
  OP_GEN_SHF_ROT       , OP_GEN_SHF_ROT       , OP_GEN_BRANCH        , OP_GEN_BRANCH        , //C3
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //C7
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_BRANCH        , OP_GEN_BRANCH        , //CB
  OP_GEN_BREAK         , OP_GEN_BREAK         , OP_GEN_BREAK         , OP_GEN_BREAK         , //CF
  OP_GEN_SHF_ROT       , OP_GEN_SHF_ROT       , OP_GEN_SHF_ROT       , OP_GEN_SHF_ROT       , //D3
  OP_GEN_ARITH_DECIMAL , OP_GEN_ARITH_DECIMAL , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //D7
  OP_X87_ARITH         , OP_X87_DATAMOV       , OP_X87_ARITH         , OP_X87_DATAMOV       , //DB
  OP_X87_ARITH         , OP_X87_DATAMOV       , OP_X87_ARITH         , OP_X87_DATAMOV       , //DF
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //E3
  OP_GEN_INOUT         , OP_GEN_INOUT         , OP_GEN_INOUT         , OP_GEN_INOUT         , //E7
  OP_GEN_BRANCH        , OP_GEN_BRANCH        , OP_GEN_BRANCH        , OP_GEN_BRANCH        , //EB
  OP_GEN_INOUT         , OP_GEN_INOUT         , OP_GEN_INOUT         , OP_GEN_INOUT         , //EF
  OP_PREFIX            , OP_GEN_BREAK         , OP_PREFIX            , OP_PREFIX            , //F3
  OP_SYSTEM            , OP_GEN_FLAG_CONTROL  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , //F7
  OP_GEN_FLAG_CONTROL  , OP_GEN_FLAG_CONTROL  , OP_GEN_FLAG_CONTROL  , OP_GEN_FLAG_CONTROL  , //FB
  OP_GEN_FLAG_CONTROL  , OP_GEN_FLAG_CONTROL  , OP_GEN_ARITH_BINARY  , OP_GEN_BRANCH        , //FF
};

const static U8 TypeOp2[256] = {
  OP_SYSTEM            , OP_SYSTEM            , OP_SYSTEM            , OP_SYSTEM            , //03
  OP_INVALID           , OP_SYSTEM            , OP_SYSTEM            , OP_SYSTEM            , //07
  OP_SYSTEM            , OP_SYSTEM            , OP_INVALID           , OP_GEN_CONTROL       , //0B
  OP_INVALID           , OP_GEN_CONTROL       , OP_INVALID           , OP_INVALID           , //0F
  OP_SSE_DATAMOV       , OP_SSE_DATAMOV       , OP_SSE_DATAMOV       , OP_SSE_DATAMOV       , //13
  OP_SSE               , OP_SSE               , OP_SSE_DATAMOV       , OP_SSE_DATAMOV       , //17
  OP_SSE               , OP_GEN_CONTROL       , OP_GEN_CONTROL       , OP_GEN_CONTROL       , //1B
  OP_GEN_CONTROL       , OP_GEN_CONTROL       , OP_GEN_CONTROL       , OP_GEN_CONTROL       , //1F
  OP_SYSTEM            , OP_SYSTEM            , OP_SYSTEM            , OP_SYSTEM            , //23
  OP_SYSTEM            , OP_INVALID           , OP_SYSTEM            , OP_INVALID           , //27
  OP_SSE_DATAMOV       , OP_SSE_DATAMOV       , OP_SSE               , OP_SSE               , //2B
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //2F
  OP_SYSTEM            , OP_SYSTEM            , OP_SYSTEM            , OP_SYSTEM            , //33
  OP_SYSTEM            , OP_SYSTEM            , OP_INVALID           , OP_INVALID           , //37
  OP_PREFIX            , OP_INVALID           , OP_PREFIX            , OP_INVALID           , //3B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //3F
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //43
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //47
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //4B
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //4F
  OP_SSE_DATAMOV       , OP_SSE               , OP_SSE               , OP_SSE               , //53
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //57
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //5B
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //5F
  OP_MMX               , OP_MMX               , OP_MMX               , OP_MMX               , //63
  OP_MMX               , OP_MMX               , OP_MMX               , OP_MMX               , //67
  OP_MMX               , OP_MMX               , OP_MMX               , OP_MMX               , //6B
  OP_INVALID           , OP_INVALID           , OP_MMX               , OP_MMX               , //6F
  OP_SSE               , OP_MMX               , OP_MMX               , OP_MMX               , //73
  OP_MMX               , OP_MMX               , OP_MMX               , OP_MMX               , //77
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //7B
  OP_INVALID           , OP_INVALID           , OP_MMX               , OP_MMX               , //7F
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //83
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //87
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //8B
  OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , OP_GEN_BRANCH_COND   , //8F
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //93
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //97
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //9B
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //9F
  OP_GEN_STACK         , OP_GEN_STACK         , OP_GEN_CONTROL       , OP_GEN_BIT           , //A3
  OP_GEN_SHF_ROT       , OP_GEN_SHF_ROT       , OP_INVALID           , OP_INVALID           , //A7
  OP_GEN_STACK         , OP_GEN_STACK         , OP_SYSTEM            , OP_GEN_BIT           , //AB
  OP_GEN_SHF_ROT       , OP_GEN_SHF_ROT       , OP_STATE_MANAGEMENT  , OP_GEN_ARITH_BINARY  , //AF
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_BIT           , //B3
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_CONVERSION    , OP_GEN_CONVERSION    , //B7
  OP_INVALID           , OP_GEN_CONTROL       , OP_GEN_BIT           , OP_GEN_BIT           , //BB
  OP_GEN_BIT           , OP_GEN_BIT           , OP_GEN_CONVERSION    , OP_GEN_CONVERSION    , //BF
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_SSE               , OP_SSE               , //C3
  OP_SSE               , OP_SSE               , OP_SSE               , OP_GEN_DATAMOV       , //C7
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //CB
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , //CF
  OP_INVALID           , OP_MMX               , OP_MMX               , OP_MMX               , //D3
  OP_SSE               , OP_MMX               , OP_INVALID           , OP_SSE               , //D7
  OP_MMX               , OP_MMX               , OP_SSE               , OP_MMX               , //DB
  OP_MMX               , OP_MMX               , OP_SSE               , OP_MMX               , //DF
  OP_SSE               , OP_MMX               , OP_SSE               , OP_MMX               , //E3
  OP_SSE               , OP_MMX               , OP_INVALID           , OP_SSE               , //E7
  OP_MMX               , OP_MMX               , OP_SSE               , OP_MMX               , //EB
  OP_MMX               , OP_MMX               , OP_SSE               , OP_MMX               , //EF
  OP_INVALID           , OP_MMX               , OP_MMX               , OP_MMX               , //F3
  OP_SSE               , OP_MMX               , OP_SSE               , OP_SSE               , //F7
  OP_MMX               , OP_MMX               , OP_MMX               , OP_SSE               , //FB
  OP_MMX               , OP_MMX               , OP_MMX               , OP_INVALID           , //FF
};

const static U8 TypeOp3_38[256] = {
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //03
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //07
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //0B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //0F
  OP_SSE               , OP_INVALID           , OP_INVALID           , OP_INVALID           , //13
  OP_SSE               , OP_SSE               , OP_INVALID           , OP_SSE               , //17
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //1B
  OP_SSE               , OP_SSE               , OP_SSE               , OP_INVALID           , //1F
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //23
  OP_SSE               , OP_SSE               , OP_INVALID           , OP_INVALID           , //27
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //2B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //2F
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //33
  OP_SSE               , OP_SSE               , OP_INVALID           , OP_SSE               , //37
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //3B
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //3F
  OP_SSE               , OP_SSE               , OP_INVALID           , OP_INVALID           , //43
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //47
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //4B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //4F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //53
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //57
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //5B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //5F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //63
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //67
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //6B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //6F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //73
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //77
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //7B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //7F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //83
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //87
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //8B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //8F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //93
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //97
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //9B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //9F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //A3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //A7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //AB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //AF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //B3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //B7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //BB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //BF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //C3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //C7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //CB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //CF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //D3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //D7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //DB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //DF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //E3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //E7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //EB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //EF
  OP_GEN_DATAMOV       , OP_GEN_DATAMOV       , OP_INVALID           , OP_INVALID           , //F3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //F7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //FB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //FF
};

const static U8 TypeOp3_3A[256] = {
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //03
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //07
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //0B
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //0F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //13
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //17
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //1B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //1F
  OP_SSE               , OP_SSE               , OP_SSE               , OP_INVALID           , //23
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //27
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //2B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //2F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //33
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //37
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //3B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //3F
  OP_SSE               , OP_SSE               , OP_SSE               , OP_INVALID           , //43
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //47
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //4B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //4F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //53
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //57
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //5B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //5F
  OP_SSE               , OP_SSE               , OP_SSE               , OP_SSE               , //63
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //67
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //6B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //6F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //73
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //77
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //7B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //7F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //83
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //87
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //8B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //8F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //93
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //97
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //9B
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //9F
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //A3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //A7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //AB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //AF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //B3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //B7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //BB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //BF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //C3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //C7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //CB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //CF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //D3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //D7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //DB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //DF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //E3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //E7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //EB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //EF
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //F3
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //F7
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //FB
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           , //FF
};

const static U8 TypeOpX[32] = {
  // escapes for F6
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_ARITH_BINARY  ,
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  ,
  // escapes for F7
  OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_LOGICAL       , OP_GEN_ARITH_BINARY  ,
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  ,
  // escapes for FE
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_INVALID           , OP_INVALID           ,
  OP_INVALID           , OP_INVALID           , OP_INVALID           , OP_INVALID           ,
  // escapes for FF
  OP_GEN_ARITH_BINARY  , OP_GEN_ARITH_BINARY  , OP_GEN_BRANCH        , OP_GEN_BRANCH        ,
  OP_GEN_BRANCH        , OP_GEN_BRANCH        , OP_GEN_STACK         , OP_INVALID           ,
};

enum Prefixes {
  ES_OVERRIDE = 0x26,
  CS_OVERRIDE = 0x2E,
  SS_OVERRIDE = 0x36,
  DS_OVERRIDE = 0x3E,
  FS_OVERRIDE = 0x64,
  GS_OVERRIDE = 0x65,
  AD_OVERRIDE = 0x67,
  WAIT_FPU    = 0x9B,
  LOCK        = 0xF0,
  REP_N_STR   = 0xF2,
  REP_STR     = 0xF3,
};

enum Opcodes {
  // 1-byte opcodes of special interest (for one reason or another)
  OP_2BYTE  = 0x0f,     // start of 2-byte opcode
  OP_OSIZE  = 0x66,     // operand size prefix
  OP_CALLF  = 0x9a,
  OP_RETNI  = 0xc2,     // ret near+immediate
  OP_RETN   = 0xc3,
  OP_ENTER  = 0xc8,
  OP_INT3   = 0xcc,
  OP_INTO   = 0xce,
  OP_CALLN  = 0xe8,
  OP_JMPF   = 0xea,
  OP_ICEBP  = 0xf1,
};

enum ExeState {
  Start             =  0,
  Pref_Op_Size      =  1,
  Pref_MultiByte_Op =  2,
  ParseFlags        =  3,
  ExtraFlags        =  4,
  ReadModRM         =  5,
  Read_OP3_38       =  6,
  Read_OP3_3A       =  7,
  ReadSIB           =  8,
  Read8             =  9,
  Read16            = 10,
  Read32            = 11,
  Read8_ModRM       = 12,
  Read16_f          = 13,
  Read32_ModRM      = 14,
  Error             = 15,
};

struct OpCache {
  U32 Op[CacheSize];
  U32 Index;
};

struct Instruction {
  U32 Data;
  U8 Prefix, Code, ModRM, SIB, REX, Flags, BytesRead, Size, Category;
  bool MustCheckREX, Decoding, o16, imm8;
};

#define CodeShift            3
#define CodeMask             (0xFF<<CodeShift)
#define ClearCodeMask        ((-1)^CodeMask)
#define PrefixMask           ((1<<CodeShift)-1)
#define OperandSizeOverride  (0x01<<(8+CodeShift))
#define MultiByteOpcode      (0x02<<(8+CodeShift))
#define PrefixREX            (0x04<<(8+CodeShift))
#define Prefix38             (0x08<<(8+CodeShift))
#define Prefix3A             (0x10<<(8+CodeShift))
#define HasExtraFlags        (0x20<<(8+CodeShift))
#define HasModRM             (0x40<<(8+CodeShift))
#define ModRMShift           (7+8+CodeShift)
#define SIBScaleShift        (ModRMShift+8-6)
#define RegDWordDisplacement (0x01<<(8+SIBScaleShift))
#define AddressMode          (0x02<<(8+SIBScaleShift))
#define TypeShift            (2+8+SIBScaleShift)
#define CategoryShift        5
#define CategoryMask         ((1<<CategoryShift)-1)
#define ModRM_mod            0xC0
#define ModRM_reg            0x38
#define ModRM_rm             0x07
#define SIB_scale            0xC0
#define SIB_index            0x38
#define SIB_base             0x07
#define REX_w                0x08

#define MinRequired          8 // minimum required consecutive valid instructions to be considered as code

inline bool IsInvalidX64Op(U8 Op){
  for (int i=0; i<19; i++){
    if (Op == InvalidX64Ops[i])
      return true;
  }
  return false;
}

inline bool IsValidX64Prefix(U8 Prefix){
  for (int i=0; i<8; i++){
    if (Prefix == X64Prefixes[i])
      return true;
  }
  return ((Prefix>=0x40 && Prefix<=0x4F) || (Prefix>=0x64 && Prefix<=0x67));
}

void ProcessMode(Instruction &Op, ExeState &State){
  if ((Op.Flags&fMODE)==fAM){
    Op.Data|=AddressMode;
    Op.BytesRead = 0;
    switch (Op.Flags&fTYPE){
      case fDR : Op.Data|=(2<<TypeShift);
      case fDA : Op.Data|=(1<<TypeShift);
      case fAD : {
        State = Read32;
        break;
      }
      case fBR : {
        Op.Data|=(2<<TypeShift);
        State = Read8;
      }
    }
  }
  else{
    switch (Op.Flags&fTYPE){
      case fBI : State = Read8; break;
      case fWI : {
        State = Read16;
        Op.Data|=(1<<TypeShift);
        Op.BytesRead = 0;
        break;
      }
      case fDI : {
        // x64 Move with 8byte immediate? [REX.W is set, opcodes 0xB8+r]
        Op.imm8=((Op.REX & REX_w)>0 && (Op.Code&0xF8)==0xB8);
        if (!Op.o16 || Op.imm8){
          State = Read32;
          Op.Data|=(2<<TypeShift);
        }
        else{
          State = Read16;
          Op.Data|=(3<<TypeShift);
        }
        Op.BytesRead = 0;
        break;
      }
      default: State = Start; /*no immediate*/
    }
  }
}

void ProcessFlags2(Instruction &Op, ExeState &State){
  //if arriving from state ExtraFlags, we've already read the ModRM byte
  if ((Op.Flags&fMODE)==fMR && State!=ExtraFlags){
    State = ReadModRM;
    return;
  }
  ProcessMode(Op, State);
}

void ProcessFlags(Instruction &Op, ExeState &State){
  if (Op.Code==OP_CALLF || Op.Code==OP_JMPF || Op.Code==OP_ENTER){
    Op.BytesRead = 0;
    State = Read16_f;
    return; //must exit, ENTER has ModRM too
  }
  ProcessFlags2(Op, State);
}

void CheckFlags(Instruction &Op, ExeState &State){
  //must peek at ModRM byte to read the REG part, so we can know the opcode
  if (Op.Flags==fMEXTRA)
    State = ExtraFlags;
  else if (Op.Flags==fERR){
    memset(&Op, 0, sizeof(Instruction));
    State = Error;
  }
  else
    ProcessFlags(Op, State);
}

void ReadFlags(Instruction &Op, ExeState &State){
  Op.Flags = Table1[Op.Code];
  Op.Category = TypeOp1[Op.Code];
  CheckFlags(Op, State);
}

void ProcessModRM(Instruction &Op, ExeState &State){
  if ((Op.ModRM & ModRM_mod)==0x40)
    State = Read8_ModRM; //register+byte displacement
  else if ((Op.ModRM & ModRM_mod)==0x80 || (Op.ModRM & (ModRM_mod|ModRM_rm))==0x05 || (Op.ModRM<0x40 && (Op.SIB & SIB_base)==0x05) ){
    State = Read32_ModRM; //register+dword displacement
    Op.BytesRead = 0;
  }
  else
    ProcessMode(Op, State);
}

void ApplyCodeAndSetFlag(Instruction &Op, U32 Flag = 0){
  Op.Data&=ClearCodeMask; \
  Op.Data|=(Op.Code<<CodeShift)|Flag;
}

inline U32 OpN(OpCache &Cache, U32 n){
  return Cache.Op[ (Cache.Index-n)&(CacheSize-1) ];
}

inline U32 OpNCateg(U32 &Mask, U32 n){
  return ((Mask>>(CategoryShift*(n-1)))&CategoryMask);
}

inline int pref(int i) { return (buf(i)==0x0f)+2*(buf(i)==0x66)+3*(buf(i)==0x67); }

// Get context at buf(i) relevant to parsing 32-bit x86 code
U32 execxt(int i, int x=0) {
  int prefix=0, opcode=0, modrm=0, sib=0;
  if (i) prefix+=4*pref(i--);
  if (i) prefix+=pref(i--);
  if (i) opcode+=buf(i--);
  if (i) modrm+=buf(i--)&(ModRM_mod|ModRM_rm);
  if (i&&((modrm&ModRM_rm)==4)&&(modrm<ModRM_mod)) sib=buf(i)&SIB_scale;
  return prefix|opcode<<4|modrm<<12|x<<20|sib<<(28-6);
}

bool exeModel(Mixer& m, bool Forced = false, ModelStats *Stats = NULL) {
  const int N1=9, N2=10;
  static ContextMap2 cm(MEM()*2, N1+N2);
  static OpCache Cache;
  static U32 StateBH[256];
  static ExeState pState = Start, State = Start;
  static Instruction Op;
  static U32 TotalOps = 0, OpMask = 0, OpCategMask = 0, Context = 0, BrkPoint = 0, BrkCtx = 0;
  static bool Valid = false;
  if (!bpos) {
    pState = State;
    U8 B = (U8)c4;
    Op.Size++;
    switch (State){
      case Start: case Error: {
        // previous code may have just been a REX prefix
        bool Skip = false;
        if (Op.MustCheckREX){
          Op.MustCheckREX = false;
          // valid x64 code?
          if (!IsInvalidX64Op(B) && !IsValidX64Prefix(B)){
            Op.REX = Op.Code;
            Op.Code = B;
            Op.Data = PrefixREX|(Op.Code<<CodeShift)|(Op.Data&PrefixMask);
            Skip = true;
          }
        }

        Op.ModRM = Op.SIB = Op.REX = Op.Flags = Op.BytesRead = 0;
        if (!Skip){
          Op.Code = B;
          // possible REX prefix?
          Op.MustCheckREX = ((Op.Code&0xF0)==0x40) && (!(Op.Decoding && ((Op.Data&PrefixMask)==1)));

          // check prefixes
          Op.Prefix = (Op.Code==ES_OVERRIDE || Op.Code==CS_OVERRIDE || Op.Code==SS_OVERRIDE || Op.Code==DS_OVERRIDE) + //invalid in x64
                      (Op.Code==FS_OVERRIDE)*2 +
                      (Op.Code==GS_OVERRIDE)*3 +
                      (Op.Code==AD_OVERRIDE)*4 +
                      (Op.Code==WAIT_FPU)*5 +
                      (Op.Code==LOCK)*6 +
                      (Op.Code==REP_N_STR || Op.Code==REP_STR)*7;

          if (!Op.Decoding){
            TotalOps+=(Op.Data!=0)-(Cache.Index && Cache.Op[ Cache.Index&(CacheSize-1) ]!=0);
            OpMask = (OpMask<<1)|(State!=Error);
            OpCategMask = (OpCategMask<<CategoryShift)|(Op.Category);
            Op.Size = 0;

            Cache.Op[ Cache.Index&(CacheSize-1) ] = Op.Data;
            Cache.Index++;

            if (!Op.Prefix)
              Op.Data = Op.Code<<CodeShift;
            else{
              Op.Data = Op.Prefix;
              Op.Category = TypeOp1[Op.Code];
              Op.Decoding = true;
              BrkCtx = hash(1+(BrkPoint = 0), Op.Prefix, OpCategMask&CategoryMask);
              break;
            }
          }
          else{
            // we only have enough bits for one prefix, so the
            // instruction will be encoded with the last one
            if (!Op.Prefix){
              Op.Data|=(Op.Code<<CodeShift);
              Op.Decoding = false;
            }
            else{
              Op.Data = Op.Prefix;
              Op.Category = TypeOp1[Op.Code];
              BrkCtx = hash(1+(BrkPoint = 1), Op.Prefix, OpCategMask&CategoryMask);
              break;
            }
          }
        }

        if ((Op.o16=(Op.Code==OP_OSIZE)))
          State = Pref_Op_Size;
        else if (Op.Code==OP_2BYTE)
          State = Pref_MultiByte_Op;
        else
          ReadFlags(Op, State);
        BrkCtx = hash(1+(BrkPoint = 2), State, Op.Code, (OpCategMask&CategoryMask), OpN(Cache,1)&((ModRM_mod|ModRM_reg|ModRM_rm)<<ModRMShift));
        break;
      }
      case Pref_Op_Size : {
        Op.Code = B;
        ApplyCodeAndSetFlag(Op, OperandSizeOverride);
        ReadFlags(Op, State);
        BrkCtx = hash(1+(BrkPoint = 3), State);
        break;
      }
      case Pref_MultiByte_Op : {
        Op.Code = B;
        Op.Data|=MultiByteOpcode;

        if (Op.Code==0x38)
          State = Read_OP3_38;
        else if (Op.Code==0x3A)
          State = Read_OP3_3A;
        else{
          ApplyCodeAndSetFlag(Op);
          Op.Flags = Table2[Op.Code];
          Op.Category = TypeOp2[Op.Code];
          CheckFlags(Op, State);
        }
        BrkCtx = hash(1+(BrkPoint = 4), State);
        break;
      }
      case ParseFlags : {
        ProcessFlags(Op, State);
        BrkCtx = hash(1+(BrkPoint = 5), State);
        break;
      }
      case ExtraFlags : case ReadModRM : {
        Op.ModRM = B;
        Op.Data|=(Op.ModRM<<ModRMShift)|HasModRM;
        Op.SIB = 0;
        if (Op.Flags==fMEXTRA){
          Op.Data|=HasExtraFlags;
          int i = ((Op.ModRM>>3)&0x07) | ((Op.Code&0x01)<<3) | ((Op.Code&0x08)<<1);
          Op.Flags = TableX[i];
          Op.Category = TypeOpX[i];
          if (Op.Flags==fERR){
            memset(&Op, 0, sizeof(Instruction));
            State = Error;
            BrkCtx = hash(1+(BrkPoint = 6), State);
            break;
          }
          ProcessFlags(Op, State);
          BrkCtx = hash(1+(BrkPoint = 7), State);
          break;
        }

        if ((Op.ModRM & ModRM_rm)==4 && Op.ModRM<ModRM_mod){
          State = ReadSIB;
          BrkCtx = hash(1+(BrkPoint = 8), State);
          break;
        }

        ProcessModRM(Op, State);
        BrkCtx = hash(1+(BrkPoint = 9), State, Op.Code );
        break;
      }
      case Read_OP3_38 : case Read_OP3_3A : {
        Op.Code = B;
        ApplyCodeAndSetFlag(Op, Prefix38<<(State-Read_OP3_38));
        if (State==Read_OP3_38){
          Op.Flags = Table3_38[Op.Code];
          Op.Category = TypeOp3_38[Op.Code];
        }
        else{
          Op.Flags = Table3_3A[Op.Code];
          Op.Category = TypeOp3_3A[Op.Code];
        }
        CheckFlags(Op, State);
        BrkCtx = hash(1+(BrkPoint = 10), State);
        break;
      }
      case ReadSIB : {
        Op.SIB = B;
        Op.Data|=((Op.SIB & SIB_scale)<<SIBScaleShift);
        ProcessModRM(Op, State);
        BrkCtx = hash(1+(BrkPoint = 11), State, Op.SIB&SIB_scale);
        break;
      }
      case Read8 : case Read16 : case Read32 : {
        if (++Op.BytesRead>=(2*(State-Read8)<<Op.imm8)){
          Op.BytesRead = 0;
          Op.imm8 = false;
          State = Start;
        }
        BrkCtx = hash(1+(BrkPoint = 12), State, Op.Flags&fMODE, Op.BytesRead, ((Op.BytesRead>1)?(buf(Op.BytesRead)<<8):0)|((Op.BytesRead)?B:0) );
        break;
      }
      case Read8_ModRM : {
        ProcessMode(Op, State);
        BrkCtx = hash(1+(BrkPoint = 13), State);
        break;
      }
      case Read16_f : {
        if (++Op.BytesRead==2){
          Op.BytesRead = 0;
          ProcessFlags2(Op, State);
        }
        BrkCtx = hash(1+(BrkPoint = 14), State);
        break;
      }
      case Read32_ModRM : {
        Op.Data|=RegDWordDisplacement;
        if (++Op.BytesRead==4){
          Op.BytesRead = 0;
          ProcessMode(Op, State);
        }
        BrkCtx = hash(1+(BrkPoint = 15), State);
        break;
      }
    }

    Valid = (TotalOps>2*MinRequired) && ((OpMask&((1<<MinRequired)-1))==((1<<MinRequired)-1));
    Context = State+16*Op.BytesRead+16*(Op.REX & REX_w);
    StateBH[Context] = (StateBH[Context]<<8)|B;

    if (Valid || Forced){
      int mask=0, count0=0;
      for (int i=0, j=0; i<N1; ++i){
        if (i>1) mask=mask*2+(buf(i-1)==0), count0+=mask&1;
        j=(i<4)?i+1:5+(i-4)*(2+(i>6));
        cm.set(hash(execxt(j, buf(1)*(j>6)), ((1<<N1)|mask)*(count0*N1/2>=i), (0x08|(blpos&0x07))*(i<4)));
      }

      cm.set(BrkCtx);

      mask = PrefixMask|(0xF8<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A;
      cm.set(hash(OpN(Cache, 1)&(mask|RegDWordDisplacement|AddressMode), State+16*Op.BytesRead, Op.Data&mask, Op.REX, Op.Category));

      mask = 0x04|(0xFE<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A|((ModRM_mod|ModRM_reg)<<ModRMShift);
      cm.set(hash(
        OpN(Cache, 1)&mask, OpN(Cache, 2)&mask, OpN(Cache, 3)&mask,
        Context+256*((Op.ModRM & ModRM_mod)==ModRM_mod),
        Op.Data&((mask|PrefixREX)^(ModRM_mod<<ModRMShift))
      ));

      mask = 0x04|CodeMask;
      cm.set(hash(OpN(Cache, 1)&mask, OpN(Cache, 2)&mask, OpN(Cache, 3)&mask, OpN(Cache, 4)&mask, (Op.Data&mask)|(State<<11)|(Op.BytesRead<<15)));

      mask = 0x04|(0xFC<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A;
      cm.set(hash(State+16*Op.BytesRead, Op.Data&mask, Op.Category*8 + (OpMask&0x07), Op.Flags, ((Op.SIB & SIB_base)==5)*4+((Op.ModRM & ModRM_reg)==ModRM_reg)*2+((Op.ModRM & ModRM_mod)==0)));

      mask = PrefixMask|CodeMask|OperandSizeOverride|MultiByteOpcode|PrefixREX|Prefix38|Prefix3A|HasExtraFlags|HasModRM|((ModRM_mod|ModRM_rm)<<ModRMShift);
      cm.set(hash(Op.Data&mask, State+16*Op.BytesRead, Op.Flags));

      mask = PrefixMask|CodeMask|OperandSizeOverride|MultiByteOpcode|Prefix38|Prefix3A|HasExtraFlags|HasModRM;
      cm.set(hash(OpN(Cache, 1)&mask, State, Op.BytesRead*2+((Op.REX&REX_w)>0), Op.Data&((U16)(mask^OperandSizeOverride))));

      mask = 0x04|(0xFE<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A|(ModRM_reg<<ModRMShift);
      cm.set(hash(OpN(Cache, 1)&mask, OpN(Cache, 2)&mask, State+16*Op.BytesRead, Op.Data&(mask|PrefixMask|CodeMask)));

      cm.set(State+16*Op.BytesRead);

      cm.set(hash(
        (0x100|B)*(Op.BytesRead>0),
        State+16*pState+256*Op.BytesRead,
        ((Op.Flags&fMODE)==fAM)*16 + (Op.REX & REX_w) + (Op.o16)*4 + ((Op.Code & 0xFE)==0xE8)*2 + ((Op.Data & MultiByteOpcode)!=0 && (Op.Code & 0xF0)==0x80)
      ));
    }
  }

  if (Valid || Forced)
    cm.mix(m);
  else{
      for (int i=0; i<(N1+N2)*8; ++i)
        m.add(0);
  }
  U8 s = ((StateBH[Context]>>(28-bpos))&0x08) |
         ((StateBH[Context]>>(21-bpos))&0x04) |
         ((StateBH[Context]>>(14-bpos))&0x02) |
         ((StateBH[Context]>>( 7-bpos))&0x01) |
         ((Op.Category==OP_GEN_BRANCH)<<4)|
         (((c0&((1<<bpos)-1))==0)<<5);

  m.set(Context*4+(s>>4), 1024);
  m.set(State*64+bpos*8+(Op.BytesRead>0)*4+(s>>4), 1024);
  m.set((BrkCtx&0x1FF)|((s&0x20)<<4), 1024);
  m.set(hash(Op.Code, State, OpN(Cache, 1)&CodeMask)&0x1FFF, 8192);
  m.set(hash(State, bpos, Op.Code, Op.BytesRead)&0x1FFF, 8192);

  if (Stats)
    (*Stats).x86_64 = Valid|(Context<<1)|(s<<9);
  return Valid;
}

void indirectModel(Mixer& m) {
  static ContextMap cm(MEM(), 12);
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

struct DMCNode {
  unsigned int nx[2];
  U8 state;
  unsigned int c0:12, c1:12;
};

void dmcModel(Mixer& m) {
  static int top=0, curr=0;
  static Array<DMCNode> t(MEM()*2);
  static StateMap sm;
  static int threshold=256;

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
   } else if (t[curr].c0<=3840) t[curr].c0+=256;
  t[curr].state=nex(t[curr].state, y);
  curr=t[curr].nx[y];

  const int pr1=sm.p(t[curr].state);
  const int n1=t[curr].c1;
  const int n0=t[curr].c0;
  const int pr2=(n1+5)*4096/(n0+n1+10);
  m.add(stretch(pr1)/4);
  m.add(stretch(pr2)/4);
}

/*
  XML Model.
  Attempts to parse the tag structure and detect specific content types.

  Changelog:
  (17/08/2017) v96: Initial release by Márcio Pais
  (17/08/2017) v97: Bug fixes (thank you Mauro Vezzosi) and improvements.
*/

struct XMLAttribute {
  U32 Name, Value, Length;
};

struct XMLContent {
  U32 Data, Length, Type;
};

struct XMLTag {
  U32 Name, Length;
  int Level;
  bool EndTag, Empty;
  XMLContent Content;
  struct XMLAttributes {
    XMLAttribute Items[4];
    U32 Index;
  } Attributes;
};

struct XMLTagCache {
  XMLTag Tags[CacheSize];
  U32 Index;
};

enum ContentFlags {
  Text        = 0x001,
  Number      = 0x002,
  Date        = 0x004,
  Time        = 0x008,
  URL         = 0x010,
  Link        = 0x020,
  Coordinates = 0x040,
  Temperature = 0x080,
  ISBN        = 0x100,
};

enum XMLState {
  None               = 0,
  ReadTagName        = 1,
  ReadTag            = 2,
  ReadAttributeName  = 3,
  ReadAttributeValue = 4,
  ReadContent        = 5,
  ReadCDATA          = 6,
  ReadComment        = 7,
};

#define DetectContent(){ \
  if ((c4&0xF0F0F0F0)==0x30303030){ \
    int i = 0, j = 0; \
    while ((i<4) && ( (j=(c4>>(8*i))&0xFF)>=0x30 && j<=0x39 )) \
      i++; \
\
    if (i==4 && ( ((c8&0xFDF0F0FD)==0x2D30302D && buf(9)>=0x30 && buf(9)<=0x39) || ((c8&0xF0FDF0FD)==0x302D302D) )) \
      (*Content).Type |= Date; \
  } \
  else if (((c8&0xF0F0FDF0)==0x30302D30 || (c8&0xF0F0F0FD)==0x3030302D) && buf(9)>=0x30 && buf(9)<=0x39){ \
    int i = 2, j = 0; \
    while ((i<4) && ( (j=(c8>>(8*i))&0xFF)>=0x30 && j<=0x39 )) \
      i++; \
\
    if (i==4 && (c4&0xF0FDF0F0)==0x302D3030) \
      (*Content).Type |= Date; \
  } \
\
  if ((c4&0xF0FFF0F0)==0x303A3030 && buf(5)>=0x30 && buf(5)<=0x39 && ((buf(6)<0x30 || buf(6)>0x39) || ((c8&0xF0F0FF00)==0x30303A00 && (buf(9)<0x30 || buf(9)>0x39)))) \
    (*Content).Type |= Time; \
\
  if ((*Content).Length>=8 && (c8&0x80808080)==0 && (c4&0x80808080)==0) \
    (*Content).Type |= Text; \
\
  if ((c8&0xF0F0FF)==0x3030C2 && (c4&0xFFF0F0FF)==0xB0303027){ \
    int i = 2; \
    while ((i<7) && buf(i)>=0x30 && buf(i)<=0x39) \
      i+=(i&1)*2+1; \
\
    if (i==10) \
      (*Content).Type |= Coordinates; \
  } \
\
  if ((c4&0xFFFFFA)==0xC2B042 && B!=0x47 && (((c4>>24)>=0x30 && (c4>>24)<=0x39) || ((c4>>24)==0x20 && (buf(5)>=0x30 && buf(5)<=0x39)))) \
    (*Content).Type |= Temperature; \
\
  if (B>=0x30 && B<=0x39) \
    (*Content).Type |= Number; \
\
  if (c4==0x4953424E && buf(5)==0x20) \
    (*Content).Type |= ISBN; \
}

void XMLModel(Mixer& m, ModelStats *Stats = NULL){
  static ContextMap cm(MEM()/4, 4);
  static XMLTagCache Cache;
  static U32 StateBH[8];
  static XMLState State = None, pState = None;
  static U32 c8, WhiteSpaceRun = 0, pWSRun = 0, IndentTab = 0, IndentStep = 2, LineEnding = 2;

  if (bpos==0) {
    U8 B = (U8)c4;
    XMLTag *pTag = &Cache.Tags[ (Cache.Index-1)&(CacheSize-1) ], *Tag = &Cache.Tags[ Cache.Index&(CacheSize-1) ];
    XMLAttribute *Attribute = &((*Tag).Attributes.Items[ (*Tag).Attributes.Index&3 ]);
    XMLContent *Content = &(*Tag).Content;
    pState = State;
    c8 = (c8<<8)|buf(5);
    if ((B==0x09 || B==0x20) && (B==(U8)(c4>>8) || !WhiteSpaceRun)){
      WhiteSpaceRun++;
      IndentTab = (B==0x09);
    }
    else{
      if ((State==None || (State==ReadContent && (*Content).Length<=LineEnding+WhiteSpaceRun)) && WhiteSpaceRun>1+IndentTab && WhiteSpaceRun!=pWSRun){
        IndentStep=abs((int)(WhiteSpaceRun-pWSRun));
        pWSRun = WhiteSpaceRun;
      }
      WhiteSpaceRun=0;
    }
    if (B==0x0A)
      LineEnding = 1+((U8)(c4>>8)==0x0D);

    switch (State){
      case None : {
        if (B==0x3C){
          State = ReadTagName;
          memset(Tag, 0, sizeof(XMLTag));
          (*Tag).Level = ((*pTag).EndTag || (*pTag).Empty)?(*pTag).Level:(*pTag).Level+1;
        }
        if ((*Tag).Level>1)
          DetectContent();

        cm.set(hash(pState, State, ((*pTag).Level+1)*IndentStep - WhiteSpaceRun));
        break;
      }
      case ReadTagName : {
        if ((*Tag).Length>0 && (B==0x09 || B==0x0A || B==0x0D || B==0x20))
          State = ReadTag;
        else if ((B==0x3A || (B>='A' && B<='Z') || B==0x5F || (B>='a' && B<='z')) || ((*Tag).Length>0 && (B==0x2D || B==0x2E || (B>='0' && B<='9')))){
          (*Tag).Length++;
          (*Tag).Name = (*Tag).Name * 263 * 32 + (B&0xDF);
        }
        else if (B == 0x3E){
          if ((*Tag).EndTag){
            State = None;
            Cache.Index++;
          }
          else
            State = ReadContent;
        }
        else if (B!=0x21 && B!=0x2D && B!=0x2F && B!=0x5B){
          State = None;
          Cache.Index++;
        }
        else if ((*Tag).Length==0){
          if (B==0x2F){
            (*Tag).EndTag = true;
            (*Tag).Level = max(0,(*Tag).Level-1);
          }
          else if (c4==0x3C212D2D){
            State = ReadComment;
            (*Tag).Level = max(0,(*Tag).Level-1);
          }
        }

        if ((*Tag).Length==1 && (c4&0xFFFF00)==0x3C2100){
          memset(Tag, 0, sizeof(XMLTag));
          State = None;
        }
        else if ((*Tag).Length==5 && c8==0x215B4344 && c4==0x4154415B){
          State = ReadCDATA;
          (*Tag).Level = max(0,(*Tag).Level-1);
        }

        int i = 1;
        do{
          pTag = &Cache.Tags[ (Cache.Index-i)&(CacheSize-1) ];
          i+=1+((*pTag).EndTag && Cache.Tags[ (Cache.Index-i-1)&(CacheSize-1) ].Name==(*pTag).Name);
        }
        while ( i<CacheSize && ((*pTag).EndTag || (*pTag).Empty) );

        cm.set(hash(pState*8+State, (*Tag).Name, (*Tag).Level, (*pTag).Name, (*pTag).Level!=(*Tag).Level ));
        break;
      }
      case ReadTag : {
        if (B==0x2F)
          (*Tag).Empty = true;
        else if (B==0x3E){
          if ((*Tag).Empty){
            State = None;
            Cache.Index++;
          }
          else
            State = ReadContent;
        }
        else if (B!=0x09 && B!=0x0A && B!=0x0D && B!=0x20){
          State = ReadAttributeName;
          (*Attribute).Name = B&0xDF;
        }
        cm.set(hash(pState, State, (*Tag).Name, B, (*Tag).Attributes.Index ));
        break;
      }
      case ReadAttributeName : {
        if ((c4&0xFFF0)==0x3D20 && (B==0x22 || B==0x27)){
          State = ReadAttributeValue;
          if ((c8&0xDFDF)==0x4852 && (c4&0xDFDF0000)==0x45460000)
            (*Content).Type |= Link;
        }
        else if (B!=0x22 && B!=0x27 && B!=0x3D)
          (*Attribute).Name = (*Attribute).Name * 263 * 32 + (B&0xDF);

        cm.set(hash(pState*8+State, (*Attribute).Name, (*Tag).Attributes.Index, (*Tag).Name, (*Content).Type ));
        break;
      }
      case ReadAttributeValue : {
        if (B==0x22 || B==0x27){
          (*Tag).Attributes.Index++;
          State = ReadTag;
        }
        else{
          (*Attribute).Value = (*Attribute).Value* 263 * 32 + (B&0xDF);
          (*Attribute).Length++;
          if ((c8&0xDFDFDFDF)==0x48545450 && ((c4>>8)==0x3A2F2F || c4==0x733A2F2F))
            (*Content).Type |= URL;
        }
        cm.set(hash(pState, State, (*Attribute).Name, (*Content).Type ));
        break;
      }
      case ReadContent : {
        if (B==0x3C){
          State = ReadTagName;
          Cache.Index++;
          memset(&Cache.Tags[ Cache.Index&(CacheSize-1) ], 0, sizeof(XMLTag));
          Cache.Tags[ Cache.Index&(CacheSize-1) ].Level = (*Tag).Level+1;
        }
        else{
          (*Content).Length++;
          (*Content).Data = (*Content).Data * 997*16 + (B&0xDF);

          DetectContent();
        }
        cm.set(hash(pState, State, (*Tag).Name, c4&0xC0FF ));
        break;
      }
      case ReadCDATA : {
        if ((c4&0xFFFFFF)==0x5D5D3E){
          State = None;
          Cache.Index++;
        }
        cm.set(hash(pState, State));
        break;
      }
      case ReadComment : {
        if ((c4&0xFFFFFF)==0x2D2D3E){
          State = None;
          Cache.Index++;
        }
        cm.set(hash(pState, State));
        break;
      }
    }

    StateBH[State] = (StateBH[State]<<8)|B;
    pTag = &Cache.Tags[ (Cache.Index-1)&(CacheSize-1) ];
    cm.set(hash(State, (*Tag).Level, pState*2+(*Tag).EndTag, (*Tag).Name));
    cm.set(hash((*pTag).Name, State*2+(*pTag).EndTag, (*pTag).Content.Type, (*Tag).Content.Type));
    cm.set(hash(State*2+(*Tag).EndTag, (*Tag).Name, (*Tag).Content.Type, c4&0xE0FF));
  }
  cm.mix(m);
  U8 s = ((StateBH[State]>>(28-bpos))&0x08) |
         ((StateBH[State]>>(21-bpos))&0x04) |
         ((StateBH[State]>>(14-bpos))&0x02) |
         ((StateBH[State]>>( 7-bpos))&0x01) |
         ((bpos)<<4);
  if (Stats)
    (*Stats).XML = (s<<3)|State;
}

U32 last_prediction = 2048;

int contextModel2() {
  static ContextMap2 cm(MEM()*31, 9);
  static TextModel textModel(MEM()*16);
  static RunContextMap rcm7(MEM()), rcm9(MEM()), rcm10(MEM());
  static Mixer m(NUM_INPUTS, 10800+1024*21+16384, NUM_SETS, 32);
  static U32 cxt[16];
  static Filetype ft2,filetype=preprocessor::DEFAULT;
  static int size=0;  // bytes remaining in block
  static int info=0;  // image width or audio type
  static ModelStats stats;

  // Parse filetype and size
  if (bpos==0) {
    --size;
    ++blpos;
    if (size==-1) info=0, ft2=(Filetype)buf(1);
		if (size==-5 && !preprocessor::HasInfo(ft2)) {
			size=buf(4)<<24|buf(3)<<16|buf(2)<<8|buf(1);
      blpos=0;
    }
    if (size==-9) {
			size=buf(8)<<24|buf(7)<<16|buf(6)<<8|buf(5);
      info=buf(4)<<24|buf(3)<<16|buf(2)<<8|buf(1);
      blpos=0;
      if (ft2==preprocessor::TEXT && info) size = info-8;
    }
    if (!blpos) filetype=ft2;
    if (size==0) filetype=preprocessor::DEFAULT;
  }

  m.update();
  m.add(64);

  int ismatch=ilog(matchModel(m));
  if (filetype==preprocessor::IMAGE1) return im1bitModel(m, info), m.p();
  if (filetype==preprocessor::IMAGE4) return im4bitModel(m, info), m.p();
  if (filetype==preprocessor::IMAGE8) return im8bitModel(m, info), m.p();
  if (filetype==preprocessor::IMAGE8GRAY) return im8bitModel(m, info, 1), m.p();
  if (filetype==preprocessor::IMAGE24) return im24bitModel(m, info), m.p();
  if (filetype==preprocessor::IMAGE32) return im24bitModel(m, info, 1), m.p();
  if ((filetype!=preprocessor::EXE && jpegModel(m)) || (size>0 && imgModel(m, &stats)) || audioModel(m, &stats))
    return m.p();

  if (bpos==0) {
    for (int i=15; i>0; --i)
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
    recordModel(m, filetype, &stats);
    recordModel1(m);
    wordModel(m);
    nestModel(m);
    indirectModel(m);
    dmcModel(m);
    XMLModel(m, &stats);
    textModel.Predict(m, buf, &stats);
    exeModel(m, true, &stats);
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

class Predictor {
  int pr;
public:
  Predictor();
  int p() const {return pr;}
  void update();
};

Predictor::Predictor(): pr(2048) {
  for (int i=0; i<1024; ++i)
    dt[i]=16384/(i+i+3);
}

void Predictor::update() {
  static APM a(256), a1(0x10000), a2(0x10000), a3(0x10000),
                      a4(0x10000), a5(0x10000), a6(0x10000);
  static U32 x5=0;

  c0+=c0+y;
  if (c0>=256) {
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

  int pr0=contextModel2();
  AddPrediction(pr0);

  pr=a.p(pr0, c0);
  AddPrediction(pr);

  int pr1=a1.p(pr0, c0+256*buf(1));
  AddPrediction(pr1);
  int pr2=a2.p(pr0, c0^(hash(buf(1), buf(2))&0xffff));
  AddPrediction(pr2);
  int pr3=a3.p(pr0, c0^(hash(buf(1), buf(2), buf(3))&0xffff));
  AddPrediction(pr3);
  pr0=(pr0+pr1+pr2+pr3+2)>>2;
  AddPrediction(pr0);

  pr1=a4.p(pr, c0+256*buf(1));
  AddPrediction(pr1);
  pr2=a5.p(pr, c0^(hash(buf(1), buf(2))&0xffff));
  AddPrediction(pr2);
  pr3=a6.p(pr, c0^(hash(buf(1), buf(2), buf(3))&0xffff));
  AddPrediction(pr3);
  pr=(pr+pr1+pr2+pr3+2)>>2;
  AddPrediction(pr);

  pr=(pr+pr0+1)>>1;
  AddPrediction(pr);
  ResetPredictions();
  last_prediction = pr;
}

Predictor paq8;
}

PAQ8::PAQ8(int memory) {
  level = memory;
  buf.setsize(MEM()*8);
}

const std::valarray<float>& PAQ8::Predict() {
  return model_predictions;
}

unsigned int PAQ8::NumOutputs() {
  return model_predictions.size();
}

void PAQ8::Perceive(int bit) {
  y = bit;
  paq8.update();
}
