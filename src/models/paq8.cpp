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
#include <unordered_map>
#include <memory>
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

namespace paq8 {
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
U8 grp0; // Quantized partial byte as ASCII group

int dt[1024];  // i -> 16K/(i+3)

struct ModelStats{
  Filetype Type;
  U64 Misses;
  struct {
    U32 length;      //used by SSE stage
    U8 expectedByte; //used by SSE stage
  } Match;
  struct {
    struct {
      U8 WW, W, NN, N, Wp1, Np1; 
    } pixels;        //used by SSE stage
    U8 plane;        //used by SSE stage
    U8 ctx;          //used by SSE stage
  } Image;
  U32 XML, x86_64, Record;
  struct {
    U8 state:3;      //unused
    U8 lastPunct:5;  //unused
    U8 wordLength:4; //unused
    U8 boolmask:4;   //unused
    U8 firstLetter;  //used by SSE stage
    U8 mask;         //used by SSE stage
  } Text;
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

#define NUM_INPUTS 1552
#define NUM_SETS 28

std::valarray<float> model_predictions(0.5, NUM_INPUTS + NUM_SETS + 11);
unsigned int prediction_index = 0;
float conversion_factor = 1.0 / 4095;

void AddPrediction(int x) {
  //printf("%d\n", prediction_index);
  model_predictions[prediction_index++] = x * conversion_factor;
}

void ResetPredictions() {
  prediction_index = 0;
}

class Mixer {
  const int N, S, init_w;
  Array<short, 16> tx;
  std::unordered_map<unsigned int, std::unique_ptr<Array<short, 16>>> wx_;

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
      int context = cxt[i];
      auto* wts = wx_[context].get();
      if (wts == nullptr) {
        wx_[context] = std::unique_ptr<Array<short, 16>>(new Array<short, 16>(N));
        wts = wx_[context].get();
        for (int i=0; i<N; ++i) (*wts)[i]=init_w;
      }
      train(&tx[0], &(*wts)[0], nx, err);
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
        int context = cxt[i];
        auto* wts = wx_[context].get();
        if (wts == nullptr) {
          wx_[context] = std::unique_ptr<Array<short, 16>>(new Array<short, 16>(N));
          wts = wx_[context].get();
          for (int i=0; i<N; ++i) (*wts)[i]=init_w;
        }
        pr[i]=squash((dot_product(&tx[0], &(*wts)[0], nx) * 9)>>9);
        mp->add(stretch(pr[i]));
      }
      mp->set(0, 1);
      return mp->p();
    }
    else {
      auto* wts = wx_[0].get();
      if (wts == nullptr) {
        wx_[0] = std::unique_ptr<Array<short, 16>>(new Array<short, 16>(N));
        wts = wx_[0].get();
        for (int i=0; i<N; ++i) (*wts)[i]=init_w;          
      }

      int z = dot_product(&tx[0], &(*wts)[0], nx);
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
    N((n+7)&-8), S(s), init_w(w),
    tx(N), cxt(S), ncxt(0), base(0), nx(0), pr(S), mp(0) {
  for (int i=0; i<S; ++i)
    pr[i]=2048;
  if (S>1) {
    mp=new Mixer(S, 1, 1, 0x7fff);
  }
}

class APM1 {
  int index;
  const int N;
  Array<U16> t;
public:
  APM1(int n);
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

APM1::APM1(int n): index(0), N(n), t(n*33) {
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
  StateMap32(const int n=256, const bool init=true);
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

StateMap32::StateMap32(const int n, const bool init): N(n), cxt(0), t(n) {
  if (init && (N==256)) {
    for (int i=0; i<N; ++i) {
      U32 n0=nex(i,2);
      U32 n1=nex(i,3);
      if (n0==0) n1*=64;
      if (n1==0) n0*=64;
      t[i] = ((n1<<16)/(n0+n1+1))<<16;
    }
  }
  else {
    for (int i=0; i<N; ++i)
      t[i]=(1u<<31)+0;  //initial p=0.5, initial count=0
  }
}

class APM : public StateMap32 {
public:
  APM(int n) : StateMap32(n*24) {
    for (int i=0; i<N; ++i) {
      int p = ((i%24*2+1)*4096)/48-2048;
      t[i] = (U32(squash(p))<<20)+6; //initial count: 6
    }
  }
  int p(int pr, int cx, const int limit=0xFF) {
    //adapt (update prediction from previous pass)
    update(limit);
    //predict
    pr = (stretch(pr)+2048)*23;
    int wt = pr&0xfff;  // interpolation weight (0..4095)
    cx = cx*24+(pr>>12);
    cxt = cx+(wt>>11);
    pr = ((t[cx]>>13)*(4096-wt)+(t[cx+1]>>13)*wt)>>19;
    return pr;
  }
};

/////////////////////////////////////////////////////////////////////////////////////////

#define PHI32 UINT32_C(0x9E3779B9) // 2654435769
#define PHI64 UINT64_C(0x9E3779B97F4A7C15) // 11400714819323198485
#define MUL64_1  UINT64_C(0x993DDEFFB1462949)
#define MUL64_2  UINT64_C(0xE9C91DC159AB0D2D)
#define MUL64_3  UINT64_C(0x83D6A14F1B0CED73)
#define MUL64_4  UINT64_C(0xA14F1B0CED5A841F)
#define MUL64_5  UINT64_C(0xC0E51314A614F4EF) 
#define MUL64_6  UINT64_C(0xDA9CC2600AE45A27)
#define MUL64_7  UINT64_C(0x826797AA04A65737)
#define MUL64_8  UINT64_C(0x2375BE54C41A08ED)
#define MUL64_9  UINT64_C(0xD39104E950564B37)
#define MUL64_10 UINT64_C(0x3091697D5E685623)
#define MUL64_11 UINT64_C(0x20EB84EE04A3C7E1)
#define MUL64_12 UINT64_C(0xF501F1D0944B2383)
#define MUL64_13 UINT64_C(0xE3E4E8AA829AB9B5)

static inline U32 finalize64(const U64 hash, const int hashbits) {
  return U32(hash>>(64-hashbits));
}

static inline U64 checksum64(const U64 hash, const int hashbits, const int checksumbits) {
  return hash>>(64-hashbits-checksumbits); 
}

static inline U64 hash(const U64 x0) {
  return (x0+1)*PHI64;
}
static inline U64 hash(const U64 x0, const U64 x1) {
  return (x0+1)*PHI64   + (x1+1)*MUL64_1;
}
static inline U64 hash(const U64 x0, const U64 x1, const U64 x2) {
  return (x0+1)*PHI64   + (x1+1)*MUL64_1 + (x2+1)*MUL64_2;
}
static inline U64 hash(const U64 x0, const U64 x1, const U64 x2, const U64 x3) {
  return (x0+1)*PHI64   + (x1+1)*MUL64_1 + (x2+1)*MUL64_2 +
         (x3+1)*MUL64_3;
}
static inline U64 hash(const U64 x0, const U64 x1, const U64 x2, const U64 x3, const U64 x4) {
  return (x0+1)*PHI64   + (x1+1)*MUL64_1 + (x2+1)*MUL64_2 +
         (x3+1)*MUL64_3 + (x4+1)*MUL64_4;
}
static inline U64 hash(const U64 x0, const U64 x1, const U64 x2, const U64 x3, const U64 x4,
                  const U64 x5) {
  return (x0+1)*PHI64   + (x1+1)*MUL64_1 + (x2+1)*MUL64_2 +
         (x3+1)*MUL64_3 + (x4+1)*MUL64_4 + (x5+1)*MUL64_5;
}
static inline U64 hash(const U64 x0, const U64 x1, const U64 x2, const U64 x3, const U64 x4,
                  const U64 x5, const U64 x6) {
  return (x0+1)*PHI64   + (x1+1)*MUL64_1 + (x2+1)*MUL64_2 +
         (x3+1)*MUL64_3 + (x4+1)*MUL64_4 + (x5+1)*MUL64_5 +
         (x6+1)*MUL64_6;
}
static inline U64 hash(const U64 x0, const U64 x1, const U64 x2, const U64 x3, const U64 x4,
                  const U64 x5, const U64 x6, const U64 x7) {
  return (x0+1)*PHI64   + (x1+1)*MUL64_1 + (x2+1)*MUL64_2 +
         (x3+1)*MUL64_3 + (x4+1)*MUL64_4 + (x5+1)*MUL64_5 +
         (x6+1)*MUL64_6 + (x7+1)*MUL64_7;
}
static inline U64 combine64(const U64 seed, const U64 x) {
  return hash(seed+x);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <int B> class BH {
  enum {M=8};  // search limit
  Array<U8> t; // elements
  const U32 mask; // size-1
  const int hashbits;
public:
  BH(int i): t(i*B), mask(i-1), hashbits(ilog2(mask+1)) {}
  U8* operator[](const U64 i);
};

template <int B>
inline U8* BH<B>::operator[](const U64 ctx) {
  const U16 chk=checksum64(ctx,hashbits,16)&0xffff;
  const U32 i=finalize64(ctx,hashbits)*M & mask;
  U8 *p;
  U16 *cp;
  int j;
  for (j=0; j<M; ++j) {
    p=&t[(i+j)*B];
    cp=(U16*)p;
    if (p[2]==0) {*cp=chk;break;} // empty slot
    if (*cp==chk) break;  // found
  }
  if (j==0) return p+1;  // front
  static U8 tmp[B];  // element to move to front
  if (j==M) {
    --j;
    memset(tmp, 0, B);
    memmove(tmp, &chk, 2);
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
private:
  Array<U8,64> t;  // storage area for n items (1 item = B bytes): 0:checksum 1:priority 2:data 3:data  ... B-1:data
  const int N;     // number of items in table
  const int mask;
  const int hashbits;
public:
  HashTable(int n): t(n), N(n), mask(N-1), hashbits(ilog2(mask+1)) {}
  U8* operator[](U64 i);
};

template <int B>
inline U8* HashTable<B>::operator[](U64 i) { //i: context selector
  U8 chk=checksum64(i,hashbits,8)&0xff; // 8-bit checksum
  i=finalize64(i,hashbits)*B & mask;    // index: force bounds
  //search for the checksum in t
  U8 *p = &t[0];
  if (p[i]==chk) return p+i+1;
  if (p[i^B]==chk) return p+(i^B)+1;
  if (p[i^(B*2)]==chk) return p+(i^(B*2))+1;
  //not found, let's overwrite the lowest priority element
  if (p[i+1]>p[(i+1)^B] || p[i+1]>p[(i+1)^(B*2)]) i^=B;
  if (p[i+1]>p[(i+1)^B^(B*2)]) i^=B^(B*2);
  memset(p+i, 0, B);
  p[i]=chk;
  return p+i+1;
}

class RunContextMap {
  BH<4> t;
  U8* cp;
public:
  RunContextMap(int m): t(m/4) {cp=t[0]+1;}
  void set(U64 cx) {
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

/*
  Map for modelling contexts of (nearly-)stationary data.
  The context is looked up directly. For each bit modelled, a 16bit prediction is stored.
  The adaptation rate is controlled by the caller, see mix().

  - BitsOfContext: How many bits to use for each context. Higher bits are discarded.
  - BitsPerContext: How many bits [1..8] of input are to be modelled for each context.
  New contexts must be set at those intervals.

  Uses (2^(BitsOfContext+1))*((2^BitsPerContext)-1) bytes of memory.
*/

class SmallStationaryContextMap {
  Array<U16> Data;
  int Context, Mask, Stride, bCount, bTotal, B;
  U16 *cp;
public:
  SmallStationaryContextMap(int BitsOfContext, int BitsPerContext = 8) : Data((1ull<<BitsOfContext)*((1ull<<BitsPerContext)-1)), Context(0), Mask((1<<BitsOfContext)-1), Stride((1<<BitsPerContext)-1), bCount(0), bTotal(BitsPerContext), B(0) {
    Reset();
    cp=&Data[0];
  }
  void set(U32 ctx) {
    Context = (ctx&Mask)*Stride;
    bCount=B=0;
  }
  void Reset() {
    for (U32 i=0; i<Data.size(); ++i)
      Data[i]=0x7FFF;
  }
  void mix(Mixer& m, const int rate = 7, const int Multiplier = 1, const int Divisor = 4) {
    *cp+=((y<<16)-(*cp)+(1<<(rate-1)))>>rate;
    B+=(y && B>0);
    cp = &Data[Context+B];
    int Prediction = (*cp)>>4;
    m.add((stretch(Prediction)*Multiplier)/Divisor);
    m.add(((Prediction-2048)*Multiplier)/(Divisor*2));
    bCount++; B+=B+1;
    if (bCount==bTotal)
      bCount=B=0;
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

  Uses (2^(BitsOfContext+2))*((2^BitsPerContext)-1) bytes of memory.
*/

class StationaryMap {
  Array<U32> Data;
  const int mask, maskbits, stride;
  int Context, bCount, bTotal, B;
  U32 *cp;
public:
  StationaryMap(int BitsOfContext, int BitsPerContext = 8, int Rate = 0): Data((1ull<<BitsOfContext)*((1ull<<BitsPerContext)-1)), mask((1<<BitsOfContext)-1), maskbits(BitsOfContext), stride((1<<BitsPerContext)-1), Context(0), bCount(0), bTotal(BitsPerContext), B(0) {
    Reset(Rate);
    cp=&Data[0];
  }
  void set_direct(U32 ctx) {
    Context = (ctx&mask)*stride;
    bCount=B=0;
  }
  void set(U64 ctx) {
    Context = (finalize64(ctx,maskbits)&mask)*stride;
    bCount=B=0;
  }
  void Reset( int Rate = 0 ){
    for (U32 i=0; i<Data.size(); ++i)
      Data[i]=(0x7FF<<20)|min(1023,Rate);
  }
  void mix(Mixer& m, const int Multiplier = 1, const int Divisor = 4, const U16 Limit = 1023) {
    // update
    U32 Count = min(min(Limit,0x3FF), ((*cp)&0x3FF)+1);
    int Prediction = (*cp)>>10, Error = (y<<22)-Prediction;
    Error = ((Error/8)*dt[Count])/1024;
    Prediction = min(0x3FFFFF,max(0,Prediction+Error));
    *cp = (Prediction<<10)|Count;
    // predict
    B+=(y && B>0);
    cp=&Data[Context+B];
    Prediction = (*cp)>>20;
    m.add((stretch(Prediction)*Multiplier)/Divisor);
    m.add(((Prediction-2048)*Multiplier)/(Divisor*2));
    bCount++; B+=B+1;
    if (bCount==bTotal)
      bCount=B=0;
  }
};

class IndirectMap {
  Array<U8> Data;
  StateMap32 Map;
  const int mask, maskbits, stride;
  int Context, bCount, bTotal, B;
  U8 *cp;
public:
  IndirectMap(int BitsOfContext, int BitsPerContext = 8): Data((1ull<<BitsOfContext)*((1ull<<BitsPerContext)-1)), mask((1<<BitsOfContext)-1), maskbits(BitsOfContext), stride((1<<BitsPerContext)-1), Context(0), bCount(0), bTotal(BitsPerContext), B(0) {
    cp=&Data[0];
  }
  void set_direct(const U32 ctx) {
    Context = (ctx&mask)*stride;
    bCount=B=0;
  }
  void set(const U64 ctx) {
    Context = (finalize64(ctx,maskbits)&mask)*stride;
    bCount=B=0;
  }
  void mix(Mixer& m, const int Multiplier = 1, const int Divisor = 4, const U16 Limit = 1023) {
    // update
    *cp = nex(*cp, y);
    // predict
    B+=(y && B>0);
    cp=&Data[Context+B];
    const U8 state = *cp;
    const int p1 = Map.p(state, Limit);
    m.add((stretch(p1)*Multiplier)/Divisor);
    m.add(((p1-2048)*Multiplier)/(Divisor*2));
    bCount++; B+=B+1;
    if (bCount==bTotal)
      bCount=B=0;
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
  Array<U16> chk;
  Array<U8*> runp;
  StateMap *sm;
  int cn;
  const U32 mask;
  const int hashbits;
  void update(U32 cx, int c);
  int mix1(Mixer& m, int cc, int bp, int c1, int y1);
public:
  ContextMap(U64 m, int c);
  ~ContextMap();
  void set(U64 cx);
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

ContextMap::ContextMap(U64 m, int c): C(c), t(m>>6), cp(c), cp0(c),
    cxt(c), chk(c), runp(c), cn(0), mask(U32(t.size()-1)), hashbits(ilog2(mask+1)) {
  sm=new StateMap[C];
  for (int i=0; i<C; ++i) {
    cp0[i]=cp[i]=&t[0].bh[0][0];
    runp[i]=cp[i]+3;
  }
}

ContextMap::~ContextMap() {
  delete[] sm;
}

inline void ContextMap::set(U64 cx) {
  cx=hash(cx, cn);
  cxt[cn]=finalize64(cx,hashbits);
  chk[cn]=checksum64(cx,hashbits,16)&0xffff;
  cn++;
}

int ContextMap::mix1(Mixer& m, int cc, int bp, int c1, int y1) {

  int result=0;
  for (int i=0; i<cn; ++i) {
    if (cp[i]) {
      int ns=nex(*cp[i], y1);
      if (ns>=204 && rnd() << ((452-ns)>>3)) ns-=4;
      *cp[i]=ns;
    }

    if (bp>1 && runp[i][0]==0)
      cp[i]=nullptr;
    else {
      switch(bp) {
        case 1: case 3: case 6: cp[i]=cp0[i]+1+(c0&1); break;
        case 4: case 7: cp[i]=cp0[i]+3+(c0&3); break;
        case 2: case 5: {
          const U16 checksum = chk[i];
          const U32 ctx = cxt[i];
          cp0[i]=cp[i]=t[(ctx+c0)&mask].get(checksum); break;
        }
        default:
        {
          const U16 checksum = chk[i];
          const U32 ctx = cxt[i];
          cp0[i]=cp[i]=t[(ctx+c0)&mask].get(checksum);
          // Update pending bit histories for bits 2-7
          if (cp0[i][3]==2) {
            const int c=cp0[i][4]+256;
            U8 *p=t[(ctx+(c>>6))&mask].get(checksum);
            p[0]=1+((c>>5)&1);
            p[1+((c>>5)&1)]=1+((c>>4)&1);
            p[3+((c>>4)&3)]=1+((c>>3)&1);
            p=t[(ctx+(c>>3))&mask].get(checksum);
            p[0]=1+((c>>2)&1);
            p[1+((c>>2)&1)]=1+((c>>1)&1);
            p[3+((c>>1)&3)]=1+(c&1);
            cp0[i][6]=0;
          }
          // Update run count of previous context
          if (runp[i][0]==0)  // new context
            runp[i][0]=2, runp[i][1]=c1;
          else if (runp[i][1]!=c1)  // different byte in context
            runp[i][0]=1, runp[i][1]=c1;
          else if (runp[i][0]<254)  // same byte in context
            runp[i][0]+=2;
          else if (runp[i][0]==255)
            runp[i][0]=128;
          runp[i]=cp0[i]+3;
        } break;
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

    const int s = cp[i]!=nullptr ?  *cp[i] : 0;
    int p1=sm[i].p(s);
    const int st=(stretch(p1)+(1<<1))>>2;
    m.add(st);
    m.add((p1-2047+(1<<2))>>3);
    const int n0=-!nex(s,2);
    const int n1=-!nex(s,3);
    m.add(st*abs(n1-n0));
    const int p0=4095-p1;
    m.add(((p1&n0)-(p0&n1)+(1<<3))>>4);
    result+=s>0;
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
  Array<U32> Contexts; // C whole byte context hashes
  Array<U16> Chk;      // C whole byte context checksums
  Array<bool> HasHistory; // True if context has a full valid byte history (i.e., seen at least 3 times)
  StateMap32 **Maps6b, **Maps8b, **Maps12b;
  U32 index; // Next context to set by set()
  const U32 mask;
  const int hashbits;
  U32 bits;
  U8 lastByte, lastBit, bitPos;
  inline void Update() {
    for (U32 i=0; i<index; i++) {
      if (BitState[i])
        *BitState[i] = nex(*BitState[i], lastBit);

      if (bitPos>1 && ByteHistory[i][0]==0)
        BitState[i] = nullptr;
      else {
        switch (bitPos) {
          case 0: {
            const U16 chk = Chk[i];
            const U32 ctx = Contexts[i];
            BitState[i] = BitState0[i] = Table[(ctx+bits)&mask].Find(chk);
            // Update pending bit histories for bits 2-7
            if (BitState0[i][3]==2) {
              const int c = BitState0[i][4]+256;
              U8 *p = Table[(ctx+(c>>6))&mask].Find(chk);
              p[0] = 1+((c>>5)&1);
              p[1+((c>>5)&1)] = 1+((c>>4)&1);
              p[3+((c>>4)&3)] = 1+((c>>3)&1);
              p = Table[(ctx+(c>>3))&mask].Find(chk);
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
            const U16 chk = Chk[i];
            const U32 ctx = Contexts[i];
            BitState[i] = BitState0[i] = Table[(ctx+bits)&mask].Find(chk);
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
  ContextMap2(const U64 Size, const U32 Count) : C(Count), Table(Size>>6), BitState(Count), BitState0(Count), ByteHistory(Count), Contexts(Count), Chk(Count), HasHistory(Count), mask(U32(Table.size()-1)), hashbits(ilog2(mask+1)) {
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
  inline void set(U64 ctx) { // set next whole byte context to ctx
    ctx=hash(ctx, index);
    Contexts[index]=finalize64(ctx,hashbits);
    Chk[index]=checksum64(ctx,hashbits,16)&0xffff;
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
      int state = (BitState[i]!=nullptr)?*BitState[i]:0;
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
      m.add((st*abs(n1-n0)));
      m.add((p1&n0)-(p0&n1));
      m.add(stretch(Maps12b[i]->p((state<<9)|(bitPos<<6)|k))>>2);
      m.add(stretch(Maps6b[i]->p((state<<3)|bitPos))>>2);
    }
    if (bitPos==7) index = 0;
    return result;
  }
};

///////////////// Ordinary Least Squares predictor /////////////////

template <typename F, typename T, const bool hasZeroMean = true>
class OLS {
  static constexpr F ftol = 1E-8;
  static constexpr F sub = F(int64_t(!hasZeroMean)<<(8*sizeof(T)-1));
private:
  int n, kmax, km, index;
  F lambda, nu;
  F *x, *w, *b;
  F **mCovariance, **mCholesky;
  int Factor() {
    // copy the matrix
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        mCholesky[i][j] = mCovariance[i][j];

    for (int i=0; i<n; i++)
      mCholesky[i][i] += nu;
    for (int i=0; i<n; i++) {
      for (int j=0; j<i; j++) {
        F sum = mCholesky[i][j];
        for (int k=0; k<j; k++)
          sum -= (mCholesky[i][k] * mCholesky[j][k]);
        mCholesky[i][j] = sum / mCholesky[j][j];
      }
      F sum = mCholesky[i][i];
      for (int k=0; k<i; k++)
        sum -= (mCholesky[i][k] * mCholesky[i][k]);
      if (sum>ftol)
        mCholesky[i][i] = sqrt(sum);
      else
        return 1;
    }
    return 0;
  }

  void Solve() {
    for (int i=0; i<n; i++) {
      F sum = b[i];
      for (int j=0; j<i; j++)
        sum -= (mCholesky[i][j] * w[j]);
      w[i] = sum / mCholesky[i][i];
    }
    for (int i=n-1; i>=0; i--) {
      F sum = w[i];
      for (int j=i+1; j<n; j++)
        sum -= (mCholesky[j][i] * w[j]);
      w[i] = sum / mCholesky[i][i];
    }
  }
public:
  OLS(int n, int kmax=1, F lambda=0.998, F nu=0.001) : n(n), kmax(kmax), lambda(lambda), nu(nu) {
    km = index = 0;
    x = new F[n], w = new F[n], b = new F[n];
    mCovariance = new F*[n], mCholesky = new F*[n];
    for (int i=0; i<n; i++) {
      x[i] = w[i] = b[i] = 0.;
      mCovariance[i] = new F[n], mCholesky[i] = new F[n];
      for (int j=0; j<n; j++)
        mCovariance[i][j] = mCholesky[i][j] = 0.;
    }
  }
  ~OLS() {
    delete[] x, delete[] w, delete[] b;
    for (int i=0; i<n; i++) {
      delete[] mCovariance[i];
      delete[] mCholesky[i];
    }
    delete[] mCovariance, delete[] mCholesky;
  }
  void Add(const T val) {
    if (index<n)
      x[index++] = F(val)-sub;
  }
  void AddFloat(const F val) {
    if (index<n)
      x[index++] = val-sub;
  }
  F Predict(const T **p) {
    F sum = 0.;
    for (int i=0; i<n; i++)
      sum += w[i] * (x[i] = F(*p[i])-sub);
    return sum+sub;
  }
  F Predict() {
    index = 0;
    F sum = 0.;
    for (int i=0; i<n; i++)
      sum += w[i] * x[i];
    return sum+sub;
  }
  void Update(const T val) {
    for (int j=0; j<n; j++)
      for (int i=0; i<n; i++)
        mCovariance[j][i] = lambda * mCovariance[j][i] + (1.0 - lambda) * (x[j] * x[i]);
    for (int i=0; i<n; i++)
      b[i] = lambda * b[i] + (1.0 - lambda) * (x[i] * (F(val)-sub));
    km++;
    if (km>=kmax) {
      if (!Factor()) Solve();
      km = 0;
    }
  }
};

////////////////////////////// Indirect Context //////////////////////////////

template <typename T>
class IndirectContext {
private:
  Array<T> data;
  T* ctx;
  U32 ctxMask, inputMask, inputBits;
public:
  IndirectContext(const int BitsPerContext, const int InputBits = 8) :
    data(1ull<<BitsPerContext),
    ctx(&data[0]),
    ctxMask((1ul<<BitsPerContext)-1),
    inputMask((1ul<<InputBits)-1),
    inputBits(InputBits)
  {
  }
  void operator+=(const U32 i) {
    (*ctx)<<=inputBits;
    (*ctx)|=i&inputMask;
  }
  void operator=(const U32 i) {
    ctx = &data[i&ctxMask];
  }
  T& operator()(void) {
    return *ctx;
  }
};

////////////////////////////// Move-to-Front list //////////////////////////////

class MTFList{
private:
  int Root, Index;
  Array<int, 16> Previous;
  Array<int, 16> Next;
public:
  MTFList(const U16 n): Root(0), Index(0), Previous(n), Next(n) {
    for (int i=0; i<n; i++) {
      Previous[i] = i-1;
      Next[i] = i+1;
    }
    Next[n-1] = -1;
  }
  inline int GetFirst(){
    return Index=Root;
  }
  inline int GetNext(){
    if (Index>=0) Index = Next[Index];
    return Index;
  }
  inline void MoveToFront(int i){
    if ((Index=i)==Root) return;
    const int p=Previous[Index], n=Next[Index];
    if (p>=0) Next[p] = Next[Index];
    if (n>=0) Previous[n] = Previous[Index];
    Previous[Root] = Index;
    Next[Index] = Root;
    Root = Index;
    Previous[Root]=-1;
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
  U64 Hash[4], Type, Language;
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
          memmove(&W->Letters[i+1], &W->Letters[i+2], W->End-i-1);
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
  (02/12/2018) v172: Sets 8 mixer contexts
*/

const U8 AsciiGroupC0[254] ={
  0, 10,
  0, 1, 10, 10,
  0, 4, 2, 3, 10, 10, 10, 10,
  0, 0, 5, 4, 2, 2, 3, 3, 10, 10, 10, 10, 10, 10, 10, 10,
  0, 0, 0, 0, 5, 5, 9, 4, 2, 2, 2, 2, 3, 3, 3, 3, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
  0, 0, 0, 0, 0, 0, 0, 0, 5, 8, 8, 5, 9, 9, 6, 5, 2, 2, 2, 2, 2, 2, 2, 8, 3, 3, 3, 3, 3, 3, 3, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 8, 8, 8, 8, 8, 5, 5, 9, 9, 9, 9, 9, 7, 8, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 8, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
};
const U8 AsciiGroup[128] = {
  0,  5,  5,  5,  5,  5,  5,  5,
  5,  5,  4,  5,  5,  4,  5,  5,
  5,  5,  5,  5,  5,  5,  5,  5,
  5,  5,  5,  5,  5,  5,  5,  5,
  6,  7,  8, 17, 17,  9, 17, 10,
  11, 12, 17, 17, 13, 14, 15, 16,
  1,  1,  1,  1,  1,  1,  1,  1,
  1,  1, 18, 19, 20, 23, 21, 22,
  23,  2,  2,  2,  2,  2,  2,  2,
  2,  2,  2,  2,  2,  2,  2,  2,
  2,  2,  2,  2,  2,  2,  2,  2,
  2,  2,  2, 24, 27, 25, 27, 26,
  27,  3,  3,  3,  3,  3,  3,  3,
  3,  3,  3,  3,  3,  3,  3,  3,
  3,  3,  3,  3,  3,  3,  3,  3,
  3,  3,  3, 28, 30, 29, 30, 30
};

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
    U64 numHashes[2]; // hashes of the last 2 numbers seen
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
    U64 asciiMask;
    U32 masks[5],
        wordLength[2];
    int UTF8Remaining;// remaining bytes for current UTF8-encoded Unicode code point (-1 if invalid byte found)
    U8 firstLetter;   // first letter of current word
    U8 firstChar;     // first character of current line
    U8 expectedDigit; // next expected digit of detected numerical sequence
    U8 prevPunct;     // most recent punctuation character seen
    Word TopicDescriptor; // last word before ':'
  } Info;
  U64 ParseCtx;
  void Update(Buf& buffer, ModelStats *Stats = nullptr);
  void SetContexts(Buf& buffer, ModelStats *Stats = nullptr);
public:
  TextModel(const U32 Size) : Map(Size, 33), Stemmers(Language::Count-1), Languages(Language::Count-1), WordPos(0x10000), State(Parse::Unknown), pState(State), Lang{ {0}, {0}, Language::Unknown, Language::Unknown }, Info{}, ParseCtx(0) {
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
    mixer.set(finalize64(hash((Lang.Id!=Language::Unknown)?1+Stemmers[Lang.Id-1]->IsVowel(buffer(1)):0, Info.masks[1]&0xFF, c0), 11), 2048);
    mixer.set(finalize64(hash(ilog2(Info.wordLength[0]+1), c0,
      (Info.lastDigit<Info.wordLength[0]+Info.wordGap)|
      ((Info.lastUpper<Info.lastLetter+Info.wordLength[1])<<1)|
      ((Info.lastPunct<Info.wordLength[0]+Info.wordGap)<<2)|
      ((Info.lastUpper<Info.wordLength[0])<<3)
    ), 11), 2048);
    mixer.set(finalize64(hash(Info.masks[1]&0x3FF, grp0, Info.lastUpper<Info.wordLength[0], Info.lastUpper<Info.lastLetter+Info.wordLength[1]), 12), 4096);
    mixer.set(finalize64(hash(Info.spaces&0x1FF, grp0,
      (Info.lastUpper<Info.wordLength[0])|
      ((Info.lastUpper<Info.lastLetter+Info.wordLength[1])<<1)|
      ((Info.lastPunct<Info.lastLetter)<<2)|
      ((Info.lastPunct<Info.wordLength[0]+Info.wordGap)<<3)|
      ((Info.lastPunct<Info.lastLetter+Info.wordLength[1]+Info.wordGap)<<4)
    ), 12), 4096);
    mixer.set(finalize64(hash(Info.firstLetter*(Info.wordLength[0]<4), min(6, Info.wordLength[0]), c0), 11), 2048);
    mixer.set(finalize64(hash((*pWord)[0], (*pWord)(0), min(4, Info.wordLength[0]), Info.lastPunct<Info.lastLetter), 11), 2048);
    mixer.set(finalize64(hash(min(4, Info.wordLength[0]), grp0,
      Info.lastUpper<Info.wordLength[0],
      (Info.nestHash>0)?Info.nestHash&0xFF:0x100|(Info.firstLetter*(Info.wordLength[0]>0 && Info.wordLength[0]<4))
    ), 12), 4096);
    mixer.set(finalize64(hash(grp0, Info.masks[4]&0x1F, (Info.masks[4]>>5)&0x1F), 13), 8192);
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

  U8 c = buffer(1), pC=tolower(c), g = (c<0x80)?AsciiGroup[c]:31;
  if(!((g<=4) && g==(Info.asciiMask&0x1f))) //repetition is allowed for groups 0..4
      Info.asciiMask = ((Info.asciiMask<<5) | g)&((U64(1)<<60)-1); //keep last 12 groups (12x5=60 bits)
  Info.masks[4] = Info.asciiMask&((1<<30)-1);
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
      Info.numHashes[0] = combine64(Info.numHashes[0], c);
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
  Info.maskPunct = (BytePos[(unsigned char)',']>BytePos[(unsigned char)'.'])|((BytePos[(unsigned char)',']>BytePos[(unsigned char)'!'])<<1)|((BytePos[(unsigned char)',']>BytePos[(unsigned char)'?'])<<2)|((BytePos[(unsigned char)',']>BytePos[(unsigned char)':'])<<3)|((BytePos[(unsigned char)',']>BytePos[(unsigned char)';'])<<4);
  if (Stats) {
    Stats->Text.state = State;
    Stats->Text.lastPunct = std::min<U32>(0x1F, Info.lastPunct);
    Stats->Text.wordLength = std::min<U32>(0xF, Info.wordLength[0]);
    Stats->Text.boolmask = (Info.lastDigit<Info.wordLength[0]+Info.wordGap)|
                          ((Info.lastUpper<Info.lastLetter+Info.wordLength[1])<<1)|
                          ((Info.lastPunct<Info.wordLength[0]+Info.wordGap)<<2)|
                          ((Info.lastUpper<Info.wordLength[0])<<3);
    Stats->Text.firstLetter = Info.firstLetter;
    Stats->Text.mask = Info.masks[1]&0xFF;
  }
}

void TextModel::SetContexts(Buf& buffer, ModelStats *Stats) {
  const U8 c = buffer(1), lc = tolower(c), m2 = Info.masks[2]&0xF, column = min(0xFF, Info.lastNewLine);;
  const U16 w = ((State==Parse::ReadingWord)?cWord->Hash[1]:pWord->Hash[1])&0xFFFF;
  const U32 h = ((State==Parse::ReadingWord)?cWord->Hash[1]:pWord->Hash[2])*271+c;
  U64 i = State<<6;

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
  Map.set(Info.masks[4]); //last 6 groups
  Map.set(hash(U32(Info.asciiMask),U32(Info.asciiMask>>32))); //last 12 groups
  Map.set(Info.asciiMask & ((1<<20)-1)); //last 4 groups
  Map.set(Info.asciiMask & ((1<<10)-1)); //last 2 groups
  Map.set(hash((Info.asciiMask>>5) &((1<<30)-1),buffer(1)));
  Map.set(hash((Info.asciiMask>>10)&((1<<30)-1),buffer(1),buffer(2)));
  Map.set(hash((Info.asciiMask>>15)&((1<<30)-1),buffer(1),buffer(2),buffer(3)));
}

class MatchModel {
private:
  enum Parameters : U32 {
    MaxLen = 0xFFFF, // longest allowed match
    MaxExtend = 0,   // longest allowed match expansion // warning: larger value -> slowdown
    MinLen = 5,      // minimum required match length
    StepSize = 2,    // additional minimum length increase per higher order hash
    DeltaLen = 5,    // minimum length to switch to delta mode
    NumCtxs = 3,     // number of contexts used
    NumHashes = 3    // number of hashes used
  };
  Array<U32> Table;
  StateMap32 **StateMaps;
  SmallStationaryContextMap **SCM;
  StationaryMap **Maps;
  IndirectContext<U8> iCtx;
  U32 hashes[NumHashes];
  U32 ctx[NumCtxs];
  U32 length;    // rebased length of match (length=1 represents the smallest accepted match length), or 0 if no match
  U32 index;     // points to next byte of match in buffer, 0 when there is no match
  const U32 mask;
  const int hashbits;
  U8 expectedByte; // prediction is based on this byte (buffer[index]), valid only when length>0
  bool delta;
  void Update(Buf& buffer, ModelStats *Stats = nullptr) {
    delta = false;
    // update hashes
    for (U32 i=0, minLen=MinLen+(NumHashes-1)*StepSize; i<NumHashes; i++, minLen-=StepSize) {
      U64 hash = 0;
      for (U32 j=minLen; j>0; j--)
        hash = combine64(hash, buffer(j));
      hashes[i] = finalize64(hash, hashbits);
    }
    // extend current match, if available
    if (length) {
      index++;
      if (length<MaxLen)
        length++;
    }
    // or find a new match, starting with the highest order hash and falling back to lower ones
    else {
      U32 minLen = MinLen+(NumHashes-1)*StepSize, bestLen = 0, bestIndex = 0;
      for (U32 i=0; i<NumHashes && length<minLen; i++, minLen-=StepSize) {
        index = Table[hashes[i]];
        if (index>0) {
          length = 0;
          while (length<(minLen+MaxExtend) && buffer(length+1)==buffer[index-length-1])
            length++;
          if (length>bestLen) {
            bestLen = length;
            bestIndex = index;
          }
        }
      }
      if (bestLen>=MinLen) {
        length = bestLen-(MinLen-1); // rebase, a length of 1 actually means a length of MinLen
        index = bestIndex;
      }
      else
        length = index = 0;
    }
    // update position information in hashtable
    for (U32 i=0; i<NumHashes; i++)
      Table[hashes[i]] = pos;
    expectedByte = buffer[index];
    iCtx+=y, iCtx=(buffer(1)<<8)|expectedByte;
    SCM[0]->set(expectedByte);
    SCM[1]->set(expectedByte);
    SCM[2]->set(pos);
    Maps[0]->set_direct((expectedByte<<8)|buffer(1));
    Maps[1]->set(hash(expectedByte, c0, buffer(1), buffer(2), min(3,(int)ilog2(length+1))));
    Maps[2]->set_direct(iCtx());
    if (Stats)
      Stats->Match.expectedByte = (length>0)?expectedByte:0;
  }
public:
  MatchModel(const U32 Size) :
    Table(Size/sizeof(U32)),
    iCtx{19,1},
    hashes{ 0 },
    ctx{ 0 },
    length(0),
    mask(Size/sizeof(U32)-1),
    hashbits(ilog2(mask+1)),
    expectedByte(0),
    delta(false)
  {
    StateMaps = new StateMap32*[NumCtxs];
    StateMaps[0] = new StateMap32(56*256);
    StateMaps[1] = new StateMap32(8*256*256+1);
    StateMaps[2] = new StateMap32(256*256);
    SCM = new SmallStationaryContextMap*[3];
    SCM[0] = new SmallStationaryContextMap(8,8);
    SCM[1] = new SmallStationaryContextMap(11,1);
    SCM[2] = new SmallStationaryContextMap(8,8);
    Maps = new StationaryMap*[3];
    Maps[0] = new StationaryMap(16,8);
    Maps[1] = new StationaryMap(22,1);
    Maps[2] = new StationaryMap(4,1);
  }
  ~MatchModel(){
    for (U32 i=0; i<NumCtxs; i++)
      delete StateMaps[i];
    delete[] StateMaps;
    for (U32 i=0; i<3; i++)
      delete SCM[i];
    delete[] SCM;
    for (U32 i=0; i<3; i++)
      delete Maps[i];
    delete[] Maps;
  }
  int Predict(Mixer& mixer, Buf& buffer, ModelStats *Stats = nullptr) {
    if (bpos==0)
      Update(buffer, Stats);
    else{
      const U8 B = c0<<(8-bpos);
      SCM[1]->set((bpos<<8)|(expectedByte^B));
      Maps[1]->set(hash(expectedByte, c0, buffer(1), buffer(2), min(3,(int)ilog2(length+1))));
      iCtx+=y, iCtx=(bpos<<16)|(buffer(1)<<8)|(expectedByte^B);
      Maps[2]->set_direct(iCtx());
    }
    const int expectedBit = (expectedByte>>(7-bpos))&1;

    if(length>0) {
      const bool isMatch = (bpos==0)?(buffer(1)==buffer[index-1]):(((expectedByte+256)>>(8-bpos))==c0); // next bit matches the prediction?
      if(!isMatch) {
        delta = (length+MinLen)>DeltaLen;
        length = 0;
      }
    }

    for (U32 i=0; i<NumCtxs; i++)
      ctx[i] = 0;
    if (length>0) {
      if (length<=16)
        ctx[0] = (length-1)*2 + expectedBit; // 0..31
      else
        ctx[0] = 24 + (min(length-1, 63)>>2)*2 + expectedBit; // 32..55
      ctx[0] = ((ctx[0]<<8) | c0);
      ctx[1] = ((expectedByte<<11) | (bpos<<8) | buffer(1)) + 1;
      const int sign = 2*expectedBit-1;
      mixer.add(sign*(min(length,32)<<5)); // +/- 32..1024
      mixer.add(sign*(ilog(length)<<2));   // +/-  0..1024
    }
    else { // no match at all or delta mode
      mixer.add(0);
      mixer.add(0);
    }

    if (delta)
      ctx[2] = (expectedByte<<8) | c0;

    for (U32 i=0; i<NumCtxs; i++) {
      const U32 c = ctx[i];
      const U32 p = StateMaps[i]->p(c);
      if (c!=0)
        mixer.add((stretch(p)+1)>>1);
      else
        mixer.add(0);
    }

    SCM[0]->mix(mixer);
    SCM[1]->mix(mixer, 6);
    SCM[2]->mix(mixer, 5);
    Maps[0]->mix(mixer, 1, 4, 255);
    Maps[1]->mix(mixer);
    Maps[2]->mix(mixer);
    
    if (Stats)
      Stats->Match.length = length;
    return length;
  }
};

class SparseMatchModel {
private:
  enum Parameters : U32 {
    MaxLen    = 0xFFFF, // longest allowed match
    MinLen    = 3,      // default minimum required match length
    NumHashes = 4,      // number of hashes used
  };
  struct sparseConfig {
    U32 offset    = 0;      // number of last input bytes to ignore when searching for a match
    U32 stride    = 1;      // look for a match only every stride bytes after the offset
    U32 deletions = 0;      // when a match is found, ignore these many initial post-match bytes, to model deletions
    U32 minLen    = MinLen;
    U32 bitMask   = 0xFF;   // match every byte according to this bit mask
  };
  Array<U32> Table;
  StationaryMap **Maps;
  IndirectContext<U8> iCtx8;
  IndirectContext<U16> iCtx16;
  MTFList list;
  sparseConfig sparse[NumHashes];
  U32 hashes[NumHashes];
  U32 hashIndex;   // index of hash used to find current match
  U32 length;      // rebased length of match (length=1 represents the smallest accepted match length), or 0 if no match
  U32 index;       // points to next byte of match in buffer, 0 when there is no match
  const U32 mask;
  const int hashbits;
  U8 expectedByte; // prediction is based on this byte (buffer[index]), valid only when length>0
  bool valid;
  void Update(Buf& buffer, ModelStats *Stats = nullptr) {
    // update sparse hashes
    for (U32 i=0; i<NumHashes; i++) {
      U64 hash = 0;
      for (U32 j=0, k=sparse[i].offset+1; j<sparse[i].minLen; j++, k+=sparse[i].stride)
        hash = combine64(hash, buffer(k)&sparse[i].bitMask);
      hashes[i] = finalize64(hash, hashbits);
    }
    // extend current match, if available
    if (length) {
      index++;
      if (length<MaxLen)
        length++;
    }
    // or find a new match
    else {     
      for (int i=list.GetFirst(); i>=0; i=list.GetNext()) {
        index = Table[hashes[i]];
        if (index>0) {
          U32 offset = sparse[i].offset+1;
          while (length<sparse[i].minLen && ((buffer(offset)^buffer[index-offset])&sparse[i].bitMask)==0) {
            length++;
            offset+=sparse[i].stride;
          }
          if (length>=sparse[i].minLen) {
            length-=(sparse[i].minLen-1);
            index+=sparse[i].deletions;
            hashIndex = i;
            list.MoveToFront(i);
            break;
          }
        }
        length = index = 0;
      }
    }
    // update position information in hashtable
    for (U32 i=0; i<NumHashes; i++)
      Table[hashes[i]] = pos;
    
    expectedByte = buffer[index];
    if (valid)
      iCtx8+=y, iCtx16+=buffer(1);
    valid = length>1; // only predict after at least one byte following the match
    if (valid) {
      Maps[0]->set(hash(expectedByte, c0, buffer(1), buffer(2), ilog2(length+1)*NumHashes+hashIndex));
      Maps[1]->set_direct((expectedByte<<8)|buffer(1));
      iCtx8=(buffer(1)<<8)|expectedByte, iCtx16=(buffer(1)<<8)|expectedByte;
      Maps[2]->set_direct(iCtx8());
      Maps[3]->set_direct(iCtx16());
    }
  }
public:
  SparseMatchModel(const U64 Size, const bool AllowBypass = false) :
    Table(Size/sizeof(U32)),
    iCtx8{19,1},
    iCtx16{16},
    list(NumHashes),
    hashes{ 0 },
    hashIndex(0),
    length(0),
    mask(Size/sizeof(U32)-1),
    hashbits(ilog2(mask+1)),
    expectedByte(0),
    valid(false)
  {
    Maps = new StationaryMap*[4];
    Maps[0] = new StationaryMap(22,1);
    Maps[1] = new StationaryMap(14,4);
    Maps[2] = new StationaryMap(8,1);
    Maps[3] = new StationaryMap(19,1);
    sparse[0].minLen=5, sparse[0].bitMask=0xDF;
    sparse[1].offset=1, sparse[1].minLen=4;
    sparse[2].stride=2, sparse[2].minLen=4, sparse[2].bitMask=0xDF;
    sparse[3].minLen=5, sparse[3].bitMask=0xF;
  }
  ~SparseMatchModel(){
    for (U32 i=0; i<4; i++)
      delete Maps[i];
    delete[] Maps;
  }
  int Predict(Mixer& mixer, Buf& buffer, ModelStats *Stats = nullptr) {
    const U8 B = c0<<(8-bpos);
    if (bpos==0)
      Update(buffer, Stats);
    else if (valid) {
      Maps[0]->set(hash(expectedByte, c0, buffer(1), buffer(2), ilog2(length+1)*NumHashes+hashIndex));
      if (bpos==4)
        Maps[1]->set_direct(0x10000|((expectedByte^U8(c0<<4))<<8)|buffer(1));
      iCtx8+=y, iCtx8=(bpos<<16)|(buffer(1)<<8)|(expectedByte^B);
      Maps[2]->set_direct(iCtx8());
      Maps[3]->set_direct((bpos<<16)|(iCtx16()^U32(B|(B<<8))));
    }

    // check if next bit matches the prediction, accounting for the required bitmask
    if (length>0 && (((expectedByte^B)&sparse[hashIndex].bitMask)>>(8-bpos))!=0)
      length = 0;

    if (valid) {
      if (length>1 && ((sparse[hashIndex].bitMask>>(7-bpos))&1)>0) {
        const int expectedBit = (expectedByte>>(7-bpos))&1;
        const int sign = 2*expectedBit-1;
        mixer.add(sign*(min(length-1, 64)<<4)); // +/- 16..1024
        mixer.add(sign*(1<<min(length-2, 3))*min(length-1, 8)<<4); // +/- 16..1024
        mixer.add(sign*512);
      }
      else {
        mixer.add(0); mixer.add(0); mixer.add(0);
      }

      for (int i=0;i<4;i++)
        Maps[i]->mix(mixer, 1, 2);
    }
    else
      for (int i=0; i<11; i++, mixer.add(0));

    mixer.set((hashIndex<<6)|(bpos<<3)|min(7, length), NumHashes*64);
    mixer.set((hashIndex<<11)|(min(7, ilog2(length+1))<<8)|(c0^(expectedByte>>(8-bpos))), NumHashes*2048);

    return length;
  }
};

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
    static U64 word0=0, word1=0, word2=0, word3=0, word4=0, word5=0;
    static U32 wrdhsh=0;
    static U64 xword0=0,xword1=0,xword2=0,cword0=0,ccword=0;
    static U64 number0=0, number1=0;
    static U32 text0=0,data0=0,type0=0;
    static U32 lastLetter=0, firstLetter=0, lastUpper=0, lastDigit=0, wordGap=0;
    static ContextMap cm(MEM()*16, 61);
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
            if (c>4) word0=combine64(word0, c);
            text0=text0*997*16+c;
            wordlen++;
            wordlen=min(wordlen,45);
            f=0;
            w=U32(word0)&(wpos.size()-1);
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
                word0=0;
                wordlen=0;
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
            number0=combine64(number0, c);
            lastDigit = 0;
        }
        else if (number0) {
            type0 = (type0<<2)|2;
            number1=number0;
            number0=0,ccword=0;
        }
        if (!((c>='a' && c<='z') ||(c>='0' && c<='9') || (c>=128 ))){
            data0^=combine64(data0, c);
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
         cm.set(hash(266,h, word1,word2,lastUpper<wordlen));
         cm.set(hash(267,text0&0xffffff,0));
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
     ,type0&0xFFF));
    }
    if (wordlen1)    cm.set(hash(col,wordlen1,above&0x5F,c4&0x5F)); else cm.set(0); //wordlist
    if (wrdhsh)  cm.set(hash(mask2&0x3F, wrdhsh&0xFFF, (0x100|firstLetter)*(wordlen<6),(wordGap>4)*2+(wordlen1>5)) ); else cm.set(0);
    if ( lastLetter<16) cm.set(hash((*pWord).Hash[2], h)); else cm.set(0);
    }
    cm.mix(m);
}

void nestModel(Mixer& m)
{
  static int ic=0, bc=0, pc=0, qc=0, lvc=0, ac=0, ec=0, uc=0, sense1=0, sense2=0, w=0;
  static unsigned int vc=0, wc=0;
  static ContextMap cm(MEM()/2, 12);

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
    U64 i=0;
    cm.set(hash(++i, (vv>0 && vv<3)?0:(lc|0x100), ic&0x3FF, ec&0x7, ac&0x7, uc ));
    cm.set(hash(++i, ic, w, ilog2(bc+1)));
    cm.set(hash(++i, (3*vc+77*pc+373*ic+qc)&0xffff));
    cm.set(hash(++i, (31*vc+27*pc+281*qc)&0xffff));
    cm.set(hash(++i, (13*vc+271*ic+qc+bc)&0xffff));
    cm.set(hash(++i, (17*pc+7*ic)&0xffff));
    cm.set(hash(++i, (13*vc+ic)&0xffff));
    cm.set(hash(++i, (vc/3+pc)&0xffff));
    cm.set(hash(++i, (7*wc+qc)&0xffff));
    cm.set(hash(++i, vc&0xffff, f4&0xf));
    cm.set(hash(++i, (3*pc)&0xffff, f4&0xf));
    cm.set(hash(++i, ic&0xffff, f4&0xf));
  }
  cm.mix(m);
}

inline U8 Clip(int Px){
  return min(0xFF,max(0,Px));
}
inline U8 Clamp4( int Px, U8 n1, U8 n2, U8 n3, U8 n4){
  return min( max(n1,max(n2,max(n3,n4))), max( min(n1,min(n2,min(n3,n4))), Px ));
}

inline U8 LogMeanDiffQt(const U8 a, const U8 b, const U8 limit = 7){
  return (a!=b)?((a>b)<<3)|min(limit,ilog2((a+b)/max(2,abs(a-b)*2)+1)):0;
}
inline U32 LogQt(const U8 Px, const U8 bits){
  return (U32(0x100|Px))>>max(0,(int)(ilog2(Px)-bits));
}

struct dBASE {
  U8 Version;
  U32 nRecords;
  U16 RecordLength, HeaderLength;
  int Start, End;
};

void recordModel(Mixer& m, Filetype filetype, ModelStats *Stats = nullptr) {
  static int cpos1[256] , cpos2[256], cpos3[256], cpos4[256];
  static int wpos1[0x10000];
  static int rlen[3] = {2,3,4}; // run length and 2 candidates
  static int rcount[2] = {0,0}; // candidate counts
  static U8 padding = 0; // detected padding byte
  static U8 N=0, NN=0, NNN=0, NNNN=0, WxNW=0;
  static int prevTransition = 0, nTransition = 0; // position of the last padding transition
  static int col = 0, mxCtx = 0, x = 0;
  static ContextMap cm(32768, 3), cn(32768/2, 3), co(32768*2, 3), cp(MEM(), 16);
  static const int nMaps = 6;
  static StationaryMap Maps[nMaps] ={ {10,8},{10,8},{8,8},{8,8},{8,8},{11,1} };
  static SmallStationaryContextMap sMap[3]{ {11, 1}, {3, 1}, {19,1} };
  static IndirectMap iMap[3]{{8,8},{8,8},{8,8}};
  static bool MayBeImg24b = false;
  static dBASE dbase {};
  static const int nIndCtxs = 5;
  static IndirectContext<U16> iCtx[nIndCtxs]{ {16,8}, {16,8}, {16,8}, {20,8}, {11,1} };

  if (!bpos) {
    int w=c4&0xffff, c=w&255, d=w>>8;
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

    col=pos%rlen[0];
    x = min(0x1F,col/max(1,rlen[0]/32));
    N = buf(rlen[0]), NN = buf(rlen[0]*2), NNN = buf(rlen[0]*3), NNNN = buf(rlen[0]*4);
    for (int i=0; i<nIndCtxs-1; iCtx[i]+=c, i++);
    iCtx[0]=(c<<8)|N;
    iCtx[1]=(buf(rlen[0]-1)<<8)|N;
    iCtx[2]=(c<<8)|buf(rlen[0]-1);
    iCtx[3]=finalize64(hash(c, N, buf(rlen[0]+1)), 20);

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
    
    if (!col)
      nTransition = 0;
    if ((((c4>>8) == SPACE*0x010101) && (c != SPACE)) || (!(c4>>8) && c && ((padding != SPACE) || (pos-prevTransition > rlen[0])))){
      prevTransition = pos;
      nTransition+=(nTransition<31);
      padding = (U8)d;
    }
    
    U64 i=0;

    cm.set(hash(++i, c<<8 | (min(255, pos-cpos1[c])>>2)));
    cm.set(hash(++i, w<<9 | llog(pos-wpos1[w])>>2));
    cm.set(hash(++i, rlen[0] | N<<10 | NN<<18));

    cn.set(hash(++i, w | rlen[0]<<16));
    cn.set(hash(++i, d | rlen[0]<<8));
    cn.set(hash(++i, c | rlen[0]<<8));

    co.set(hash(++i, c<<8|min(255, pos-cpos1[c])));
    co.set(hash(++i, c<<17 | d<<9|llog(pos-wpos1[w])>>2));
    co.set(hash(++i, c<<8 | N));

    cp.set(hash(++i, rlen[0] | N<<10 | col<<18));
    cp.set(hash(++i, rlen[0] | c<<10 | col<<18));
    cp.set(hash(++i, col | rlen[0]<<12));

    if (rlen[0]>8){
      cp.set(hash(++i, min(min(0xFF,rlen[0]),pos-prevTransition), min(0x3FF,col), (w&0xF0F0)|(w==((padding<<8)|padding)), nTransition ) );
      cp.set(hash(++i, w, (buf(rlen[0]+1)==padding && N==padding), col/max(1,rlen[0]/32) ) );
    }
    else
      cp.set(0), cp.set(0);

    cp.set(hash(++i, N|((NN&0xF0)<<4)|((NNN&0xE0)<<7)|((NNNN&0xE0)<<10)|((col/max(1,rlen[0]/16))<<18) ));
    cp.set(hash(++i, (N&0xF8)|((NN&0xF8)<<8)|(col<<16) ));
    cp.set(hash(++i, N, NN));

    cp.set(hash(++i, col, iCtx[0]()));
    cp.set(hash(++i, col, iCtx[1]()));
    cp.set(hash(++i, col, iCtx[0]()&0xFF, iCtx[1]()&0xFF));

    cp.set(hash(++i, iCtx[2]()));
    cp.set(hash(++i, iCtx[3]()));
    cp.set(hash(++i, iCtx[1]()&0xFF, iCtx[3]()&0xFF));

    cp.set(hash(++i, N, (WxNW=c^buf(rlen[0]+1))));
    cp.set(hash(++i, (Stats!=nullptr && Stats->Match.length>0)?Stats->Match.expectedByte:0x100|U8(iCtx[1]()), N, WxNW));

    int k=0x300;
    if (MayBeImg24b)
      k = (col%3)<<8, Maps[0].set_direct(Clip(((U8)(c4>>16))+c-(c4>>24))|k);
    else
      Maps[0].set_direct(Clip(c*2-d)|k);
    Maps[1].set_direct(Clip(c+N-buf(rlen[0]+1))|k);
    Maps[2].set_direct(Clip(N+NN-NNN));
    Maps[3].set_direct(Clip(N*2-NN));
    Maps[4].set_direct(Clip(N*3-NN*3+NNN));
    iMap[0].set_direct(N+NN-NNN);
    iMap[1].set_direct(N*2-NN);
    iMap[2].set_direct(N*3-NN*3+NNN);

    // update last context positions
    cpos4[c]=cpos3[c];
    cpos3[c]=cpos2[c];
    cpos2[c]=cpos1[c];
    cpos1[c]=pos;
    wpos1[w]=pos;

    mxCtx = (rlen[0]>128)?(min(0x7F,col/max(1,rlen[0]/128))):col;
  }
  U8 B = c0<<(8-bpos);
  U32 ctx = (N^B)|(bpos<<8);
  iCtx[nIndCtxs-1]+=y, iCtx[nIndCtxs-1]=ctx;
  Maps[nMaps-1].set_direct(ctx);
  sMap[0].set(ctx);
  sMap[1].set(iCtx[nIndCtxs-1]());
  sMap[2].set((ctx<<8)|WxNW);
  
  cm.mix(m);
  cn.mix(m);
  co.mix(m);
  cp.mix(m);
  for (int i=0;i<nMaps;i++)
    Maps[i].mix(m, 1, 3);
  for (int i=0; i<3; i++)
    iMap[i].mix(m, 1, 3, 255);
  sMap[0].mix(m, 6, 1, 3);
  sMap[1].mix(m, 6, 1, 3);
  sMap[2].mix(m, 5, 1, 2);

  m.set( (rlen[0]>2)*( (bpos<<7)|mxCtx ), 1024 );
  m.set( ((N^B)>>4)|(x<<4), 512 );
  m.set( (grp0<<5)|x, 11*32);
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

void linearPredictionModel(Mixer& m) {
  static const int nOLS=3, nLnrPrd=nOLS+2;
  static SmallStationaryContextMap sMap[nLnrPrd]{ {11,1},{11,1},{11,1},{11,1},{11,1} };
  static OLS<double, U8> ols[nOLS]{ {32, 4, 0.995}, {32, 4, 0.995}, {32, 4, 0.995} };
  static U8 prd[nLnrPrd]{ 0 };

  if (bpos==0) {
    const U8 W=buf(1), WW=buf(2), WWW=buf(3);
    int i=0;
    for (; i<nOLS; i++)
      ols[i].Update(W);
    for (i=1; i<=32; i++) {
      ols[0].Add(buf(i));
      ols[1].Add(buf(i*2-1));
      ols[2].Add(buf(i*2));
    }
    for (i=0; i<nOLS; i++)
      prd[i]=Clip(floor(ols[i].Predict()));
    prd[i++]=Clip(W*2-WW);
    prd[i  ]=Clip(W*3-WW*3+WWW);
  }
  const U8 B=c0<<(8-bpos);
  for (int i=0; i<nLnrPrd; i++) {
    sMap[i].set((prd[i]-B)*8+bpos);
    sMap[i].mix(m, 6, 1, 2);
  }
}

void sparseModel(Mixer& m, int seenbefore, int howmany) {
  static ContextMap cm(MEM()*2, 40+2);
  if (bpos==0) {
    U64 i=0;
    cm.set(hash(++i,seenbefore));
    cm.set(hash(++i,howmany));
    cm.set(hash(++i,buf(1)|buf(5)<<8));
    cm.set(hash(++i,buf(1)|buf(6)<<8));
    cm.set(hash(++i,buf(3)|buf(6)<<8));
    cm.set(hash(++i,buf(4)|buf(8)<<8));
    cm.set(hash(++i,buf(1)|buf(3)<<8|buf(5)<<16));
    cm.set(hash(++i,buf(2)|buf(4)<<8|buf(6)<<16));
    cm.set(hash(++i,c4&0x00f0f0ff));
    cm.set(hash(++i,c4&0x00ff00ff));
    cm.set(hash(++i,c4&0xff0000ff));
    cm.set(hash(++i,c4&0x00f8f8f8));
    cm.set(hash(++i,c4&0xf8f8f8f8));
    cm.set(hash(++i,f4&0x00000fff));
    cm.set(hash(++i,f4));
    cm.set(hash(++i,c4&0x00e0e0e0));
    cm.set(hash(++i,c4&0xe0e0e0e0));
    cm.set(hash(++i,c4&0x810000c1));
    cm.set(hash(++i,c4&0xC3CCC38C));
    cm.set(hash(++i,c4&0x0081CC81));
    cm.set(hash(++i,c4&0x00c10081));
    for (int j=1; j<8; ++j) {
      cm.set(hash(++i,seenbefore|buf(j)<<8));
      cm.set(hash(++i,(buf(j+2)<<8)|buf(j+1)));
      cm.set(hash(++i,(buf(j+3)<<8)|buf(j+1)));
    }
  }
  cm.mix(m);
}

U32 x4=0;
void sparseModel1(Mixer& m, int seenbefore, int howmany) {
   static ContextMap cm(MEM()*4, 31);
    static SmallStationaryContextMap scm1(7,8), scm2(8,8), scm3(4,8),
     scm4(6,8), scm5(4,8),scm6(4,8), scma(7,8);
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
  static ContextMap cm(MEM(), 3);
  if( bpos == 0 ){
    static int pos00=0,pos20=0,posnl=0;
    int c=c4&0xff;
    if(c==0x00)pos00=pos;
    if(c==0x20)pos20=pos;
    if(c==0xff||c=='\r'||c=='\n')posnl=pos;
    U64 i=0;
    cm.set(hash(++i, min(pos-pos00,255) | c<<8));
    cm.set(hash(++i, min(pos-pos20,255) | c<<8));
    cm.set(hash(++i, min(pos-posnl,255) | c<<8));
  }
  cm.mix(m);
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
  const int S=14; // number of contexts
  static U8* cp[S];
  static StateMap sm[S];
  static StateMap32 map(16);
  static U8 WW=0, W=0, NWW=0, NW=0, N=0, NE=0, NEE=0, NNWW = 0, NNW=0, NN=0, NNE=0, NNEE=0;
  static int col=0, line=0, run=0, prevColor=0, px=0;
  if (!cp[0]){
    for (int i=0;i<S;i++)
      cp[i]=t[263*i]+1;
  }
  for (int i=0;i<S;i++)
    *cp[i]=nex(*cp[i],y);

  if (!bpos || bpos==4){
    WW=W, NWW=NW, NW=N, N=NE, NE=NEE, NNWW=NWW, NNW=NN, NN=NNE, NNE=NNEE;
    if (!bpos)
      W=c4&0xF, NEE=buf(w-1)>>4, NNEE=buf(w*2-1)>>4;
    else
      W=c0&0xF, NEE=buf(w-1)&0xF, NNEE=buf(w*2-1)&0xF;
    run=(W!=WW || !col)?(prevColor=WW,0):min(0xFFF,run+1);
    px=1;
    U64 i=0; cp[i]=t[hash(i,W,NW,N)];
        i++; cp[i]=t[hash(i,N, min(0xFFF, col/8))];
        i++; cp[i]=t[hash(i,W,NW,N,NN,NE)];
        i++; cp[i]=t[hash(i,W, N, NE+NNE*16, NEE+NNEE*16)];
        i++; cp[i]=t[hash(i,W, N, NW+NNW*16, NWW+NNWW*16)];
        i++; cp[i]=t[hash(i,W, ilog2(run+1), prevColor, col/max(1,w/2) )];
        i++; cp[i]=t[hash(i,NE, min(0x3FF, (col+line)/max(1,w*8)))];
        i++; cp[i]=t[hash(i,NW, (col-line)/max(1,w*8))];
        i++; cp[i]=t[hash(i,WW*16+W,NN*16+N,NNWW*16+NW)];
        i++; cp[i]=t[hash(i,N,NN)];
        i++; cp[i]=t[hash(i,W,WW)];
        i++; cp[i]=t[hash(i,W,NE)];
        i++; cp[i]=t[hash(i,WW,NN,NEE)];
        i++; cp[i]=t[-1];
    ++col;
    col*=col<w*2;
    line+=(!col);
  }
  else{
    px+=px+y;
    int j=(y+1)<<(bpos&3);
    for (int i=0;i<S;i++)
      cp[i]+=j;
  }

  // predict
  for (int i=0; i<S; i++) {
    const U8 s = *cp[i];
    const int n0=-!nex(s, 2), n1=-!nex(s, 3);
    const int p1 = sm[i].p(s);
    const int st = stretch(p1)>>1;
    m.add(st);
    m.add((p1-2047)>>2);
    m.add(st*abs(n1-n0));
  }
  m.add(stretch(map.p(px))>>1);

  m.set(W*16+px, 256);
  m.set(min(31,col/max(1,w/16))+N*32, 512);
  m.set((bpos&3)+4*W+64*min(7,ilog2(run+1)), 512);
  m.set(W+NE*16+(bpos&3)*256, 1024);
  m.set(px, 16);
  m.set(0,1);
}

void im8bitModel(Mixer& m, int w, ModelStats *Stats = nullptr, int gray = 0) {
  static const int nOLS = 5;
  static const int nMaps0 = 2;
  static const int nMaps1 = 55;
  static const int nMaps = nMaps0 + nMaps1 + nOLS;
  static const int nPltMaps = 4;
  static ContextMap cm(MEM()*4, 48 + nPltMaps);
  static StationaryMap Map[nMaps] = {{ 0,8}, {15,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}};
  static SmallStationaryContextMap pltMap[nPltMaps] = { {11,1},{11,1},{11,1},{11,1} };
  static IndirectContext<U8> iCtx[nPltMaps] = { {16,8}, {16,8}, {16,8}, {16,8} };
  //pixel neighborhood
  static U8 WWWWWW, WWWWW, WWWW, WWW, WW, W;
  static U8 NWWWW, NWWW, NWW, NW, N, NE, NEE, NEEE, NEEEE;
  static U8 NNWWW, NNWW, NNW, NN, NNE, NNEE, NNEEE;
  static U8 NNNWW, NNNW, NNN, NNNE, NNNEE;
  static U8 NNNNW, NNNN, NNNNE;
  static U8 NNNNN;
  static U8 NNNNNN;
  static int ctx, lastPos=0, col=0, x=0, line=0;
  static int columns[2] = {1,1}, column[2];
  static U8 MapCtxs[nMaps1] = { 0 }, pOLS[nOLS] = { 0 };
  static const double lambda[nOLS] ={ 0.996, 0.87, 0.93, 0.8, 0.9 };
  static const int num[nOLS] ={ 32, 12, 15, 10, 14 };
  static OLS<double, U8> ols[nOLS] = { 
    {num[0], 1, lambda[0]},
    {num[1], 1, lambda[1]},
    {num[2], 1, lambda[2]},
    {num[3], 1, lambda[3]},
    {num[4], 1, lambda[4]}
  };
  static const U8 *ols_ctx1[32] = { &WWWWWW, &WWWWW, &WWWW, &WWW, &WW, &W, &NWWWW, &NWWW, &NWW, &NW, &N, &NE, &NEE, &NEEE, &NEEEE, &NNWWW, &NNWW, &NNW, &NN, &NNE, &NNEE, &NNEEE, &NNNWW, &NNNW, &NNN, &NNNE, &NNNEE, &NNNNW, &NNNN, &NNNNE, &NNNNN, &NNNNNN };
  static const U8 *ols_ctx2[12] = { &WWW, &WW, &W, &NWW, &NW, &N, &NE, &NEE, &NNW, &NN, &NNE, &NNN }; 
  static const U8 *ols_ctx3[15] = { &N, &NE, &NEE, &NEEE, &NEEEE, &NN, &NNE, &NNEE, &NNEEE, &NNN, &NNNE, &NNNEE, &NNNN, &NNNNE, &NNNNN };
  static const U8 *ols_ctx4[10] = { &N, &NE, &NEE, &NEEE, &NN, &NNE, &NNEE, &NNN, &NNNE, &NNNN };
  static const U8 *ols_ctx5[14] = { &WWWW, &WWW, &WW, &W, &NWWW, &NWW, &NW, &N, &NNWW, &NNW, &NN, &NNNW, &NNN, &NNNN };
  static const U8 **ols_ctxs[nOLS] = { &ols_ctx1[0], &ols_ctx2[0], &ols_ctx3[0], &ols_ctx4[0], &ols_ctx5[0] };
  // Select nearby pixels as context
  if (!bpos) {
    if (pos!=lastPos+1){
      x = line = 0;
      columns[0] = max(1,w/max(1,ilog2(w)*2));
      columns[1] = max(1,columns[0]/max(1,ilog2(columns[0])));
    }
    else{
      ++x;
      x*=x<w;
      line+=(x==0);
    }
    lastPos = pos;
    column[0]=x/columns[0];
    column[1]=x/columns[1];
    U64 i=0;
    WWWWW=buf(5), WWWW=buf(4), WWW=buf(3), WW=buf(2), W=buf(1);
    NWWWW=buf(w+4), NWWW=buf(w+3), NWW=buf(w+2), NW=buf(w+1), N=buf(w), NE=buf(w-1), NEE=buf(w-2), NEEE=buf(w-3), NEEEE=buf(w-4);
    NNWWW=buf(w*2+3), NNWW=buf(w*2+2), NNW=buf(w*2+1), NN=buf(w*2), NNE=buf(w*2-1), NNEE=buf(w*2-2), NNEEE=buf(w*2-3);
    NNNWW=buf(w*3+2), NNNW=buf(w*3+1), NNN=buf(w*3), NNNE=buf(w*3-1), NNNEE=buf(w*3-2);
    NNNNW=buf(w*4+1), NNNN=buf(w*4), NNNNE=buf(w*4-1);
    NNNNN=buf(w*5);
    NNNNNN=buf(w*6);

    int j = 0;
    MapCtxs[j++] = Clamp4(W+N-NW,W,NW,N,NE);
    MapCtxs[j++] = Clip(W+N-NW);
    MapCtxs[j++] = Clamp4(W+NE-N,W,NW,N,NE);
    MapCtxs[j++] = Clip(W+NE-N);
    MapCtxs[j++] = Clamp4(N+NW-NNW,W,NW,N,NE);
    MapCtxs[j++] = Clip(N+NW-NNW);
    MapCtxs[j++] = Clamp4(N+NE-NNE,W,N,NE,NEE);
    MapCtxs[j++] = Clip(N+NE-NNE);
    MapCtxs[j++] = (W+NEE)/2;
    MapCtxs[j++] = Clip(N*3-NN*3+NNN);
    MapCtxs[j++] = Clip(W*3-WW*3+WWW);
    MapCtxs[j++] = (W+Clip(NE*3-NNE*3+buf(w*3-1)))/2;
    MapCtxs[j++] = (W+Clip(NEE*3-buf(w*2-3)*3+buf(w*3-4)))/2;
    MapCtxs[j++] = Clip(NN+buf(w*4)-buf(w*6));
    MapCtxs[j++] = Clip(WW+buf(4)-buf(6));
    MapCtxs[j++] = Clip((buf(w*5)-6*buf(w*4)+15*NNN-20*NN+15*N+Clamp4(W*2-NWW,W,NW,N,NN))/6);
    MapCtxs[j++] = Clip((-3*WW+8*W+Clamp4(NEE*3-NNEE*3+buf(w*3-2),NE,NEE,buf(w-3),buf(w-4)))/6);
    MapCtxs[j++] = Clip(NN+NW-buf(w*3+1));
    MapCtxs[j++] = Clip(NN+NE-buf(w*3-1));
    MapCtxs[j++] = Clip((W*2+NW)-(WW+2*NWW)+buf(w+3));
    MapCtxs[j++] = Clip(((NW+NWW)/2)*3-buf(w*2+3)*3+(buf(w*3+4)+buf(w*3+5))/2);
    MapCtxs[j++] = Clip(NEE+NE-buf(w*2-3));
    MapCtxs[j++] = Clip(NWW+WW-buf(w+4));
    MapCtxs[j++] = Clip(((W+NW)*3-NWW*6+buf(w+3)+buf(w*2+3))/2);
    MapCtxs[j++] = Clip((NE*2+NNE)-(NNEE+buf(w*3-2)*2)+buf(w*4-3));
    MapCtxs[j++] = buf(w*6);
    MapCtxs[j++] = (buf(w-4)+buf(w-6))/2;
    MapCtxs[j++] = (buf(4)+buf(6))/2;
    MapCtxs[j++] = (W+N+buf(w-5)+buf(w-7))/4;
    MapCtxs[j++] = Clip(buf(w-3)+W-NEE);
    MapCtxs[j++] = Clip(4*NNN-3*buf(w*4));
    MapCtxs[j++] = Clip(N+NN-NNN);
    MapCtxs[j++] = Clip(W+WW-WWW);
    MapCtxs[j++] = Clip(W+NEE-NE);
    MapCtxs[j++] = Clip(WW+NEE-N);
    MapCtxs[j++] = (Clip(W*2-NW)+Clip(W*2-NWW)+N+NE)/4;
    MapCtxs[j++] = Clamp4(N*2-NN,W,N,NE,NEE);
    MapCtxs[j++] = (N+NNN)/2;
    MapCtxs[j++] = Clip(NN+W-NNW);
    MapCtxs[j++] = Clip(NWW+N-NNWW);
    MapCtxs[j++] = Clip((4*WWW-15*WW+20*W+Clip(NEE*2-NNEE))/10);
    MapCtxs[j++] = Clip((buf(w*3-3)-4*NNEE+6*NE+Clip(W*3-NW*3+NNW))/4);
    MapCtxs[j++] = Clip((N*2+NE)-(NN+2*NNE)+buf(w*3-1));
    MapCtxs[j++] = Clip((NW*2+NNW)-(NNWW+buf(w*3+2)*2)+buf(w*4+3));
    MapCtxs[j++] = Clip(NNWW+W-buf(w*2+3));
    MapCtxs[j++] = Clip((-buf(w*4)+5*NNN-10*NN+10*N+Clip(W*4-NWW*6+buf(w*2+3)*4-buf(w*3+4)))/5);
    MapCtxs[j++] = Clip(NEE+Clip(buf(w-3)*2-buf(w*2-4))-buf(w-4));
    MapCtxs[j++] = Clip(NW+W-NWW);
    MapCtxs[j++] = Clip((N*2+NW)-(NN+2*NNW)+buf(w*3+1));
    MapCtxs[j++] = Clip(NN+Clip(NEE*2-buf(w*2-3))-NNE);
    MapCtxs[j++] = Clip((-buf(4)+5*WWW-10*WW+10*W+Clip(NE*2-NNE))/5);
    MapCtxs[j++] = Clip((-buf(5)+4*buf(4)-5*WWW+5*W+Clip(NE*2-NNE))/4);
    MapCtxs[j++] = Clip((WWW-4*WW+6*W+Clip(NE*3-NNE*3+buf(w*3-1)))/4);
    MapCtxs[j++] = Clip((-NNEE+3*NE+Clip(W*4-NW*6+NNW*4-buf(w*3+1)))/3);
    MapCtxs[j++] = ((W+N)*3-NW*2)/4;
    for (j=0; j<nOLS; j++) {
      ols[j].Update(W);
      pOLS[j] = Clip(int(floor(ols[j].Predict(ols_ctxs[j]))));
    }
    for (j=0; j<nPltMaps; j++)
      iCtx[j]+=W;
    iCtx[0]=W|(NE<<8);
    iCtx[1]=W|(N<<8);
    iCtx[2]=W|(WW<<8);
    iCtx[3]=N|(NN<<8);

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
      cm.set(hash(++i, W, N, NW));
      cm.set(hash(++i, N, NN, NNN));
      cm.set(hash(++i, W, NE, NEE));
      cm.set(hash(++i, W, NW, N, NE));
      cm.set(hash(++i, N, NE, NN, NNE));
      cm.set(hash(++i, N, NW, NNW, NN));
      cm.set(hash(++i, W, WW, NWW, NW));
      cm.set(hash(++i, W, NW, N, WW, NWW));
      cm.set(hash(++i, column[0]));
      cm.set(hash(++i, N, column[1] ));
      cm.set(hash(++i, W, column[1] ));
      cm.set(++i);
      for (int j=0; j<nPltMaps; j++)
        cm.set(hash(++i, iCtx[j]()));

      ctx = min(0x1F,x/min(0x20,columns[0]));
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

      ctx = min(0x1F,x/max(1,w/min(32,columns[0])))|( ( ((abs(W-N)*16>W+N)<<1)|(abs(N-NW)>8) )<<5 )|((W+N)&0x180);
    }
    if (Stats) {
      Stats->Image.pixels.W = W;
      Stats->Image.pixels.N = N;
      Stats->Image.pixels.NN = NN;
      Stats->Image.pixels.WW = WW;
      Stats->Image.ctx = ctx>>gray;
    }
  }
  
  U8 B=(c0<<(8-bpos));
  int i=1;
  Map[i++].set_direct((((U8)(Clip(W+N-NW)-B))*8+bpos)|(LogMeanDiffQt(Clip(N+NE-NNE),Clip(N+NW-NNW))<<11));

  for (int j=0; j<nMaps1; i++, j++)
    Map[i].set_direct((MapCtxs[j]-B)*8+bpos);

  for (int j=0; i<nMaps; i++, j++)
    Map[i].set_direct((pOLS[j]-B)*8+bpos);

  cm.mix(m);
  if (gray){
    for (int i=0;i<nMaps;i++)
      Map[i].mix(m);
  }
  else {
    for (int i=0; i<nPltMaps; i++) {
      pltMap[i].set((bpos<<8)|iCtx[i]());
      pltMap[i].mix(m);
    }
  }
  col=(col+1)&7;
  m.set(ctx, 2048);
  m.set(col, 8);
  m.set((N+W)>>4, 32);
  m.set(c0, 256);
  m.set( ((abs((int)(W-N))>4)<<9)|((abs((int)(N-NE))>4)<<8)|((abs((int)(W-NW))>4)<<7)|((W>N)<<6)|((N>NE)<<5)|((W>NW)<<4)|((W>WW)<<3)|((N>NN)<<2)|((NW>NNWW)<<1)|(NE>NNEE), 1024 );
  m.set(min(63,column[0]), 64);
  m.set(min(127,column[1]), 128);
  m.set(min(255,(x+line)/32), 256);
}

void im24bitModel(Mixer& m, int w, ModelStats *Stats = nullptr, int alpha=0) {
  static const int nMaps0 = 18;
  static const int nMaps1 = 76;
  static const int nOLS = 6;
  static const int nMaps = nMaps0+nMaps1+nOLS;
  static const int nSCMaps = 59;
  static ContextMap cm(MEM()*4, 47);
  static SmallStationaryContextMap SCMap[nSCMaps] = { {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                                      {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                                      {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                                      {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                                      {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                                      {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                                      {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                                      {11,1}, {11,1}, { 0,8} };
  static StationaryMap Map[nMaps] ={ { 8,8}, { 8,8}, { 8,8}, { 2,8}, { 0,8}, {15,1}, {15,1}, {15,1}, {15,1}, {15,1},
                                     {17,1}, {17,1}, {17,1}, {17,1}, {13,1}, {13,1}, {13,1}, {13,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1},
                                     {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}, {11,1}};
  //pixel neighborhood
  static U8 WWWWWW, WWWWW, WWWW, WWW, WW, W;
  static U8 NWWWW, NWWW, NWW, NW, N, NE, NEE, NEEE, NEEEE;
  static U8 NNWWW, NNWW, NNW, NN, NNE, NNEE, NNEEE;
  static U8 NNNWW, NNNW, NNN, NNNE, NNNEE;
  static U8 NNNNW, NNNN, NNNNE;
  static U8 NNNNN;
  static U8 NNNNNN;
  static U8 WWp1, Wp1, p1, NWp1, Np1, NEp1, NNp1;
  static U8 WWp2, Wp2, p2, NWp2, Np2, NEp2, NNp2;
  static int color = -1, stride = 3;
  static int ctx[2], padding, lastPos, x = 0, line = 0;
  static int columns[2] = {1,1}, column[2];
  static U8 MapCtxs[nMaps1] = { 0 }, SCMapCtxs[nSCMaps-1] = { 0 }, pOLS[nOLS] = { 0 };
  static const double lambda[nOLS] ={ 0.98, 0.87, 0.9, 0.8, 0.9, 0.7 };
  static const int num[nOLS] ={ 32, 12, 15, 10, 14, 8 };
  static OLS<double, U8> ols[nOLS][4] = { 
    {{num[0], 1, lambda[0]}, {num[0], 1, lambda[0]}, {num[0], 1, lambda[0]}, {num[0], 1, lambda[0]}},
    {{num[1], 1, lambda[1]}, {num[1], 1, lambda[1]}, {num[1], 1, lambda[1]}, {num[1], 1, lambda[1]}},
    {{num[2], 1, lambda[2]}, {num[2], 1, lambda[2]}, {num[2], 1, lambda[2]}, {num[2], 1, lambda[2]}},
    {{num[3], 1, lambda[3]}, {num[3], 1, lambda[3]}, {num[3], 1, lambda[3]}, {num[3], 1, lambda[3]}},
    {{num[4], 1, lambda[4]}, {num[4], 1, lambda[4]}, {num[4], 1, lambda[4]}, {num[4], 1, lambda[4]}},
    {{num[5], 1, lambda[5]}, {num[5], 1, lambda[5]}, {num[5], 1, lambda[5]}, {num[5], 1, lambda[5]}}
  };
  static const U8 *ols_ctx1[32] = { &WWWWWW, &WWWWW, &WWWW, &WWW, &WW, &W, &NWWWW, &NWWW, &NWW, &NW, &N, &NE, &NEE, &NEEE, &NEEEE, &NNWWW, &NNWW, &NNW, &NN, &NNE, &NNEE, &NNEEE, &NNNWW, &NNNW, &NNN, &NNNE, &NNNEE, &NNNNW, &NNNN, &NNNNE, &NNNNN, &NNNNNN };
  static const U8 *ols_ctx2[12] = { &WWW, &WW, &W, &NWW, &NW, &N, &NE, &NEE, &NNW, &NN, &NNE, &NNN }; 
  static const U8 *ols_ctx3[15] = { &N, &NE, &NEE, &NEEE, &NEEEE, &NN, &NNE, &NNEE, &NNEEE, &NNN, &NNNE, &NNNEE, &NNNN, &NNNNE, &NNNNN };
  static const U8 *ols_ctx4[10] = { &N, &NE, &NEE, &NEEE, &NN, &NNE, &NNEE, &NNN, &NNNE, &NNNN };
  static const U8 *ols_ctx5[14] = { &WWWW, &WWW, &WW, &W, &NWWW, &NWW, &NW, &N, &NNWW, &NNW, &NN, &NNNW, &NNN, &NNNN };
  static const U8 *ols_ctx6[ 8] = { &WWW, &WW, &W, &NNN, &NN, &N, &p1, &p2 };
  static const U8 **ols_ctxs[nOLS] = { &ols_ctx1[0], &ols_ctx2[0], &ols_ctx3[0], &ols_ctx4[0], &ols_ctx5[0], &ols_ctx6[0] };

  // Select nearby pixels as context
  if (!bpos) {
    if ((color < 0) || (pos-lastPos != 1)){
      stride = 3+alpha;
      padding = w%stride;
      x = line = 0;
      columns[0] = max(1,w/max(1,ilog2(w)*3));
      columns[1] = max(1,columns[0]/max(1,ilog2(columns[0])));
    }
    lastPos = pos;
    ++x;
    x*=x<w;
    line+=(x==0);
    if (x+padding<w) {
      ++color;
      color*=color<stride;
    } else {
      color=(padding>0)*(stride+1);
    }

    column[0]=x/columns[0];
    column[1]=x/columns[1];

    WWWWWW=buf(6*stride), WWWWW=buf(5*stride), WWWW=buf(4*stride), WWW=buf(3*stride), WW=buf(2*stride), W=buf(stride);
    NWWWW=buf(w+4*stride), NWWW=buf(w+3*stride), NWW=buf(w+2*stride), NW=buf(w+stride), N=buf(w), NE=buf(w-stride), NEE=buf(w-2*stride), NEEE=buf(w-3*stride), NEEEE=buf(w-4*stride);
    NNWWW=buf(w*2+stride*3), NNWW=buf((w+stride)*2), NNW=buf(w*2+stride), NN=buf(w*2), NNE=buf(w*2-stride), NNEE=buf((w-stride)*2), NNEEE=buf(w*2-stride*3);
    NNNWW=buf(w*3+stride*2), NNNW=buf(w*3+stride), NNN=buf(w*3), NNNE=buf(w*3-stride), NNNEE=buf(w*3-stride*2);
    NNNNW=buf(w*4+stride), NNNN=buf(w*4), NNNNE=buf(w*4-stride);
    NNNNN=buf(w*5);
    NNNNNN=buf(w*6);
    WWp1=buf(stride*2+1), Wp1=buf(stride+1), p1=buf(1), NWp1=buf(w+stride+1), Np1=buf(w+1), NEp1=buf(w-stride+1), NNp1=buf(w*2+1);
    WWp2=buf(stride*2+2), Wp2=buf(stride+2), p2=buf(2), NWp2=buf(w+stride+2), Np2=buf(w+2), NEp2=buf(w-stride+2), NNp2=buf(w*2+2);

    int j = 0;
    MapCtxs[j++] = Clamp4(N+p1-Np1, W, NW, N, NE);
    MapCtxs[j++] = Clamp4(N+p2-Np2, W, NW, N, NE);
    MapCtxs[j++] = (W+Clamp4(NE*3-NNE*3+NNNE, W, N, NE, NEE))/2;
    MapCtxs[j++] = Clamp4((W+Clip(NE*2-NNE))/2, W, NW, N, NE);
    MapCtxs[j++] = (W+NEE)/2;
    MapCtxs[j++] = Clip((WWW-4*WW+6*W+Clip(NE*4-NNE*6+NNNE*4-NNNNE))/4);
    MapCtxs[j++] = Clip((-WWWW+5*WWW-10*WW+10*W+Clamp4(NE*4-NNE*6+NNNE*4-NNNNE, N, NE, NEE, NEEE))/5);
    MapCtxs[j++] = Clip((-4*WW+15*W+10*Clip(NE*3-NNE*3+NNNE)-Clip(NEEE*3-NNEEE*3+buf(w*3-3*stride)))/20);
    MapCtxs[j++] = Clip((-3*WW+8*W+Clamp4(NEE*3-NNEE*3+NNNEE, NE, NEE, NEEE, NEEEE))/6);
    MapCtxs[j++] = Clip((W+Clip(NE*2-NNE))/2+p1-(Wp1+Clip(NEp1*2-buf(w*2-stride+1)))/2);
    MapCtxs[j++] = Clip((W+Clip(NE*2-NNE))/2+p2-(Wp2+Clip(NEp2*2-buf(w*2-stride+2)))/2);
    MapCtxs[j++] = Clip((-3*WW+8*W+Clip(NEE*2-NNEE))/6+p1-(-3*WWp1+8*Wp1+Clip(buf(w-stride*2+1)*2-buf(w*2-stride*2+1)))/6);
    MapCtxs[j++] = Clip((-3*WW+8*W+Clip(NEE*2-NNEE))/6+p2-(-3*WWp2+8*Wp2+Clip(buf(w-stride*2+2)*2-buf(w*2-stride*2+2)))/6);
    MapCtxs[j++] = Clip((W+NEE)/2+p1-(Wp1+buf(w-stride*2+1))/2);
    MapCtxs[j++] = Clip((W+NEE)/2+p2-(Wp2+buf(w-stride*2+2))/2);
    MapCtxs[j++] = Clip((WW+Clip(NEE*2-NNEE))/2+p1-(WWp1+Clip(buf(w-stride*2+1)*2-buf(w*2-stride*2+1)))/2);
    MapCtxs[j++] = Clip((WW+Clip(NEE*2-NNEE))/2+p2-(WWp2+Clip(buf(w-stride*2+2)*2-buf(w*2-stride*2+2)))/2);
    MapCtxs[j++] = Clip(WW+NEE-N+p1-Clip(WWp1+buf(w-stride*2+1)-Np1));
    MapCtxs[j++] = Clip(WW+NEE-N+p2-Clip(WWp2+buf(w-stride*2+2)-Np2));
    MapCtxs[j++] = Clip(W+N-NW);
    MapCtxs[j++] = Clip(W+N-NW+p1-Clip(Wp1+Np1-NWp1));
    MapCtxs[j++] = Clip(W+N-NW+p2-Clip(Wp2+Np2-NWp2));
    MapCtxs[j++] = Clip(W+NE-N);
    MapCtxs[j++] = Clip(N+NW-NNW);
    MapCtxs[j++] = Clip(N+NW-NNW+p1-Clip(Np1+NWp1-buf(w*2+stride+1)));
    MapCtxs[j++] = Clip(N+NW-NNW+p2-Clip(Np2+NWp2-buf(w*2+stride+2)));
    MapCtxs[j++] = Clip(N+NE-NNE);
    MapCtxs[j++] = Clip(N+NE-NNE+p1-Clip(Np1+NEp1-buf(w*2-stride+1)));
    MapCtxs[j++] = Clip(N+NE-NNE+p2-Clip(Np2+NEp2-buf(w*2-stride+2)));
    MapCtxs[j++] = Clip(N+NN-NNN);
    MapCtxs[j++] = Clip(N+NN-NNN+p1-Clip(Np1+NNp1-buf(w*3+1)));
    MapCtxs[j++] = Clip(N+NN-NNN+p2-Clip(Np2+NNp2-buf(w*3+2)));
    MapCtxs[j++] = Clip(W+WW-WWW);
    MapCtxs[j++] = Clip(W+WW-WWW+p1-Clip(Wp1+WWp1-buf(stride*3+1)));
    MapCtxs[j++] = Clip(W+WW-WWW+p2-Clip(Wp2+WWp2-buf(stride*3+2)));
    MapCtxs[j++] = Clip(W+NEE-NE);
    MapCtxs[j++] = Clip(W+NEE-NE+p1-Clip(Wp1+buf(w-stride*2+1)-NEp1));
    MapCtxs[j++] = Clip(W+NEE-NE+p2-Clip(Wp2+buf(w-stride*2+2)-NEp2));
    MapCtxs[j++] = Clip(NN+p1-NNp1);
    MapCtxs[j++] = Clip(NN+p2-NNp2);
    MapCtxs[j++] = Clip(NN+W-NNW);
    MapCtxs[j++] = Clip(NN+W-NNW+p1-Clip(NNp1+Wp1-buf(w*2+stride+1)));
    MapCtxs[j++] = Clip(NN+W-NNW+p2-Clip(NNp2+Wp2-buf(w*2+stride+2)));
    MapCtxs[j++] = Clip(NN+NW-NNNW);
    MapCtxs[j++] = Clip(NN+NW-NNNW+p1-Clip(NNp1+NWp1-buf(w*3+stride+1)));
    MapCtxs[j++] = Clip(NN+NW-NNNW+p2-Clip(NNp2+NWp2-buf(w*3+stride+2)));
    MapCtxs[j++] = Clip(NN+NE-NNNE);
    MapCtxs[j++] = Clip(NN+NE-NNNE+p1-Clip(NNp1+NEp1-buf(w*3-stride+1)));
    MapCtxs[j++] = Clip(NN+NE-NNNE+p2-Clip(NNp2+NEp2-buf(w*3-stride+2)));
    MapCtxs[j++] = Clip(NN+NNNN-NNNNNN);
    MapCtxs[j++] = Clip(NN+NNNN-NNNNNN+p1-Clip(NNp1+buf(w*4+1)-buf(w*6+1)));
    MapCtxs[j++] = Clip(NN+NNNN-NNNNNN+p2-Clip(NNp2+buf(w*4+2)-buf(w*6+2)));
    MapCtxs[j++] = Clip(WW+p1-WWp1);
    MapCtxs[j++] = Clip(WW+p2-WWp2);
    MapCtxs[j++] = Clip(WW+WWWW-WWWWWW);
    MapCtxs[j++] = Clip(WW+WWWW-WWWWWW+p1-Clip(WWp1+buf(stride*4+1)-buf(stride*6+1)));
    MapCtxs[j++] = Clip(WW+WWWW-WWWWWW+p2-Clip(WWp2+buf(stride*4+2)-buf(stride*6+2)));
    MapCtxs[j++] = Clip(N*2-NN+p1-Clip(Np1*2-NNp1));
    MapCtxs[j++] = Clip(N*2-NN+p2-Clip(Np2*2-NNp2));
    MapCtxs[j++] = Clip(W*2-WW+p1-Clip(Wp1*2-WWp1));
    MapCtxs[j++] = Clip(W*2-WW+p2-Clip(Wp2*2-WWp2));
    MapCtxs[j++] = Clip(N*3-NN*3+NNN);
    MapCtxs[j++] = Clamp4(N*3-NN*3+NNN, W, NW, N, NE);
    MapCtxs[j++] = Clamp4(W*3-WW*3+WWW, W, NW, N, NE);
    MapCtxs[j++] = Clamp4(N*2-NN, W, NW, N, NE);
    MapCtxs[j++] = Clip((NNNNN-6*NNNN+15*NNN-20*NN+15*N+Clamp4(W*4-NWW*6+NNWWW*4-buf(w*3+4*stride), W, NW, N, NN))/6);
    MapCtxs[j++] = Clip((buf(w*3-3*stride)-4*NNEE+6*NE+Clip(W*4-NW*6+NNW*4-NNNW))/4);
    MapCtxs[j++] = Clip(((N+3*NW)/4)*3-((NNW+NNWW)/2)*3+(NNNWW*3+buf(w*3+3*stride))/4);
    MapCtxs[j++] = Clip((W*2+NW)-(WW+2*NWW)+NWWW);
    MapCtxs[j++] = (Clip(W*2-NW)+Clip(W*2-NWW)+N+NE)/4;
    MapCtxs[j++] = NNNNNN;
    MapCtxs[j++] = (NEEEE+buf(w-6*stride))/2;
    MapCtxs[j++] = (WWWWWW+WWWW)/2;
    MapCtxs[j++] = ((W+N)*3-NW*2)/4;
    MapCtxs[j++] = N;
    MapCtxs[j++] = NN;
    j = 0;
    SCMapCtxs[j++] = N+p1-Np1;
    SCMapCtxs[j++] = N+p2-Np2;
    SCMapCtxs[j++] = W+p1-Wp1;
    SCMapCtxs[j++] = W+p2-Wp2;
    SCMapCtxs[j++] = NW+p1-NWp1;
    SCMapCtxs[j++] = NW+p2-NWp2;
    SCMapCtxs[j++] = NE+p1-NEp1;
    SCMapCtxs[j++] = NE+p2-NEp2;
    SCMapCtxs[j++] = NN+p1-NNp1;
    SCMapCtxs[j++] = NN+p2-NNp2;
    SCMapCtxs[j++] = WW+p1-WWp1;
    SCMapCtxs[j++] = WW+p2-WWp2;
    SCMapCtxs[j++] = W+N-NW;
    SCMapCtxs[j++] = W+N-NW+p1-Wp1-Np1+NWp1;
    SCMapCtxs[j++] = W+N-NW+p2-Wp2-Np2+NWp2;
    SCMapCtxs[j++] = W+NE-N;
    SCMapCtxs[j++] = W+NE-N+p1-Wp1-NEp1+Np1;
    SCMapCtxs[j++] = W+NE-N+p2-Wp2-NEp2+Np2;
    SCMapCtxs[j++] = W+NEE-NE;
    SCMapCtxs[j++] = W+NEE-NE+p1-Wp1-buf(w-stride*2+1)+NEp1;
    SCMapCtxs[j++] = W+NEE-NE+p2-Wp2-buf(w-stride*2+2)+NEp2;
    SCMapCtxs[j++] = N+NN-NNN;
    SCMapCtxs[j++] = N+NN-NNN+p1-Np1-NNp1+buf(w*3+1);
    SCMapCtxs[j++] = N+NN-NNN+p2-Np2-NNp2+buf(w*3+2);
    SCMapCtxs[j++] = N+NE-NNE;
    SCMapCtxs[j++] = N+NE-NNE+p1-Np1-NEp1+buf(w*2-stride+1);
    SCMapCtxs[j++] = N+NE-NNE+p2-Np2-NEp2+buf(w*2-stride+2);
    SCMapCtxs[j++] = N+NW-NNW;
    SCMapCtxs[j++] = N+NW-NNW+p1-Np1-NWp1+buf(w*2+stride+1);
    SCMapCtxs[j++] = N+NW-NNW+p2-Np2-NWp2+buf(w*2+stride+2);
    SCMapCtxs[j++] = NE+NW-NN;
    SCMapCtxs[j++] = NE+NW-NN+p1-NEp1-NWp1+NNp1;
    SCMapCtxs[j++] = NE+NW-NN+p2-NEp2-NWp2+NNp2;
    SCMapCtxs[j++] = NW+W-NWW;
    SCMapCtxs[j++] = NW+W-NWW+p1-NWp1-Wp1+buf(w+stride*2+1);
    SCMapCtxs[j++] = NW+W-NWW+p2-NWp2-Wp2+buf(w+stride*2+2);
    SCMapCtxs[j++] = W*2-WW;
    SCMapCtxs[j++] = W*2-WW+p1-Wp1*2+WWp1;
    SCMapCtxs[j++] = W*2-WW+p2-Wp2*2+WWp2;
    SCMapCtxs[j++] = N*2-NN;
    SCMapCtxs[j++] = N*2-NN+p1-Np1*2+NNp1;
    SCMapCtxs[j++] = N*2-NN+p2-Np2*2+NNp2;
    SCMapCtxs[j++] = NW*2-NNWW;
    SCMapCtxs[j++] = NW*2-NNWW+p1-NWp1*2+buf(w*2+stride*2+1);
    SCMapCtxs[j++] = NW*2-NNWW+p2-NWp2*2+buf(w*2+stride*2+2);
    SCMapCtxs[j++] = NE*2-NNEE;
    SCMapCtxs[j++] = NE*2-NNEE+p1-NEp1*2+buf(w*2-stride*2+1);
    SCMapCtxs[j++] = NE*2-NNEE+p2-NEp2*2+buf(w*2-stride*2+2);
    SCMapCtxs[j++] = N*3-NN*3+NNN+p1-Np1*3+NNp1*3-buf(w*3+1);
    SCMapCtxs[j++] = N*3-NN*3+NNN+p2-Np2*3+NNp2*3-buf(w*3+2);
    SCMapCtxs[j++] = N*3-NN*3+NNN;
    SCMapCtxs[j++] = (W+NE*2-NNE)/2;
    SCMapCtxs[j++] = (W+NE*3-NNE*3+NNNE)/2;
    SCMapCtxs[j++] = (W+NE*2-NNE)/2+p1-(Wp1+NEp1*2-buf(w*2-stride+1))/2;
    SCMapCtxs[j++] = (W+NE*2-NNE)/2+p2-(Wp2+NEp2*2-buf(w*2-stride+2))/2;
    SCMapCtxs[j++] = NNE+NE-NNNE;
    SCMapCtxs[j++] = NNE+W-NN;
    SCMapCtxs[j++] = NNW+W-NNWW;
    j = 0;
    for (int k=(color>0)?color-1:stride-1; j<nOLS; j++) {
      pOLS[j] = Clip(int(floor(ols[j][color].Predict(ols_ctxs[j]))));
      ols[j][k].Update(p1);
    }
    int mean=W+NW+N+NE;
    const int var=(W*W+NW*NW+N*N+NE*NE-mean*mean/4)>>2;
    mean>>=2;
    const int logvar=ilog(var);

    ctx[0] = (min(color,stride-1)<<9)|((abs(W-N)>3)<<8)|((W>N)<<7)|((W>NW)<<6)|((abs(N-NW)>3)<<5)|((N>NW)<<4)|((abs(N-NE)>3)<<3)|((N>NE)<<2)|((W>WW)<<1)|(N>NN);
    ctx[1] = ((LogMeanDiffQt(buf(1),Clip(buf(w+1)+buf(w-stride+1)-buf(w*2-stride+1)))>>1)<<5)|((LogMeanDiffQt(Clip(N+NE-NNE),Clip(N+NW-NNW))>>1)<<2)|min(color,stride-1);

    U64 i=0;
    cm.set(hash(++i,        (N+1)>>1, LogMeanDiffQt(N, Clip(NN*2-NNN))));
    cm.set(hash(++i,        (W+1)>>1, LogMeanDiffQt(W, Clip(WW*2-WWW))));
    cm.set(hash(++i,        Clamp4(W+N-NW, W, NW, N, NE), LogMeanDiffQt(Clip(N+NE-NNE), Clip(N+NW-NNW))));
    cm.set(hash(++i,        (NNN+N+4)/8, Clip(N*3-NN*3+NNN)>>1));
    cm.set(hash(++i,        (WWW+W+4)/8, Clip(W*3-WW*3+WWW)>>1));
    cm.set(hash(++i, color, (W+Clip(NE*3-NNE*3+NNNE))/4, LogMeanDiffQt(N, (NW+NE)/2)));
    cm.set(hash(++i, color, Clip((-WWWW+5*WWW-10*WW+10*W+Clamp4(NE*4-NNE*6+NNNE*4-NNNNE, N, NE, NEE, NEEE))/5)/4));
    cm.set(hash(++i,        Clip(NEE+N-NNEE), LogMeanDiffQt(W, Clip(NW+NE-NNE))));
    cm.set(hash(++i,        Clip(NN+W-NNW), LogMeanDiffQt(W, Clip(NNW+WW-NNWW))));
    cm.set(hash(++i, color, p1));
    cm.set(hash(++i, color, p2));
    cm.set(hash(++i, color, Clip(W+N-NW)/2, Clip(W+p1-Wp1)/2));
    cm.set(hash(++i,        Clip(N*2-NN)/2, LogMeanDiffQt(N, Clip(NN*2-NNN))));
    cm.set(hash(++i,        Clip(W*2-WW)/2, LogMeanDiffQt(W, Clip(WW*2-WWW))));
    cm.set(hash(++i,        Clamp4(N*3-NN*3+NNN, W, NW, N, NE)/2));
    cm.set(hash(++i,        Clamp4(W*3-WW*3+WWW, W, N, NE, NEE)/2));
    cm.set(hash(++i, color, LogMeanDiffQt(W, Wp1), Clamp4((p1*W)/(Wp1<1?1:Wp1), W, N, NE, NEE))); //using max(1,Wp1) results in division by zero in VC2015
    cm.set(hash(++i, color, Clamp4(N+p2-Np2, W, NW, N, NE)));
    cm.set(hash(++i, color, Clip(W+N-NW), column[0]));
    cm.set(hash(++i, color, Clip(N*2-NN), LogMeanDiffQt(W, Clip(NW*2-NNW))));
    cm.set(hash(++i, color, Clip(W*2-WW), LogMeanDiffQt(N, Clip(NW*2-NWW))));
    cm.set(hash(++i,        (W+NEE)/2, LogMeanDiffQt(W, (WW+NE)/2)));
    cm.set(hash(++i,        (Clamp4(Clip(W*2-WW)+Clip(N*2-NN)-Clip(NW*2-NNWW), W, NW, N, NE))));
    cm.set(hash(++i, color, W, p2));
    cm.set(hash(++i,        N, NN&0x1F, NNN&0x1F));
    cm.set(hash(++i,        W, WW&0x1F, WWW&0x1F));
    cm.set(hash(++i, color, N, column[0]));
    cm.set(hash(++i, color, Clip(W+NEE-NE), LogMeanDiffQt(W, Clip(WW+NE-N))));
    cm.set(hash(++i,        NN, NNNN&0x1F, NNNNNN&0x1F, column[1]));
    cm.set(hash(++i,        WW, WWWW&0x1F, WWWWWW&0x1F, column[1]));
    cm.set(hash(++i,        NNN, NNNNNN&0x1F, buf(w*9)&0x1F, column[1]));
    cm.set(hash(++i, color, column[1]));

    cm.set(hash(++i, color, W, LogMeanDiffQt(W, WW)));
    cm.set(hash(++i, color, W, p1));
    cm.set(hash(++i, color, W/4, LogMeanDiffQt(W, p1), LogMeanDiffQt(W, p2)));
    cm.set(hash(++i, color, N, LogMeanDiffQt(N, NN)));
    cm.set(hash(++i, color, N, p1));
    cm.set(hash(++i, color, N/4, LogMeanDiffQt(N, p1), LogMeanDiffQt(N, p2)));
    cm.set(hash(++i, color, (W+N)>>3, p1>>4, p2>>4));
    cm.set(hash(++i, color, p1/2, p2/2));
    cm.set(hash(++i, color, W, p1-Wp1));
    cm.set(hash(++i, color, W+p1-Wp1));
    cm.set(hash(++i, color, N, p1-Np1));
    cm.set(hash(++i, color, N+p1-Np1));
    cm.set(hash(++i, buf(w*3-stride), buf(w*3-2*stride)));
    cm.set(hash(++i, buf(w*3+stride), buf(w*3+2*stride)));
    cm.set(hash(++i, color, mean, logvar>>4));

    i=0;
    Map[i++].set_direct((W&0xC0)|((N&0xC0)>>2)|((WW&0xC0)>>4)|(NN>>6));
    Map[i++].set_direct((N&0xC0)|((NN&0xC0)>>2)|((NE&0xC0)>>4)|(NEE>>6));
    Map[i++].set_direct(buf(1));
    Map[i++].set_direct(min(color, stride-1));
    if (Stats) {
        Stats->Image.plane = std::min<int>(color, stride-1);
        Stats->Image.pixels.W = W;
        Stats->Image.pixels.N = N;
        Stats->Image.pixels.NN = NN;
        Stats->Image.pixels.WW = WW;
        Stats->Image.pixels.Wp1 = Wp1;
        Stats->Image.pixels.Np1 = Np1;
        Stats->Image.ctx = ctx[0]>>3;
      }
  }
  U8 B=(c0<<(8-bpos));
  int i=5;

  Map[i++].set_direct((((U8)(Clip(W+N-NW)-B))*8+bpos)|(LogMeanDiffQt(Clip(N+NE-NNE), Clip(N+NW-NNW))<<11));
  Map[i++].set_direct((((U8)(Clip(N*2-NN)-B))*8+bpos)|(LogMeanDiffQt(W, Clip(NW*2-NNW))<<11));
  Map[i++].set_direct((((U8)(Clip(W*2-WW)-B))*8+bpos)|(LogMeanDiffQt(N, Clip(NW*2-NWW))<<11));
  Map[i++].set_direct((((U8)(Clip(W+N-NW)-B))*8+bpos)|(LogMeanDiffQt(p1, Clip(Wp1+Np1-NWp1))<<11));
  Map[i++].set_direct((((U8)(Clip(W+N-NW)-B))*8+bpos)|(LogMeanDiffQt(p2, Clip(Wp2+Np2-NWp2))<<11));
  Map[i++].set(hash(W-B, N-B)*8+bpos);
  Map[i++].set(hash(W-B, WW-B)*8+bpos);
  Map[i++].set(hash(N-B, NN-B)*8+bpos);
  Map[i++].set(hash(Clip(N+NE-NNE)-B, Clip(N+NW-NNW)-B)*8+bpos);
  Map[i++].set_direct((min(color, stride-1)<<11)|(((U8)(Clip(N+p1-Np1)-B))*8+bpos));
  Map[i++].set_direct((min(color, stride-1)<<11)|(((U8)(Clip(N+p2-Np2)-B))*8+bpos));
  Map[i++].set_direct((min(color, stride-1)<<11)|(((U8)(Clip(W+p1-Wp1)-B))*8+bpos));
  Map[i++].set_direct((min(color, stride-1)<<11)|(((U8)(Clip(W+p2-Wp2)-B))*8+bpos));
  
  for (int j=0; j<nMaps1; i++, j++)
    Map[i].set_direct((MapCtxs[j]-B)*8+bpos);

  for (int j=0; i<nMaps; i++, j++)
    Map[i].set_direct((pOLS[j]-B)*8+bpos);

  for (int i=0; i<nSCMaps-1; i++)
    SCMap[i].set((SCMapCtxs[i]-B)*8+bpos);
 
  // Predict next bit
  cm.mix(m);
  for (int i=0;i<nMaps;i++)
    Map[i].mix(m,1,3);
  for (int i=0;i<nSCMaps;i++)
    SCMap[i].mix(m,9,1,3);
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
  m.set(finalize64(hash(LogMeanDiffQt(W,WW,5), LogMeanDiffQt(N,NN,5), LogMeanDiffQt(W,N,5), ilog2(W), color),13), 8192);
  m.set(finalize64(hash(ctx[0], column[0]/8),13), 8192);
  m.set(finalize64(hash(LogQt(N,5), LogMeanDiffQt(N,NN,3), c0),13), 8192);
  m.set(finalize64(hash(LogQt(W,5), LogMeanDiffQt(W,WW,3), c0),13), 8192);
  m.set(min(255,(x+line)/32), 256);
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
          gray&=((B==(gray&0xFF))*0x1FF);\
      }\
  }\
}

struct TGAImage{
  U32 Header, IdLength, Bpp, ImgType, MapSize, Width, Height;
};
struct BMPImage{
  U32 Header, Offset, Bpp, Size, Palette, HdrLess, Width, Height, BitMask;
};

int imgModel(Mixer& m, ModelStats *Stats = nullptr) {
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
      case 8 : im8bitModel(m, w, Stats, gray); if (Stats) Stats->Type=(gray)?preprocessor::IMAGE8GRAY:preprocessor::IMAGE8; break;
      default: im24bitModel(m, w, Stats, alpha); if (Stats) Stats->Type=(alpha)?preprocessor::IMAGE32:preprocessor::IMAGE24;
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

inline int signedClip8(const int i) {
  return max(-128, min(127, i));
}

inline U32 SQR(U32 x) {
  return x*x;
}

void audio8bModel(Mixer& m, int info, ModelStats *Stats) {
  static const int nOLS=8, nLnrPrd=nOLS+3;
  static SmallStationaryContextMap sMap1b[nLnrPrd][3]{
    {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}},
    {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}},
    {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}}, {{11,1}, {11,1}, {11,1}}
  };
  static OLS<double, int8_t> ols[nOLS][2]{
    {{128, 24, 0.9975}, {128, 24, 0.9975}},
    {{90, 30, 0.9965}, {90, 30, 0.9965}},
    {{90, 31, 0.996}, {90, 31, 0.996}},
    {{90, 32, 0.995}, {90, 32, 0.995}},
    {{90, 33, 0.995}, {90, 33, 0.995}},
    {{90, 34, 0.9985}, {90, 34, 0.9985}},
    {{28, 4, 0.98}, {28, 4, 0.98}},
    {{28, 3, 0.992}, {28, 3, 0.992}}
  };
  static int prd[nLnrPrd][2][2]{}, residuals[nLnrPrd][2]{};
  static int stereo=0, ch=0, rpos=0, lastPos=0;
  static U32 mask=0, errLog=0, mxCtx=0;

  const int8_t B = c0<<(8-bpos);
  if (bpos==0) {
    rpos = (pos==lastPos+1)?rpos+1:0;
    lastPos = pos;
    if (rpos==0) {
      stereo = (info&1);
      mask = 0;
      Stats->Record = ((stereo+1)<<16)|(Stats->Record&0xFFFF);
      wmode=info;
    }
    ch=(stereo)?blpos&1:0;
    const int8_t s = int(((info&4)>0)?buf(1)^128:buf(1))-128;
    const int pCh = ch^stereo;
    int i = 0;
    for (errLog=0; i<nOLS; i++) {
      ols[i][pCh].Update(s);
      residuals[i][pCh] = s-prd[i][pCh][0];
      const U32 absResidual = (U32)abs(residuals[i][pCh]);
      mask+=mask+(absResidual>4);
      errLog+=SQR(absResidual);
    }
    for (; i<nLnrPrd; i++)
      residuals[i][pCh] = s-prd[i][pCh][0];
    errLog = min(0xF, ilog2(errLog));
    mxCtx = ilog2(min(0x1F, BitCount(mask)))*2+ch;

    int k1=90, k2=k1-12*stereo;
    for (int j=(i=1); j<=k1; j++, i+=1<<((j>8)+(j>16)+(j>64))) ols[1][ch].Add(X1(i));
    for (int j=(i=1); j<=k2; j++, i+=1<<((j>5)+(j>10)+(j>17)+(j>26)+(j>37))) ols[2][ch].Add(X1(i));
    for (int j=(i=1); j<=k2; j++, i+=1<<((j>3)+(j>7)+(j>14)+(j>20)+(j>33)+(j>49))) ols[3][ch].Add(X1(i));
    for (int j=(i=1); j<=k2; j++, i+=1+(j>4)+(j>8)) ols[4][ch].Add(X1(i));
    for (int j=(i=1); j<=k1; j++, i+=2+((j>3)+(j>9)+(j>19)+(j>36)+(j>61))) ols[5][ch].Add(X1(i));
    if (stereo) {
      for (i=1; i<=k1-k2; i++) {
        const double s = (double)X2(i);
        ols[2][ch].AddFloat(s);
        ols[3][ch].AddFloat(s);
        ols[4][ch].AddFloat(s);
      }
    }
    k1=28, k2=k1-6*stereo;
    for (i=1; i<=k2; i++) {
      const double s = (double)X1(i);
      ols[0][ch].AddFloat(s);
      ols[6][ch].AddFloat(s);
      ols[7][ch].AddFloat(s);
    }
    for (; i<=96; i++) ols[0][ch].Add(X1(i));
    if (stereo) {
      for (i=1; i<=k1-k2; i++) {
        const double s = (double)X2(i);
        ols[0][ch].AddFloat(s);
        ols[6][ch].AddFloat(s);
        ols[7][ch].AddFloat(s);
      }
      for (; i<=32; i++) ols[0][ch].Add(X2(i));
    }
    else
      for (; i<=128; i++) ols[0][ch].Add(X1(i));

    for (i=0; i<nOLS; i++) {
      prd[i][ch][0] = signedClip8(floor(ols[i][ch].Predict()));
      prd[i][ch][1] = signedClip8(prd[i][ch][0]+residuals[i][pCh]);
    }
    prd[i++][ch][0] = signedClip8(X1(1)*2-X1(2));
    prd[i++][ch][0] = signedClip8(X1(1)*3-X1(2)*3+X1(3));
    prd[i  ][ch][0] = signedClip8(X1(1)*4-X1(2)*6+X1(3)*4-X1(4));
    for (i=nOLS; i<nLnrPrd; i++)
      prd[i][ch][1] = signedClip8(prd[i][ch][0]+residuals[i][pCh]);
  }
  for (int i=0; i<nLnrPrd; i++) {
    const U32 ctx = (prd[i][ch][0]-B)*8+bpos;
    sMap1b[i][0].set(ctx);
    sMap1b[i][1].set(ctx);
    sMap1b[i][2].set((prd[i][ch][1]-B)*8+bpos);
    sMap1b[i][0].mix(m, 6, 1, 2+(i>=nOLS));
    sMap1b[i][1].mix(m, 9, 1, 2+(i>=nOLS));
    sMap1b[i][2].mix(m, 7, 1, 3);
  }
  m.set((errLog<<8)|c0, 4096);
  m.set((U8(mask)<<3)|(ch<<2)|(bpos>>1), 2048);
  m.set((mxCtx<<7)|(buf(1)>>1), 1280);
  m.set((errLog<<4)|(ch<<3)|bpos, 256);
  m.set(mxCtx, 10);
}

void wavModel(Mixer& m, int info, ModelStats *Stats = nullptr) {
  static int pr[3][2], n[2], counter[2];
  static double F[49][49][2],L[49][49];
  static int rpos=0, lastPos=0;
  int j,k,l,i=0;
  long double sum;
  const double a=0.996,a2=1/a;
  static SmallStationaryContextMap scm1(8,8), scm2(8,8), scm3(8,8), scm4(8,8), scm5(8,8), scm6(8,8), scm7(8,8);
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

int audioModel(Mixer& m, ModelStats *Stats = nullptr) {
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
      else if (p==(int)(16+length) && (m4(8)!=0x666d7420 || ((WAV.Chunk=i4(4)-16)&0xFFFFFFFD)!=0)){ // "fmt ", chunk size=16 or 18. should be first chunk in the file, but sometimes it's not
        length = ((i4(4)+1)&(-2)) + 8; // word aligned
        WAV.Header*=!(m4(8)==0x666d7420 && (i4(4)&0xFFFFFFFD)!=16); // was "fmt " chunk, but not 16 or 18 bytes
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

  if (info){
    if (((info-1)&2)==0)
      audio8bModel(m, info-1, Stats);
    else
      wavModel(m, info-1, Stats);
    recordModel(m, preprocessor::AUDIO, Stats);
  }

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
  static int rs1=0;  // last RS code
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
  static Array<U64> cxt(N);  // context hashes
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
    U64 n=hc*32;
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
    cxt[13]=hash(++n, rs1, prev_coef/42, prev_coef2/34, lcp[0]/60,lcp[2]/14,lcp[1]/60,lcp[3]/14);
    cxt[14]=hash(++n, mcupos&63, column>>1);
    cxt[15]=hash(++n, column>>3, min(5+2*(!comp),zu+zv), lcp[0]/10,lcp[2]/40,lcp[1]/10,lcp[3]/40);
    cxt[16]=hash(++n, ssum>>3, mcupos&63);
    cxt[17]=hash(++n, rs1, mcupos&63, run_pred[1]);
    cxt[18]=hash(++n, coef, ssum2>>5, adv_pred[3]/30, (comp)?hash(prev_coef/22,prev_coef2/50):ssum/((mcupos&0x3F)+1));
    cxt[19]=hash(++n, lcp[0]/40, lcp[1]/40, adv_pred[1]/28, (comp)?prev_coef/40+((prev_coef2/40)<<20):lcp[4]/22, min(7,zu+zv), ssum/(2*(zu+zv)+1) );
    cxt[20]=hash(++n, zv, cbuf[cpos-blockN[mcupos>>6]], adv_pred[2]/28, run_pred[2]);
    cxt[21]=hash(++n, zu, cbuf[cpos-blockW[mcupos>>6]], adv_pred[0]/28, run_pred[0]);
    cxt[22]=hash(++n, adv_pred[2]/7, run_pred[2]);
    cxt[23]=hash(  n, adv_pred[0]/7, run_pred[0]);
    cxt[24]=hash(  n, adv_pred[1]/7, run_pred[1]);
    cxt[25]=hash(++n, zv, lcp[1]/14, adv_pred[2]/16, run_pred[5]);
    cxt[26]=hash(++n, zu, lcp[0]/14, adv_pred[0]/16, run_pred[3]);
    cxt[27]=hash(++n, lcp[0]/14, lcp[1]/14, adv_pred[3]/16);
    cxt[28]=hash(++n, coef, prev_coef/10, prev_coef2/20);
    cxt[29]=hash(++n, coef, ssum>>2, prev_coef_rs);
    cxt[30]=hash(++n, coef, adv_pred[1]/17, lcp[(zu<zv)]/24,lcp[2]/20,lcp[3]/24);
    cxt[31]=hash(++n, coef, adv_pred[3]/11, lcp[(zu<zv)]/50,lcp[2+3*(zu*zv>1)]/50,lcp[3+3*(zu*zv>1)]/50);
  }

  // Predict next bit
  m1.add(128);
  jassert(hbcount<=2);
  int p;
  switch(hbcount)
  {
    case 0: for (int i=0; i<N; ++i){ cp[i]=t[cxt[i]]+1, p=sm[i].p(*cp[i]); m.add((p-2048)>>2); m1.add(p=stretch(p)); m.add(p);} break;
    case 1: { int hc=1+(huffcode&1)*3; for (int i=0; i<N; ++i){ cp[i]+=hc, p=sm[i].p(*cp[i]); m.add((p-2048)>>2); m1.add(p=stretch(p)); m.add(p); }} break;
    default: { int hc=1+(huffcode&1); for (int i=0; i<N; ++i){ cp[i]+=hc, p=sm[i].p(*cp[i]); m.add((p-2048)>>2); m1.add(p=stretch(p)); m.add(p); }} break;
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

bool exeModel(Mixer& m, bool Forced = false, ModelStats *Stats = nullptr) {
  const int N1=10, N2=10;
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
      int mask=0, count0=0, i=0;
      for (int j=0; i<N1; ++i){
        if (i>1) mask=mask*2+(buf(i-1)==0), count0+=mask&1;
        j=(i<4)?i+1:5+(i-4)*(2+(i>6));
        cm.set(hash(i, execxt(j, buf(1)*(j>6)), ((1<<N1)|mask)*(count0*N1/2>=i), (0x08|(blpos&0x07))*(i<4)));
      }

      cm.set(BrkCtx);

      mask = PrefixMask|(0xF8<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A;
      cm.set(hash(++i, OpN(Cache, 1)&(mask|RegDWordDisplacement|AddressMode), State+16*Op.BytesRead, Op.Data&mask, Op.REX, Op.Category));

      mask = 0x04|(0xFE<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A|((ModRM_mod|ModRM_reg)<<ModRMShift);
      cm.set(hash(++i,
        OpN(Cache, 1)&mask, OpN(Cache, 2)&mask, OpN(Cache, 3)&mask,
        Context+256*((Op.ModRM & ModRM_mod)==ModRM_mod),
        Op.Data&((mask|PrefixREX)^(ModRM_mod<<ModRMShift))
      ));

      mask = 0x04|CodeMask;
      cm.set(hash(++i, OpN(Cache, 1)&mask, OpN(Cache, 2)&mask, OpN(Cache, 3)&mask, OpN(Cache, 4)&mask, (Op.Data&mask)|(State<<11)|(Op.BytesRead<<15)));

      mask = 0x04|(0xFC<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A;
      cm.set(hash(++i, State+16*Op.BytesRead, Op.Data&mask, Op.Category*8 + (OpMask&0x07), Op.Flags, ((Op.SIB & SIB_base)==5)*4+((Op.ModRM & ModRM_reg)==ModRM_reg)*2+((Op.ModRM & ModRM_mod)==0)));

      mask = PrefixMask|CodeMask|OperandSizeOverride|MultiByteOpcode|PrefixREX|Prefix38|Prefix3A|HasExtraFlags|HasModRM|((ModRM_mod|ModRM_rm)<<ModRMShift);
      cm.set(hash(++i, Op.Data&mask, State+16*Op.BytesRead, Op.Flags));

      mask = PrefixMask|CodeMask|OperandSizeOverride|MultiByteOpcode|Prefix38|Prefix3A|HasExtraFlags|HasModRM;
      cm.set(hash(++i, OpN(Cache, 1)&mask, State, Op.BytesRead*2+((Op.REX&REX_w)>0), Op.Data&((U16)(mask^OperandSizeOverride))));

      mask = 0x04|(0xFE<<CodeShift)|MultiByteOpcode|Prefix38|Prefix3A|(ModRM_reg<<ModRMShift);
      cm.set(hash(++i, OpN(Cache, 1)&mask, OpN(Cache, 2)&mask, State+16*Op.BytesRead, Op.Data&(mask|PrefixMask|CodeMask)));

      cm.set(hash(++i, State+16*Op.BytesRead));

      cm.set(hash(++i,
        (0x100|B)*(Op.BytesRead>0),
        State+16*pState+256*Op.BytesRead,
        ((Op.Flags&fMODE)==fAM)*16 + (Op.REX & REX_w) + (Op.o16)*4 + ((Op.Code & 0xFE)==0xE8)*2 + ((Op.Data & MultiByteOpcode)!=0 && (Op.Code & 0xF0)==0x80)
      ));
    }
  }

  if (Valid || Forced)
    cm.mix(m);
  else{
      for (int i=0; i<(N1+N2)*7; ++i)
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
  m.set(finalize64(hash(Op.Code, State, OpN(Cache, 1)&CodeMask),13), 8192);
  m.set(finalize64(hash(State, bpos, Op.Code, Op.BytesRead),13), 8192);
  m.set(finalize64(hash(State, (bpos<<2)|(c0&3), OpCategMask&CategoryMask, ((Op.Category==OP_GEN_BRANCH)<<2)|(((Op.Flags&fMODE)==fAM)<<1)|(Op.BytesRead>0)),13), 8192);

  if (Stats)
    (*Stats).x86_64 = Valid|(Context<<1)|(s<<9);
  return Valid;
}

void indirectModel(Mixer& m) {
  static ContextMap cm(MEM(), 15);
  static U32 t1[256];
  static U16 t2[0x10000];
  static U16 t3[0x8000];
  static U16 t4[0x8000];
  static IndirectContext<U32> iCtx(16);

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
    const U8 pc=tolower(U8(c4>>8));
    iCtx+=(c=tolower(c)), iCtx=(pc<<8)|c;
    const U32 ctx0=iCtx(), mask=(U8(t1[c])==U8(t2[d]))|
                               ((U8(t1[c])==U8(t3[d2]))<<1)|
                               ((U8(t1[c])==U8(t4[d3]))<<2)|
                               ((U8(t1[c])==U8(ctx0))<<3);
    U64 i=0;
    cm.set(hash(++i,t));
    cm.set(hash(++i,t0));
    cm.set(hash(++i,ta));
    cm.set(hash(++i,tc));
    cm.set(hash(++i,t&0xff00, mask));
    cm.set(hash(++i,t0&0xff0000));
    cm.set(hash(++i,ta&0xff0000));
    cm.set(hash(++i,tc&0xff0000));
    cm.set(hash(++i,t&0xffff));
    cm.set(hash(++i,t0&0xffffff));
    cm.set(hash(++i,ta&0xffffff));
    cm.set(hash(++i,tc&0xffffff));
    cm.set(hash(++i, ctx0&0xff, c));
    cm.set(hash(++i, ctx0&0xffff));
    cm.set(hash(++i, ctx0&0x7f7fff));
  }
  cm.mix(m);
}

//////////////////////////// dmcModel //////////////////////////

// Model using DMC (Dynamic Markov Compression).
//
// The bitwise context is represented by a state graph.
//
// See the original paper: http://webhome.cs.uvic.ca/~nigelh/Publications/DMC.pdf
// See the original DMC implementation: http://maveric0.uwaterloo.ca/ftp/dmc/
//
// Main differences:
// - Instead of floats we use fixed point arithmetic.
// - The threshold for cloning a state increases gradually as memory is used up.
// - For probability estimation each state maintains both a 0,1 count ("c0" and "c1") 
//   and a bit history ("state"). The 0,1 counts are updated adaptively favoring newer events.
//   The bit history state is mapped to a probability adaptively using a StateMap.
// - The predictions of multiple "dmcModel"s are combined and stiblilized in "dmcForest". See below.

struct DMCNode { // 12 bytes
  // c0,c1: adaptive counts of zeroes and ones; 
  //        fixed point numbers with 6 integer and 10 fractional bits, i.e. scaling factor = 1024;
  //        thus the values 0 .. 65535 represent real counts of 0.0 .. 63.999
  // nx0, nx1: indexes of next DMC nodes in the state graph
  // state: bit history state - as in a contextmap
public:
  U16 c0, c1;
private:
  U32 _nx0, _nx1; // packed: their higher 28 bits are nx0, nx1; the lower 4+4 bits give the bit history state byte

public:
  U8   get_state() const         {return U8(((_nx0&0xf)<<4) | (_nx1&0xf));}
  void set_state(const U8 state) {_nx0=(_nx0&0xfffffff0) | (state>>4); _nx1=(_nx1&0xfffffff0) | (state&0xf);}
  U32  get_nx0() const           {return _nx0>>4;}
  void set_nx0(const U32 nx0)    {_nx0=(_nx0&0xf) | (nx0<<4);}
  U32  get_nx1() const           {return _nx1>>4;}
  void set_nx1(const U32 nx1)    {_nx1=(_nx1&0xf) | (nx1<<4);}
};

#define DMC_NODES_BASE (255*256) // = 65280
#define DMC_NODES_MAX  ((U64(1)<<31)/sizeof(DMCNode)) // = 178 956 970

class dmcModel {
private:
  Array<DMCNode> t;     // state graph
  StateMap32 sm;          // statemap for bit history states
  U32 top, curr;        // index of first unallocated node (i.e. number of allocated nodes); index of current node
  U32 threshold;        // cloning threshold parameter: fixed point number like c0,c1
  U32 threshold_fine;   // "threshold" scaled by 11 bits used for increasing the threshold in finer steps
  U32 extra;            // this value is used for approximating stategraph maturity level when the state graph is already full 
                        // this is the number of skipped cloning events when the counts were already large enough (>1.0)

  // helper function: adaptively increment a counter
  U32 increment_counter (const U32 x, const U32 increment) const { // x is a fixed point number as c0,c1 ; "increment"  is 0 or 1
    return (((x<<6)-x)>>6)+(increment<<10); // x * (1-1/64) + increment
  }

public: 
  dmcModel(const U64 dmc_nodes, const U32 th_start) : t(min(dmc_nodes+DMC_NODES_BASE,DMC_NODES_MAX)), sm() {resetstategraph(th_start);}

  // Initialize the state graph to a bytewise order 1 model
  // See an explanation of the initial structure in:
  // http://wing.comp.nus.edu.sg/~junping/docs/njp-icita2005.pdf
  void resetstategraph(const U32 th_start) {
    top=curr=extra=0;
    threshold=th_start;
    threshold_fine=th_start<<11;
    for (int j=0; j<256; ++j) { //256 trees
      for (int i=0; i<255; ++i) { //255 nodes in each tree
        if (i<127) { //internal tree nodes
          t[top].set_nx0(top+i+1); // left node 
          t[top].set_nx1(top+i+2); // right node
        }
        else { // 128 leaf nodes - they each references a root node of tree(i)
          int linked_tree_root=(i-127)*2*255;
          t[top].set_nx0(linked_tree_root);     // left node  -> root of tree 0,2,4... 
          t[top].set_nx1(linked_tree_root+255); // right node -> root of tree 1,3,5...
        }
        t[top].c0=t[top].c1 = th_start<1024 ? 2048 : 512; // 2.0  0.5
        t[top].set_state(0);
        top++;
      }
    }
  }

  //update stategraph
  void update() {

    U32 c0=t[curr].c0;
    U32 c1=t[curr].c1;
    const U32 n = y ==0 ? c0 : c1;

    // update counts, state
    t[curr].c0=increment_counter(c0,1-y);
    t[curr].c1=increment_counter(c1,y);
    t[curr].set_state(nex(t[curr].get_state(), y));

    // clone next state when threshold is reached
    if(n>threshold) {

      const U32 next = y==0 ? t[curr].get_nx0() : t[curr].get_nx1();
      c0=t[next].c0;
      c1=t[next].c1;
      const U32 nn=c0+c1;

      if(nn>n+threshold) {
        if(top!=t.size()) { // state graph is not yet full, let's clone
          U32 c0_top=U64(c0)*U64(n)/U64(nn);
          U32 c1_top=U64(c1)*U64(n)/U64(nn);
          c0-=c0_top;
          c1-=c1_top;

          t[top].c0=c0_top;
          t[top].c1=c1_top;
          t[next].c0=c0;
          t[next].c1=c1;

          t[top].set_nx0(t[next].get_nx0());
          t[top].set_nx1(t[next].get_nx1());
          t[top].set_state(t[next].get_state());
          if(y==0) t[curr].set_nx0(top);
          else t[curr].set_nx1(top);

          ++top;

          if(threshold<8*1024)threshold=(++threshold_fine)>>11;
        }
        else // state graph was full
          extra += nn>>10;
      }
    }

    if(y==0) curr=t[curr].get_nx0();
    else     curr=t[curr].get_nx1();
  }

  bool isfull() const {return extra>>7 > U32(t.size());}
  int pr1() const {
    const U32 n0=t[curr].c0+1;
    const U32 n1=t[curr].c1+1;
    return (n1<<12)/(n0+n1);
  }
  int pr2() {
    const U8 state=t[curr].get_state();
    return sm.p(state,256); // 64-512 are all fine
  }
  int st() {
    update();
    return stretch(pr1()) + stretch(pr2()); // average the predictions for stability
  }
};

// This class solves two problems of the DMC model
// 1) The DMC model is a memory hungry algorithm. In theory it works best when it can clone
//    nodes forever. But when the state graph is full you can't clone nodes anymore. 
//    You can either i) reset the model (the state graph) and start over
//    or ii) you can keep updating the counts forever in the already fixed state graph. Both
//    choices are troublesome: i) resetting the model degrades the predictive power significantly
//    until the graph becomes large enough again and ii) a fixed structure can't adapt anymore.
//    To solve this issue:
//    Ten models with different parameters work in tandem. Only eight of the ten models
//    are reset periodically. Due to their different cloning threshold parameters and 
//    different state graph sizes they are reset at different points in time. 
//    The remaining two models (having the highest threshold and largest stategraph) are
//    never reset and are beneficial for semi-stationary files.
// 2) The DMC model is sensitive to the cloning threshold parameter. Some files prefer
//    a smaller threshold other files prefer a larger threshold.
//    The difference in terms of compression is significant.
//    To solve this issue DMC models with different thresholds are used and their 
//    predictions are combined.
//
//    Disadvantages: with the same memory requirements we have less number of nodes
//    in each model. Also keeping more models updated at all times requires more
//    calculations and more memory access than updating one model only.
//    Advantage: more stable and better compression - even with reduced number of nodes.
//
// Further notes: 
//    Extremely small initial threshold parameters (i) help the state graph become large faster
//    and model longer input bit sequences sooner. Moreover (ii) when using a small threshold 
//    parameter the split counts c0 and c1 will be small after cloning, and after updating them
//    with 0 and 1 the prediction p=c1/(c0+c1) will be biased towards these latest events.

class dmcForest {
private:
  const U32 MODELS = 10; // 8 fast and 2 slow models
  U32 dmcparams [10] = {2,32, 64,4, 128,8, 256,16, 1024,1536};
  U64 dmcmem [10]    = {6,10, 11,7,  12,8,  13, 9,    2,   2};
  Array<dmcModel*> dmcmodels;
public:
  dmcForest():dmcmodels(MODELS) {
    for(int i=MODELS-1;i>=0;i--) 
      dmcmodels[i]=new dmcModel((MEM()>>2)/dmcmem[i],dmcparams[i]);
  }
  ~dmcForest() {
    for(int i=MODELS-1;i>=0;i--)
      delete dmcmodels[i];
  }

  // update and predict
  void mix(Mixer &m) {
    int i=MODELS;
    // the slow models predict individually
    m.add(dmcmodels[--i]->st()>>3);
    m.add(dmcmodels[--i]->st()>>3);
    // the fast models are combined for better stability
    while(i>0) {
      const int pr1=dmcmodels[--i]->st();
      const int pr2=dmcmodels[--i]->st();
      m.add((pr1+pr2)>>4);
    }

    // reset models when their structure can't adapt anymore
    // the two slow models are never reset
    if(bpos==0)
      for(int i=MODELS-3;i>=0;i--) 
        if(dmcmodels[i]->isfull())
          dmcmodels[i]->resetstategraph(dmcparams[i]);
  }
};

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

void XMLModel(Mixer& m, ModelStats *Stats = nullptr){
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
    U64 i=64;
    cm.set(hash(++i,State, (*Tag).Level, pState*2+(*Tag).EndTag, (*Tag).Name));
    cm.set(hash(++i,(*pTag).Name, State*2+(*pTag).EndTag, (*pTag).Content.Type, (*Tag).Content.Type));
    cm.set(hash(++i,State*2+(*Tag).EndTag, (*Tag).Name, (*Tag).Content.Type, c4&0xE0FF));
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

int contextModel2(ModelStats *Stats) {
  static ContextMap2 cm(MEM()*16, 10);
  static TextModel textModel(MEM()*16);
  static MatchModel matchModel(MEM()*2);
  static SparseMatchModel sparseMatchModel(MEM()/2);
  static dmcForest dmcforest;
  static RunContextMap rcm7(MEM()), rcm9(MEM()), rcm10(MEM());
  static StateMap32 StateMaps[2]={{256},{256*256}};
  static Mixer m(NUM_INPUTS, 77472, NUM_SETS, 32);
  static U32 cxt[16];
  static Filetype ft2,filetype=preprocessor::DEFAULT;
  static int size=0;  // bytes remaining in block
  static int info=0;  // image width or audio type

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
    if (Stats)
      Stats->Type = filetype;
  }

  m.update();
  m.add(64);
  
  if (bpos==0) {
    const U8 B=c4&0xFF;
    cxt[15]=(isalpha(B))?combine64(cxt[15], tolower(B)):0;
    cm.set(cxt[15]);
    for (int i=14; i>0; --i)
      cxt[i]=combine64(cxt[i-1],B);
    for (int i=0; i<7; ++i)
      cm.set(cxt[i]);
    rcm7.set(cxt[7]);
    cm.set(cxt[8]);
    rcm9.set(cxt[10]);
    rcm10.set(cxt[12]);
    cm.set(cxt[14]);
  }
  m.add((stretch(StateMaps[0].p(c0))+1)>>1);
  m.add((stretch(StateMaps[1].p(c0|(buf(1)<<8)))+1)>>1);
  int order=cm.mix(m);
  rcm7.mix(m);
  rcm9.mix(m);
  rcm10.mix(m);

  int ismatch=ilog(matchModel.Predict(m, buf, Stats));
  if (filetype==preprocessor::IMAGE1) return im1bitModel(m, info), m.p();
  if (filetype==preprocessor::IMAGE4) return im4bitModel(m, info), m.p();
  if (filetype==preprocessor::IMAGE8) return im8bitModel(m, info, Stats), m.p();
  if (filetype==preprocessor::IMAGE8GRAY) return im8bitModel(m, info, Stats, 1), m.p();
  if (filetype==preprocessor::IMAGE24) return im24bitModel(m, info, Stats), m.p();
  if (filetype==preprocessor::IMAGE32) return im24bitModel(m, info, Stats, 1), m.p();
  if ((filetype!=preprocessor::EXE && jpegModel(m)) || (size>0 && imgModel(m, Stats)) || audioModel(m, Stats))
    return m.p();

  sparseMatchModel.Predict(m, buf, Stats);
  sparseModel(m,ismatch,order);
  sparseModel1(m,ismatch,order);
  distanceModel(m);
  picModel(m);
  recordModel(m, filetype, Stats);
  recordModel1(m);
  wordModel(m);
  nestModel(m);
  indirectModel(m);
  dmcforest.mix(m);
  XMLModel(m, Stats);
  textModel.Predict(m, buf, Stats);
  exeModel(m, true, Stats);
  linearPredictionModel(m);
  
  m.set((max(0, order-3)<<3)|bpos, 64);
  order=max(0,order-5);

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
  struct {
    APM APMs[4];
    APM1 APM1s[3];
  } Text;
  struct {
    struct {
      APM APMs[4];
      APM1 APM1s[2];
    } Color, Palette;
    struct {
      APM APMs[3];
    } Gray;
  } Image;
  struct {
    APM1 APM1s[7];
  } Generic;
  ModelStats stats;
public:
  Predictor();
  int p() const {return pr;}
  void update();
};

Predictor::Predictor():
  pr(2048),
  Text{{{0x10000}, {0x10000}, {0x10000}, {0x10000}}, {{0x10000}, {0x10000}, {0x10000}}},
  Image{ 
    {{{0x1000}, {0x10000}, {0x10000}, {0x10000}}, {{0x10000}, {0x10000}}}, // color
    {{{0x1000}, {0x10000}, {0x10000}, {0x10000}}, {{0x10000}, {0x10000}}}, // palette
    {{{0x1000}, {0x10000}, {0x10000}}} //gray
  },
  Generic{{{0x2000}, {0x10000}, {0x10000}, {0x10000}, {0x10000}, {0x10000}, {0x10000}}}
  {
  memset(&stats, 0, sizeof(ModelStats));
  for (int i=0; i<1024; ++i)
    dt[i]=16384/(i+i+3);
}

void Predictor::update() {
  static U32 x5=0;

  c0+=c0+y;
  stats.Misses+=stats.Misses+((pr>>11)!=y);
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
  grp0 = (bpos>0)?AsciiGroupC0[(1<<bpos)-2+(c0&((1<<bpos)-1))]:0;
  
  int pr0=contextModel2(&stats);
  AddPrediction(pr0);
  
  int pr1, pr2, pr3;
  switch (stats.Type) {
    case preprocessor::TEXT: {
      int limit=0x3FF>>((blpos<0xFFF)*2);
      pr  = Text.APMs[0].p(pr0, (c0<<8)|(stats.Text.mask&0xF)|((stats.Misses&0xF)<<4), limit); AddPrediction(pr);
      pr1 = Text.APMs[1].p(pr0, finalize64(hash(bpos, stats.Misses&3, c4&0xffff, stats.Text.mask>>4),16), limit); AddPrediction(pr1);
      pr2 = Text.APMs[2].p(pr0, finalize64(hash(c0, stats.Match.expectedByte, std::min<U32>(3, ilog2(stats.Match.length+1))),16), limit); AddPrediction(pr2);
      pr3 = Text.APMs[3].p(pr0, finalize64(hash(c0, c4&0xffff, stats.Text.firstLetter),16), limit); AddPrediction(pr3);

      pr0 = (pr0+pr1+pr2+pr3+2)>>2; AddPrediction(pr0);

      pr1 = Text.APM1s[0].p(pr0, finalize64(hash(stats.Match.expectedByte, std::min<U32>(3, ilog2(stats.Match.length+1)), c4&0xff),16)); AddPrediction(pr1);
      pr2 = Text.APM1s[1].p(pr, finalize64(hash(c0, c4&0x00ffffff),16), 6); AddPrediction(pr2);
      pr3 = Text.APM1s[2].p(pr, finalize64(hash(c0, c4&0xffffff00),16), 6); AddPrediction(pr3);
      
      pr = (pr+pr1+pr2+pr3+2)>>2; AddPrediction(pr);
      pr = (pr+pr0+1)>>1; AddPrediction(pr);
      break;
    }
    case preprocessor::IMAGE24: case preprocessor::IMAGE32: {
      int limit=0x3FF>>((blpos<0xFFF)*4);
      pr  = Image.Color.APMs[0].p(pr0, (c0<<4)|(stats.Misses&0xF), limit); AddPrediction(pr);
      pr1 = Image.Color.APMs[1].p(pr0, finalize64(hash(c0, stats.Image.pixels.W, stats.Image.pixels.WW),16), limit); AddPrediction(pr1);
      pr2 = Image.Color.APMs[2].p(pr0, finalize64(hash(c0, stats.Image.pixels.N, stats.Image.pixels.NN),16), limit); AddPrediction(pr2);
      pr3 = Image.Color.APMs[3].p(pr0, (c0<<8)|stats.Image.ctx, limit); AddPrediction(pr3);

      pr0 = (pr0+pr1+pr2+pr3+2)>>2; AddPrediction(pr0);

      pr1 = Image.Color.APM1s[0].p(pr, finalize64(hash(c0, stats.Image.pixels.W, (c4&0xff)-stats.Image.pixels.Wp1, stats.Image.plane),16)); AddPrediction(pr1);
      pr2 = Image.Color.APM1s[1].p(pr, finalize64(hash(c0, stats.Image.pixels.N, (c4&0xff)-stats.Image.pixels.Np1, stats.Image.plane),16)); AddPrediction(pr2);

      pr=(pr*2+pr1*3+pr2*3+4)>>3; AddPrediction(pr);
      pr = (pr+pr0+1)>>1; AddPrediction(pr);
      break;
    }
    case preprocessor::IMAGE8GRAY: {
      int limit=0x3FF>>((blpos<0xFFF)*4);
      pr  = Image.Gray.APMs[0].p(pr0, (c0<<4)|(stats.Misses&0xF), limit); AddPrediction(pr);
      pr1 = Image.Gray.APMs[1].p(pr, (c0<<8)|stats.Image.ctx, limit); AddPrediction(pr1);
      pr2 = Image.Gray.APMs[2].p(pr0, bpos|(stats.Image.ctx&0xF8)|(stats.Match.expectedByte<<8), limit); AddPrediction(pr2);

      pr0 = (2*pr0+pr1+pr2+2)>>2; AddPrediction(pr0);
      pr = (pr+pr0+1)>>1; AddPrediction(pr);
      break;
    }
    case preprocessor::IMAGE8: {
      int limit=0x3FF>>((blpos<0xFFF)*4);
      pr  = Image.Palette.APMs[0].p(pr0, (c0<<4)|(stats.Misses&0xF), limit); AddPrediction(pr);
      pr1 = Image.Palette.APMs[1].p(pr0, finalize64(hash(c0, stats.Image.pixels.W, stats.Image.pixels.N),16), limit); AddPrediction(pr1);
      pr2 = Image.Palette.APMs[2].p(pr0, finalize64(hash(c0, stats.Image.pixels.N, stats.Image.pixels.NN),16), limit); AddPrediction(pr2);
      pr3 = Image.Palette.APMs[3].p(pr0, finalize64(hash(c0, stats.Image.pixels.W, stats.Image.pixels.WW),16), limit); AddPrediction(pr3);

      pr0 = (pr0+pr1+pr2+pr3+2)>>2; AddPrediction(pr0);
      
      pr1 = Image.Palette.APM1s[0].p(pr0, finalize64(hash(c0, stats.Match.expectedByte, stats.Image.pixels.N),16), 5); AddPrediction(pr1);
      pr2 = Image.Palette.APM1s[1].p(pr , finalize64(hash(c0, stats.Image.pixels.W, stats.Image.pixels.N),16), 6); AddPrediction(pr2);

      pr = (pr*2+pr1+pr2+2)>>2; AddPrediction(pr);
      pr = (pr+pr0+1)>>1; AddPrediction(pr);
      break;
    }
    default: {
      pr  = Generic.APM1s[0].p(pr0, (std::min<U32>(3, ilog2(stats.Match.length+1))<<11)|(c0<<3)|(stats.Misses&0x7)); AddPrediction(pr);
      const U16 ctx1 = c0 | buf(1)<<8;
      const U16 ctx2 = c0 ^ finalize64(hash(c4&0xffff),16);
      const U16 ctx3 = c0 ^ finalize64(hash(c4&0xffffff),16);
      pr1 = Generic.APM1s[1].p(pr0, ctx1); AddPrediction(pr1);
      pr2 = Generic.APM1s[2].p(pr0, ctx2); AddPrediction(pr2);
      pr3 = Generic.APM1s[3].p(pr0, ctx3); AddPrediction(pr3);

      pr0 = (pr0+pr1+pr2+pr3+2)>>2;

      pr1 = Generic.APM1s[4].p(pr, (stats.Match.expectedByte<<8)|buf(1)); AddPrediction(pr1);
      pr2 = Generic.APM1s[5].p(pr, ctx2); AddPrediction(pr2);
      pr3 = Generic.APM1s[6].p(pr, ctx3); AddPrediction(pr3);

      pr = (pr+pr1+pr2+pr3+2)>>2; AddPrediction(pr);
      pr = (pr+pr0+1)>>1; AddPrediction(pr);
    }
  }
  
  ResetPredictions();
  last_prediction = pr;
}

}

PAQ8::PAQ8(int memory) {
  paq8::level = memory;
  paq8::buf.setsize(paq8::MEM()*8);
  predictor_.reset(new paq8::Predictor());
}

const std::valarray<float>& PAQ8::Predict() {
  return paq8::model_predictions;
}

unsigned int PAQ8::NumOutputs() {
  return paq8::model_predictions.size();
}

void PAQ8::Perceive(int bit) {
  paq8::y = bit;
  predictor_->update();
}
