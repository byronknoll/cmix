// This code is a hybrid of paq8l and paq8pxd_8 (released by Kaido Orav).

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
#include <assert.h>

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

// min, max functions
#ifndef WINDOWS
inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}
#endif

// Error handler: print message if any, and exit
void quit(const char* message=0) {
  throw message;
}

// strings are equal ignoring case?
int equals(const char* a, const char* b) {
  assert(a && b);
  while (*a && *b) {
    int c1=*a;
    if (c1>='A'&&c1<='Z') c1+='a'-'A';
    int c2=*b;
    if (c2>='A'&&c2<='Z') c2+='a'-'A';
    if (c1!=c2) return 0;
    ++a;
    ++b;
  }
  return *a==*b;
}

//////////////////////////// Array ////////////////////////////

// Array<T, ALIGN> a(n); creates n elements of T initialized to 0 bits.
// Constructors for T are not called.
// Indexing is bounds checked if assertions are on.
// a.size() returns n.
// a.resize(n) changes size to n, padding with 0 bits or truncating.
// a.push_back(x) appends x and increases size by 1, reserving up to size*2.
// a.pop_back() decreases size by 1, does not free memory.
// Copy and assignment are not supported.
// Memory is aligned on a ALIGN byte boundary (power of 2), default is none.

template <class T, int ALIGN=0> class Array {
private:
  int n;     // user size
  int reserved;  // actual size
  char *ptr; // allocated memory, zeroed
  T* data;   // start of n elements of aligned data
  void create(int i);  // create with size i
public:
  explicit Array(int i=0) {create(i);}
  ~Array();
  T& operator[](int i) {
    return data[i];
  }
  const T& operator[](int i) const {
    return data[i];
  }
  int size() const {return n;}
  void resize(int i);  // change size to i
  void pop_back() {if (n>0) --n;}  // decrement size
  void push_back(const T& x);  // increment size, append x
private:
  Array(const Array&);  // no copy or assignment
  Array& operator=(const Array&);
};

template<class T, int ALIGN> void Array<T, ALIGN>::resize(int i) {
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

template<class T, int ALIGN> void Array<T, ALIGN>::create(int i) {
  n=reserved=i;
  if (i<=0) {
    data=0;
    ptr=0;
    return;
  }
  const int sz=ALIGN+n*sizeof(T);
  ptr = (char*)calloc(sz, 1);
  if (!ptr) quit("Out of memory");
  data = (ALIGN ? (T*)(ptr+ALIGN-(((long)ptr)&(ALIGN-1))) : (T*)ptr);
  assert((char*)data>=ptr && (char*)data<=ptr+ALIGN);
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

// A tiny subset of std::string
// size() includes NUL terminator.

class String: public Array<char> {
public:
  const char* c_str() const {return &(*this)[0];}
  void operator=(const char* s) {
    resize(strlen(s)+1);
    strcpy(&(*this)[0], s);
  }
  void operator+=(const char* s) {
    assert(s);
    pop_back();
    while (*s) push_back(*s++);
    push_back(0);
  }
  String(const char* s=""): Array<char>(1) {
    (*this)+=s;
  }
};


//////////////////////////// rnd ///////////////////////////////

// 32-bit pseudo random number generator
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

// Buf(n) buf; creates an array of n bytes (must be a power of 2).
// buf[i] returns a reference to the i'th byte with wrap (no out of bounds).
// buf(i) returns i'th byte back from pos (i > 0) 
// buf.size() returns n.

int pos;  // Number of input bytes in buf (not wrapped)

class Buf {
  Array<U8> b;
public:
  Buf(int i=0): b(i) {}
  void setsize(int i) {
    if (!i) return;
    assert(i>0 && (i&(i-1))==0);
    b.resize(i);
  }
  U8& operator[](int i) {
    return b[i&(b.size()-1)];
  }
  int operator()(int i) const {
    assert(i>0);
    return b[(pos-i)&(b.size()-1)];
  }
  int size() const {
    return b.size();
  }
};

/////////////////////// Global context /////////////////////////

int level=DEFAULT_OPTION;  // Compression level 0 to 9
#define MEM (0x10000<<level)
int y=0;  // Last bit, 0 or 1, set by encoder

// Global context set by Predictor and available to all models.
int c0=1; // Last 0-7 bits of the partial byte with a leading 1 bit (1-255)
U32 c4=0; // Last 4 whole bytes, packed.  Last byte is bits 0-7.
int bpos=0; // bits in c0 (0 to 7)
Buf buf;  // Rotating input queue set by Predictor

///////////////////////////// ilog //////////////////////////////

// ilog(x) = round(log2(x) * 16), 0 <= x < 64K
class Ilog {
  Array<U8> t;
public:
  int operator()(U16 x) const {return t[x];}
  Ilog();
} ilog;

// Compute lookup table by numerical integration of 1/x
Ilog::Ilog(): t(65536) {
  U32 x=14155776;
  for (int i=2; i<65536; ++i) {
    x+=774541002/(i*2-1);  // numerator is 2^29/ln 2
    t[i]=x>>24;
  }
}

// llog(x) accepts 32 bits
inline int llog(U32 x) {
  if (x>=0x1000000)
    return 256+ilog(x>>16);
  else if (x>=0x10000)
    return 128+ilog(x>>8);
  else
    return ilog(x);
}

///////////////////////// state table ////////////////////////

// State table:
//   nex(state, 0) = next state if bit y is 0, 0 <= state < 256
//   nex(state, 1) = next state if bit y is 1
//   nex(state, 2) = number of zeros in bit history represented by state
//   nex(state, 3) = number of ones represented
//
// States represent a bit history within some context.
// State 0 is the starting state (no bits seen).
// States 1-30 represent all possible sequences of 1-4 bits.
// States 31-252 represent a pair of counts, (n0,n1), the number
//   of 0 and 1 bits respectively.  If n0+n1 < 16 then there are
//   two states for each pair, depending on if a 0 or 1 was the last
//   bit seen.
// If n0 and n1 are too large, then there is no state to represent this
// pair, so another state with about the same ratio of n0/n1 is substituted.
// Also, when a bit is observed and the count of the opposite bit is large,
// then part of this count is discarded to favor newer data over old.

static const U8 State_table[256][4]={
  {  1,  2, 0, 0},{  3,  5, 1, 0},{  4,  6, 0, 1},{  7, 10, 2, 0}, // 0-3
  {  8, 12, 1, 1},{  9, 13, 1, 1},{ 11, 14, 0, 2},{ 15, 19, 3, 0}, // 4-7
  { 16, 23, 2, 1},{ 17, 24, 2, 1},{ 18, 25, 2, 1},{ 20, 27, 1, 2}, // 8-11
  { 21, 28, 1, 2},{ 22, 29, 1, 2},{ 26, 30, 0, 3},{ 31, 33, 4, 0}, // 12-15
  { 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1}, // 16-19
  { 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2}, // 20-23
  { 34, 37, 2, 2},{ 34, 37, 2, 2},{ 36, 39, 1, 3},{ 36, 39, 1, 3}, // 24-27
  { 36, 39, 1, 3},{ 36, 39, 1, 3},{ 38, 40, 0, 4},{ 41, 43, 5, 0}, // 28-31
  { 42, 45, 4, 1},{ 42, 45, 4, 1},{ 44, 47, 3, 2},{ 44, 47, 3, 2}, // 32-35
  { 46, 49, 2, 3},{ 46, 49, 2, 3},{ 48, 51, 1, 4},{ 48, 51, 1, 4}, // 36-39
  { 50, 52, 0, 5},{ 53, 43, 6, 0},{ 54, 57, 5, 1},{ 54, 57, 5, 1}, // 40-43
  { 56, 59, 4, 2},{ 56, 59, 4, 2},{ 58, 61, 3, 3},{ 58, 61, 3, 3}, // 44-47
  { 60, 63, 2, 4},{ 60, 63, 2, 4},{ 62, 65, 1, 5},{ 62, 65, 1, 5}, // 48-51
  { 50, 66, 0, 6},{ 67, 55, 7, 0},{ 68, 57, 6, 1},{ 68, 57, 6, 1}, // 52-55
  { 70, 73, 5, 2},{ 70, 73, 5, 2},{ 72, 75, 4, 3},{ 72, 75, 4, 3}, // 56-59
  { 74, 77, 3, 4},{ 74, 77, 3, 4},{ 76, 79, 2, 5},{ 76, 79, 2, 5}, // 60-63
  { 62, 81, 1, 6},{ 62, 81, 1, 6},{ 64, 82, 0, 7},{ 83, 69, 8, 0}, // 64-67
  { 84, 71, 7, 1},{ 84, 71, 7, 1},{ 86, 73, 6, 2},{ 86, 73, 6, 2}, // 68-71
  { 44, 59, 5, 3},{ 44, 59, 5, 3},{ 58, 61, 4, 4},{ 58, 61, 4, 4}, // 72-75
  { 60, 49, 3, 5},{ 60, 49, 3, 5},{ 76, 89, 2, 6},{ 76, 89, 2, 6}, // 76-79
  { 78, 91, 1, 7},{ 78, 91, 1, 7},{ 80, 92, 0, 8},{ 93, 69, 9, 0}, // 80-83
  { 94, 87, 8, 1},{ 94, 87, 8, 1},{ 96, 45, 7, 2},{ 96, 45, 7, 2}, // 84-87
  { 48, 99, 2, 7},{ 48, 99, 2, 7},{ 88,101, 1, 8},{ 88,101, 1, 8}, // 88-91
  { 80,102, 0, 9},{103, 69,10, 0},{104, 87, 9, 1},{104, 87, 9, 1}, // 92-95
  {106, 57, 8, 2},{106, 57, 8, 2},{ 62,109, 2, 8},{ 62,109, 2, 8}, // 96-99
  { 88,111, 1, 9},{ 88,111, 1, 9},{ 80,112, 0,10},{113, 85,11, 0}, // 100-103
  {114, 87,10, 1},{114, 87,10, 1},{116, 57, 9, 2},{116, 57, 9, 2}, // 104-107
  { 62,119, 2, 9},{ 62,119, 2, 9},{ 88,121, 1,10},{ 88,121, 1,10}, // 108-111
  { 90,122, 0,11},{123, 85,12, 0},{124, 97,11, 1},{124, 97,11, 1}, // 112-115
  {126, 57,10, 2},{126, 57,10, 2},{ 62,129, 2,10},{ 62,129, 2,10}, // 116-119
  { 98,131, 1,11},{ 98,131, 1,11},{ 90,132, 0,12},{133, 85,13, 0}, // 120-123
  {134, 97,12, 1},{134, 97,12, 1},{136, 57,11, 2},{136, 57,11, 2}, // 124-127
  { 62,139, 2,11},{ 62,139, 2,11},{ 98,141, 1,12},{ 98,141, 1,12}, // 128-131
  { 90,142, 0,13},{143, 95,14, 0},{144, 97,13, 1},{144, 97,13, 1}, // 132-135
  { 68, 57,12, 2},{ 68, 57,12, 2},{ 62, 81, 2,12},{ 62, 81, 2,12}, // 136-139
  { 98,147, 1,13},{ 98,147, 1,13},{100,148, 0,14},{149, 95,15, 0}, // 140-143
  {150,107,14, 1},{150,107,14, 1},{108,151, 1,14},{108,151, 1,14}, // 144-147
  {100,152, 0,15},{153, 95,16, 0},{154,107,15, 1},{108,155, 1,15}, // 148-151
  {100,156, 0,16},{157, 95,17, 0},{158,107,16, 1},{108,159, 1,16}, // 152-155
  {100,160, 0,17},{161,105,18, 0},{162,107,17, 1},{108,163, 1,17}, // 156-159
  {110,164, 0,18},{165,105,19, 0},{166,117,18, 1},{118,167, 1,18}, // 160-163
  {110,168, 0,19},{169,105,20, 0},{170,117,19, 1},{118,171, 1,19}, // 164-167
  {110,172, 0,20},{173,105,21, 0},{174,117,20, 1},{118,175, 1,20}, // 168-171
  {110,176, 0,21},{177,105,22, 0},{178,117,21, 1},{118,179, 1,21}, // 172-175
  {110,180, 0,22},{181,115,23, 0},{182,117,22, 1},{118,183, 1,22}, // 176-179
  {120,184, 0,23},{185,115,24, 0},{186,127,23, 1},{128,187, 1,23}, // 180-183
  {120,188, 0,24},{189,115,25, 0},{190,127,24, 1},{128,191, 1,24}, // 184-187
  {120,192, 0,25},{193,115,26, 0},{194,127,25, 1},{128,195, 1,25}, // 188-191
  {120,196, 0,26},{197,115,27, 0},{198,127,26, 1},{128,199, 1,26}, // 192-195
  {120,200, 0,27},{201,115,28, 0},{202,127,27, 1},{128,203, 1,27}, // 196-199
  {120,204, 0,28},{205,115,29, 0},{206,127,28, 1},{128,207, 1,28}, // 200-203
  {120,208, 0,29},{209,125,30, 0},{210,127,29, 1},{128,211, 1,29}, // 204-207
  {130,212, 0,30},{213,125,31, 0},{214,137,30, 1},{138,215, 1,30}, // 208-211
  {130,216, 0,31},{217,125,32, 0},{218,137,31, 1},{138,219, 1,31}, // 212-215
  {130,220, 0,32},{221,125,33, 0},{222,137,32, 1},{138,223, 1,32}, // 216-219
  {130,224, 0,33},{225,125,34, 0},{226,137,33, 1},{138,227, 1,33}, // 220-223
  {130,228, 0,34},{229,125,35, 0},{230,137,34, 1},{138,231, 1,34}, // 224-227
  {130,232, 0,35},{233,125,36, 0},{234,137,35, 1},{138,235, 1,35}, // 228-231
  {130,236, 0,36},{237,125,37, 0},{238,137,36, 1},{138,239, 1,36}, // 232-235
  {130,240, 0,37},{241,125,38, 0},{242,137,37, 1},{138,243, 1,37}, // 236-239
  {130,244, 0,38},{245,135,39, 0},{246,137,38, 1},{138,247, 1,38}, // 240-243
  {140,248, 0,39},{249,135,40, 0},{250, 69,39, 1},{ 80,251, 1,39}, // 244-247
  {140,252, 0,40},{249,135,41, 0},{250, 69,40, 1},{ 80,251, 1,40}, // 248-251
  {140,252, 0,41}};  // 252, 253-255 are reserved

#define nex(state,sel) State_table[state][sel]

///////////////////////////// Squash //////////////////////////////

// return p = 1/(1 + exp(-d)), d scaled by 8 bits, p scaled by 12 bits
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

// Inverse of squash. d = ln(p/(1-p)), d scaled by 8 bits, p by 12 bits.
// d has range -2047 to 2047 representing -8 to 8.  p has range 0 to 4095.

class Stretch {
  Array<short> t;
public:
  Stretch();
  int operator()(int p) const {
    assert(p>=0 && p<4096);
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

// Mixer m(N, M, S=1, w=0) combines models using M neural networks with
//   N inputs each, of which up to S may be selected.  If S > 1 then
//   the outputs of these neural networks are combined using another
//   neural network (with parameters S, 1, 1).  If S = 1 then the
//   output is direct.  The weights are initially w (+-32K).
//   It is used as follows:
// m.update() trains the network where the expected output is the
//   last bit (in the global variable y).
// m.add(stretch(p)) inputs prediction from one of N models.  The
//   prediction should be positive to predict a 1 bit, negative for 0,
//   nominally +-256 to +-2K.  The maximum allowed value is +-32K but
//   using such large values may cause overflow if N is large.
// m.set(cxt, range) selects cxt as one of 'range' neural networks to
//   use.  0 <= cxt < range.  Should be called up to S times such
//   that the total of the ranges is <= M.
// m.p() returns the output prediction that the next bit is 1 as a
//   12 bit number (0 to 4095).

// dot_product returns dot product t*w of n elements.  n is rounded
// up to a multiple of 8.  Result is scaled down by 8 bits.

#if !defined(__GNUC__)

#if (2 == _M_IX86_FP) // 2 if /arch:SSE2 was used.
# define __SSE2__
#elif (1 == _M_IX86_FP) // 1 if /arch:SSE was used.
# define __SSE__
#endif

#endif /* __GNUC__ */

/**
 * Vector product a*b of n signed words, returning signed integer scaled down by 8 bits.
 * n is rounded up to a multiple of 8.
 */
static int dot_product (const short* const t, const short* const w, int n);

/**
 * Train n neural network weights w[n] on inputs t[n] and err.
 * w[i] += ((t[i]*2*err)+(1<<16))>>17 bounded to +- 32K.
 * n is rounded up to a multiple of 8.
 */
static void train (const short* const t, short* const w, int n, const int e);

#if defined(__SSE2__) // faster
#include <emmintrin.h>
#define OPTIMIZE "SSE2-"

static int dot_product (const short* const t, const short* const w, int n) {
  assert(n == ((n + 7) & -8));
  __m128i sum = _mm_setzero_si128 ();
  while ((n -= 8) >= 0) { // Each loop sums eight products
    __m128i tmp = _mm_madd_epi16 (*(__m128i *) &t[n], *(__m128i *) &w[n]); // t[n] * w[n] + t[n+1] * w[n+1]
    tmp = _mm_srai_epi32 (tmp, 8); //                                        (t[n] * w[n] + t[n+1] * w[n+1]) >> 8
    sum = _mm_add_epi32 (sum, tmp); //                                sum += (t[n] * w[n] + t[n+1] * w[n+1]) >> 8
  }
  sum = _mm_add_epi32 (sum, _mm_srli_si128 (sum, 8)); // Add eight sums together ...
  sum = _mm_add_epi32 (sum, _mm_srli_si128 (sum, 4));
  return _mm_cvtsi128_si32 (sum); //                     ...  and scale back to integer
}

static void train (const short* const t, short* const w, int n, const int e) {
  assert(n == ((n + 7) & -8));
  if (e) {
    const __m128i one = _mm_set1_epi16 (1);
    const __m128i err = _mm_set1_epi16 (short(e));
    while ((n -= 8) >= 0) { // Each iteration adjusts eight weights
      __m128i tmp = _mm_adds_epi16 (*(__m128i *) &t[n], *(__m128i *) &t[n]); // t[n] * 2
      tmp = _mm_mulhi_epi16 (tmp, err); //                                     (t[n] * 2 * err) >> 16
      tmp = _mm_adds_epi16 (tmp, one); //                                     ((t[n] * 2 * err) >> 16) + 1
      tmp = _mm_srai_epi16 (tmp, 1); //                                      (((t[n] * 2 * err) >> 16) + 1) >> 1
      tmp = _mm_adds_epi16 (tmp, *(__m128i *) &w[n]); //                    ((((t[n] * 2 * err) >> 16) + 1) >> 1) + w[n]
      *(__m128i *) &w[n] = tmp; //                                          save the new eight weights, bounded to +- 32K
    }
  }
}

#elif defined(__SSE__) // fast
#include <xmmintrin.h>
#define OPTIMIZE "SSE-"

static sint32 dot_product (const sint16* const t, const sint16* const w, sint32 n) {
  assert(n == ((n + 7) & -8));
  __m64 sum = _mm_setzero_si64 ();
  while ((n -= 8) >= 0) { // Each loop sums eight products
    __m64 tmp = _mm_madd_pi16 (*(__m64 *) &t[n], *(__m64 *) &w[n]); //   t[n] * w[n] + t[n+1] * w[n+1]
    tmp = _mm_srai_pi32 (tmp, 8); //                                    (t[n] * w[n] + t[n+1] * w[n+1]) >> 8
    sum = _mm_add_pi32 (sum, tmp); //                            sum += (t[n] * w[n] + t[n+1] * w[n+1]) >> 8

    tmp = _mm_madd_pi16 (*(__m64 *) &t[n + 4], *(__m64 *) &w[n + 4]); // t[n+4] * w[n+4] + t[n+5] * w[n+5]
    tmp = _mm_srai_pi32 (tmp, 8); //                                    (t[n+4] * w[n+4] + t[n+5] * w[n+5]) >> 8
    sum = _mm_add_pi32 (sum, tmp); //                            sum += (t[n+4] * w[n+4] + t[n+5] * w[n+5]) >> 8
  }
  sum = _mm_add_pi32 (sum, _mm_srli_si64 (sum, 32)); // Add eight sums together ...
  const sint32 retval = _mm_cvtsi64_si32 (sum); //                     ...  and scale back to integer
  _mm_empty(); // Empty the multimedia state
  return retval;
}

static void train (const sint16* const t, sint16* const w, sint32 n, const sint32 e) {
  assert(n == ((n + 7) & -8));
  if (e) {
    const __m64 one = _mm_set1_pi16 (1);
    const __m64 err = _mm_set1_pi16 (sint16(e));
    while ((n -= 8) >= 0) { // Each iteration adjusts eight weights
      __m64 tmp = _mm_adds_pi16 (*(__m64 *) &t[n], *(__m64 *) &t[n]); //   t[n] * 2
      tmp = _mm_mulhi_pi16 (tmp, err); //                                 (t[n] * 2 * err) >> 16
      tmp = _mm_adds_pi16 (tmp, one); //                                 ((t[n] * 2 * err) >> 16) + 1
      tmp = _mm_srai_pi16 (tmp, 1); //                                  (((t[n] * 2 * err) >> 16) + 1) >> 1
      tmp = _mm_adds_pi16 (tmp, *(__m64 *) &w[n]); //                  ((((t[n] * 2 * err) >> 16) + 1) >> 1) + w[n]
      *(__m64 *) &w[n] = tmp; //                                       save the new four weights, bounded to +- 32K

      tmp = _mm_adds_pi16 (*(__m64 *) &t[n + 4], *(__m64 *) &t[n + 4]); // t[n+4] * 2
      tmp = _mm_mulhi_pi16 (tmp, err); //                                 (t[n+4] * 2 * err) >> 16
      tmp = _mm_adds_pi16 (tmp, one); //                                 ((t[n+4] * 2 * err) >> 16) + 1
      tmp = _mm_srai_pi16 (tmp, 1); //                                  (((t[n+4] * 2 * err) >> 16) + 1) >> 1
      tmp = _mm_adds_pi16 (tmp, *(__m64 *) &w[n + 4]); //              ((((t[n+4] * 2 * err) >> 16) + 1) >> 1) + w[n]
      *(__m64 *) &w[n + 4] = tmp; //                                   save the new four weights, bounded to +- 32K
    }
    _mm_empty(); // Empty the multimedia state
  }
}
#else
// dot_product returns dot product t*w of n elements.  n is rounded
// up to a multiple of 8.  Result is scaled down by 8 bits.
int dot_product(short *t, short *w, int n) {
  int sum=0;
  n=(n+7)&-8;
  for (int i=0; i<n; i+=2)
    sum+=(t[i]*w[i]+t[i+1]*w[i+1]) >> 8;
  return sum;
}

// Train neural network weights w[n] given inputs t[n] and err.
// w[i] += t[i]*err, i=0..n-1.  t, w, err are signed 16 bits (+- 32K).
// err is scaled 16 bits (representing +- 1/2).  w[i] is clamped to +- 32K
// and rounded.  n is rounded up to a multiple of 8.

void train(short *t, short *w, int n, int err) {
  n=(n+7)&-8;
  for (int i=0; i<n; ++i) {
    int wt=w[i]+(((t[i]*err*2>>16)+1)>>1);
    if (wt<-32768) wt=-32768;
    if (wt>32767) wt=32767;
    w[i]=wt;
  }
}
#endif // slow!

std::vector<float> model_predictions(1155, 0.5);
unsigned int prediction_index = 0;
float conversion_factor = 1.0 / 4095;

class Mixer {
  const int N, M, S;   // max inputs, max contexts, max context sets
  Array<short, 16> tx; // N inputs from add()
  Array<short, 16> wx; // N*M weights
  Array<int> cxt;  // S contexts
  int ncxt;        // number of contexts (0 to S)
  int base;        // offset of next context
  int nx;          // Number of inputs in tx, 0 to N
  Array<int> pr;   // last result (scaled 12 bits)
  Mixer* mp;       // points to a Mixer to combine results
public:
  Mixer(int n, int m, int s=1, int w=0);

  // Adjust weights to minimize coding cost of last prediction
  void update() {
    for (int i=0; i<ncxt; ++i) {
      int err=((y<<12)-pr[i])*7;
      assert(err>=-32768 && err<32768);
      train(&tx[0], &wx[cxt[i]*N], nx, err);
    }
    nx=base=ncxt=0;
  }

  // Input x (call up to N times)
  void add(int x) {
    assert(nx<N);
    model_predictions[prediction_index] = squash(x) * conversion_factor;
    ++prediction_index;
    tx[nx++]=x;
  }

  // Set a context (call S times, sum of ranges <= M)
  void set(int cx, int range) {
    assert(range>=0);
    assert(ncxt<S);
    assert(cx>=0);
    assert(base+cx<M);
    cxt[ncxt++]=base+cx;
    base+=range;
  }

  // predict next bit
  int p() {
    while (nx&7) tx[nx++]=0;  // pad
    if (mp) {  // combine outputs
      // printf("%d\n", prediction_index);
      prediction_index = 0;
      mp->update();
      for (int i=0; i<ncxt; ++i) {
        pr[i]=squash(dot_product(&tx[0], &wx[cxt[i]*N], nx)>>5);
        mp->add(stretch(pr[i]));
      }
      mp->set(0, 1);
      return mp->p();
    }
    else {  // S=1 context
      return pr[0]=squash(dot_product(&tx[0], &wx[0], nx)>>8);
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
  assert(n>0 && N>0 && (N&7)==0 && M>0);
  for (int i=0; i<S; ++i)
    pr[i]=2048;
  for (int i=0; i<N*M; ++i)
    wx[i]=w;
  if (S>1) mp=new Mixer(S, 1, 1, 0x7fff);
}

//////////////////////////// APM //////////////////////////////

// APM maps a probability and a context into a new probability
// that bit y will next be 1.  After each guess it updates
// its state to improve future guesses.  Methods:
//
// APM a(N) creates with N contexts, uses 66*N bytes memory.
// a.p(pr, cx, rate=7) returned adjusted probability in context cx (0 to
//   N-1).  rate determines the learning rate (smaller = faster, default 7).
//   Probabilities are scaled 12 bits (0-4095).

class APM {
  int index;     // last p, context
  const int N;   // number of contexts
  Array<U16> t;  // [N][33]:  p, context -> p
public:
  APM(int n);
  int p(int pr=2048, int cxt=0, int rate=7) {
    assert(pr>=0 && pr<4096 && cxt>=0 && cxt<N && rate>0 && rate<32);
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

// A StateMap maps a nonstationary counter state to a probability.
// After each mapping, the mapping is adjusted to improve future
// predictions.  Methods:
//
// sm.p(cx) converts state cx (0-255) to a probability (0-4095).

// Counter state -> probability * 4096
class StateMap {
protected:
  int cxt;  // context
  Array<U16> t; // 256 states -> probability * 64K
public:
  StateMap();
  int p(int cx) {
    assert(cx>=0 && cx<t.size());
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

// Hash 2-5 ints.
inline U32 hash(U32 a, U32 b, U32 c=0xffffffff, U32 d=0xffffffff,
    U32 e=0xffffffff) {
  U32 h=a*200002979u+b*30005491u+c*50004239u+d*70004807u+e*110002499u;
  return h^h>>9^a>>2^b>>3^c>>4^d>>5^e>>6;
}

///////////////////////////// BH ////////////////////////////////

// A BH maps a 32 bit hash to an array of B bytes (checksum and B-2 values)
//
// BH bh(N); creates N element table with B bytes each.
//   N must be a power of 2.  The first byte of each element is
//   reserved for a checksum to detect collisions.  The remaining
//   B-1 bytes are values, prioritized by the first value.  This
//   byte is 0 to mark an unused element.
//   
// bh[i] returns a pointer to the i'th element, such that
//   bh[i][0] is a checksum of i, bh[i][1] is the priority, and
//   bh[i][2..B-1] are other values (0-255).
//   The low lg(n) bits as an index into the table.
//   If a collision is detected, up to M nearby locations in the same
//   cache line are tested and the first matching checksum or
//   empty element is returned.
//   If no match or empty element is found, then the lowest priority
//   element is replaced.

// 2 byte checksum with LRU replacement (except last 2 by priority)
template <int B> class BH {
  enum {M=8};  // search limit
  Array<U8, 64> t; // elements
  U32 n; // size-1
public:
  BH(int i): t(i*B), n(i-1) {
    assert(B>=2 && i>0 && (i&(i-1))==0); // size a power of 2?
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
    if (*cp==chk) break;  // found
  }
  if (j==0) return p+1;  // front
  static U8 tmp[B];  // element to move to front
  if (j==M) {
    --j;
    memset(tmp, 0, B);
    //*(U16*)tmp=chk;
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
//
// A ContextMap maps contexts to a bit histories and makes predictions
// to a Mixer.  Methods common to all classes:
//
// ContextMap cm(M, C); creates using about M bytes of memory (a power
//   of 2) for C contexts.
// cm.set(cx);  sets the next context to cx, called up to C times
//   cx is an arbitrary 32 bit value that identifies the context.
//   It should be called before predicting the first bit of each byte.
// cm.mix(m) updates Mixer m with the next prediction.  Returns 1
//   if context cx is found, else 0.  Then it extends all the contexts with
//   global bit y.  It should be called for every bit:
//
//     if (bpos==0) 
//       for (int i=0; i<C; ++i) cm.set(cxt[i]);
//     cm.mix(m);
//
// The different types are as follows:
//
// - RunContextMap.  The bit history is a count of 0-255 consecutive
//     zeros or ones.  Uses 4 bytes per whole byte context.  C=1.
//     The context should be a hash.
// - SmallStationaryContextMap.  0 <= cx < M/512.
//     The state is a 16-bit probability that is adjusted after each
//     prediction.  C=1.
// - ContextMap.  For large contexts, C >= 1.  Context need not be hashed.

// Predict to mixer m from bit history state s, using sm to map s to
// a probability.
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

// A RunContextMap maps a context into the next byte and a repeat
// count up to M.  Size should be a power of 2.  Memory usage is 3M/4.
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

// Context is looked up directly.  m=size is power of 2 in bytes.
// Context should be < m/512.  High bits are discarded.
class SmallStationaryContextMap {
  Array<U16> t;
  int cxt;
  U16 *cp;
public:
  SmallStationaryContextMap(int m): t(m/2), cxt(0) {
    assert((m/2&m/2-1)==0); // power of 2?
    for (int i=0; i<t.size(); ++i)
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

// Context map for large contexts.  Most modeling uses this type of context
// map.  It includes a built in RunContextMap to predict the last byte seen
// in the same context, and also bit-level contexts that map to a bit
// history state.
//
// Bit histories are stored in a hash table.  The table is organized into
// 64-byte buckets alinged on cache page boundaries.  Each bucket contains
// a hash chain of 7 elements, plus a 2 element queue (packed into 1 byte) 
// of the last 2 elements accessed for LRU replacement.  Each element has
// a 2 byte checksum for detecting collisions, and an array of 7 bit history
// states indexed by the last 0 to 2 bits of context.  The buckets are indexed
// by a context ending after 0, 2, or 5 bits of the current byte.  Thus, each
// byte modeled results in 3 main memory accesses per context, with all other
// accesses to cache.
//
// On bits 0, 2 and 5, the context is updated and a new bucket is selected.
// The most recently accessed element is tried first, by comparing the
// 16 bit checksum, then the 7 elements are searched linearly.  If no match
// is found, then the element with the lowest priority among the 5 elements 
// not in the LRU queue is replaced.  After a replacement, the queue is
// emptied (so that consecutive misses favor a LFU replacement policy).
// In all cases, the found/replaced element is put in the front of the queue.
//
// The priority is the state number of the first element (the one with 0
// additional bits of context).  The states are sorted by increasing n0+n1
// (number of bits seen), implementing a LFU replacement policy.
//
// When the context ends on a byte boundary (bit 0), only 3 of the 7 bit
// history states are used.  The remaining 4 bytes implement a run model
// as follows: <count:7,d:1> <b1> <unused> <unused> where <b1> is the last byte
// seen, possibly repeated.  <count:7,d:1> is a 7 bit count and a 1 bit
// flag (represented by count * 2 + d).  If d=0 then <count> = 1..127 is the 
// number of repeats of <b1> and no other bytes have been seen.  If d is 1 then 
// other byte values have been seen in this context prior to the last <count> 
// copies of <b1>.
//
// As an optimization, the last two hash elements of each byte (representing
// contexts with 2-7 bits) are not updated until a context is seen for
// a second time.  This is indicated by <count,d> = <1,0> (2).  After update,
// <count,d> is updated to <2,0> or <1,1> (4 or 3).

class ContextMap {
  const int C;  // max number of contexts
  class E {  // hash element, 64 bytes
    U16 chk[7];  // byte context checksums
    U8 last;     // last 2 accesses (0-6) in low, high nibble
  public:
    U8 bh[7][7]; // byte context, 3-bit context -> bit history state
      // bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
      // bh[][0] is also a replacement priority, 0 = empty
    U8* get(U16 chk);  // Find element (0-6) matching checksum.
      // If not found, insert or replace lowest priority (not last).
  };
  Array<E, 64> t;  // bit histories for bits 0-1, 2-4, 5-7
    // For 0-1, also contains a run count in bh[][4] and value in bh[][5]
    // and pending update count in bh[7]
  Array<U8*> cp;   // C pointers to current bit history
  Array<U8*> cp0;  // First element of 7 element array containing cp[i]
  Array<U32> cxt;  // C whole byte contexts (hashes)
  Array<U8*> runp; // C [0..3] = count, value, unused, unused
  StateMap *sm;    // C maps of state -> p
  int cn;          // Next context to set by set()
  void update(U32 cx, int c);  // train model that context cx predicts c
  int mix1(Mixer& m, int cc, int bp, int c1, int y1);
    // mix() with global context passed as arguments to improve speed.
public:
  ContextMap(int m, int c=1);  // m = memory in bytes, a power of 2, C = c
  ~ContextMap();
  void set(U32 cx, int next=-1);   // set next whole byte context to cx
    // if next is 0 then set order does not matter
  int mix(Mixer& m) {return mix1(m, c0, bpos, buf(1), y);}
};

// Find or create hash element matching checksum ch
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

// Construct using m bytes of memory for c contexts
ContextMap::ContextMap(int m, int c): C(c), t(m>>6), cp(c), cp0(c),
    cxt(c), runp(c), cn(0) {
  assert(m>=64 && (m&m-1)==0);  // power of 2?
  assert(sizeof(E)==64);
  sm=new StateMap[C];
  for (int i=0; i<C; ++i) {
    cp0[i]=cp[i]=&t[0].bh[0][0];
    runp[i]=cp[i]+3;
  }
}

ContextMap::~ContextMap() {
  delete[] sm;
}

// Set the i'th context to cx
inline void ContextMap::set(U32 cx, int next) {
  int i=cn++;
  i&=next;
  assert(i>=0 && i<C);
  cx=cx*987654323+i;  // permute (don't hash) cx to spread the distribution
  cx=cx<<16|cx>>16;
  cxt[i]=cx*123456791+i;
}

// Update the model with bit y1, and predict next bit to mixer m.
// Context: cc=c0, bp=bpos, c1=buf(1), y1=y.
int ContextMap::mix1(Mixer& m, int cc, int bp, int c1, int y1) {

  // Update model with y
  int result=0;
  for (int i=0; i<cn; ++i) {
    if (cp[i]) {
      assert(cp[i]>=&t[0].bh[0][0] && cp[i]<=&t[t.size()-1].bh[6][6]);
      assert((long(cp[i])&63)>=15);
      int ns=nex(*cp[i], y1);
      if (ns>=204 && rnd() << ((452-ns)>>3)) ns-=4;  // probabilistic increment
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
        if (runp[i][0]==0)  // new context
          runp[i][0]=2, runp[i][1]=c1;
        else if (runp[i][1]!=c1)  // different byte in context
          runp[i][0]=1, runp[i][1]=c1;
        else if (runp[i][0]<254)  // same byte in context
          runp[i][0]+=2;
        else if (runp[i][0]==255)
          runp[i][0]=128;
        runp[i]=cp0[i]+3;
      }
    }

    // predict from last byte in context
    int rc=runp[i][0];  // count*2, +1 if 2 different bytes seen
    if ((runp[i][1]+256)>>(8-bp)==cc) {
      int b=((runp[i][1]>>(7-bp))&1)*2-1;  // predicted bit + for 1, - for 0
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

//////////////////////////// Models //////////////////////////////

// All of the models below take a Mixer as a parameter and write
// predictions to it.

//////////////////////////// matchModel ///////////////////////////

// matchModel() finds the longest matching context and returns its length

int matchModel(Mixer& m) {
  const int MAXLEN=65534;  // longest allowed match + 1
  static Array<int> t(MEM);  // hash table of pointers to contexts
  static int h=0;  // hash of last 7 bytes
  static int ptr=0;  // points to next byte of match if any
  static int len=0;  // length of match, or 0 if no match
  static int result=0;
  
  static SmallStationaryContextMap scm1(0x20000);

  if (!bpos) {
    h=(h*997*8+buf(1)+1)&(t.size()-1);  // update context hash
    if (len) ++len, ++ptr;
    else {  // find match
      ptr=t[h];
      if (ptr && pos-ptr<buf.size())
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

// Model a 1728 by 2376 2-color CCITT bitmap image, left to right scan,
// MSB first (216 bytes per row, 513216 bytes total).  Insert predictions
// into m.

void picModel(Mixer& m) {
  static U32 r0, r1, r2, r3;  // last 4 rows, bit 8 is over current pixel
  static Array<U8> t(0x10200);  // model: cxt -> state
  const int N=3;  // number of contexts
  static int cxt[N];  // contexts
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

// Model English text (words and columns/end of line)
static U32 frstchar=0,spafdo=0,spaces=0,spacecount=0, words=0,wordcount=0,wordlen=0,wordlen1=0;
void wordModel(Mixer& m) {
	static U32 word0=0, word1=0, word2=0, word3=0, word4=0, word5=0;  // hashes
	static U32 xword0=0,xword1=0,xword2=0,cword0=0,ccword=0;
	static U32 number0=0, number1=0;  // hashes
	static U32 text0=0;  // hash stream of letters
	static ContextMap cm(MEM*31, 45);
	static int nl1=-3, nl=-2;  // previous, current newline position
	static U32 mask = 0;
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
			word0=word0*263*32+c;
			text0=text0*997*16+c;
			wordlen++;
			wordlen=min(wordlen,45);
			f=0;
		}
		else {
			if (word0) {
				word5=word4*23;
				word4=word3*19;
				word3=word2*17;
				word2=word1*13;
				word1=word0*11;
				wordlen1=wordlen;
				if (c==':') cword0=word0;
				if (c==']'&& (frstchar!=':' || frstchar!='*')) xword0=word0;
			//	if (c==0x27) xword0=word0;
			    if (c=='=') cword0=word0;
				ccword=0;
				word0=wordlen=0;
				if((c=='.'||c=='!'||c=='?') && buf(2)!=10) f=1; 
				
			}
			if ((c4&0xFFFF)==0x3D3D) xword1=word1,xword2=word2; // '=='
				if ((c4&0xFFFF)==0x2727) xword1=word1,xword2=word2; // ''
			if (c==32 || c==10 ) { ++spaces, ++spacecount; if (c==10 ) nl1=nl, nl=pos-1;}
			else if (c=='.' || c=='!' || c=='?' || c==',' || c==';' || c==':') spafdo=0,ccword=c*31; 
			else { ++spafdo; spafdo=min(63,spafdo); }
		}
		if (c>='0' && c<='9') {
			number0=number0*263*32+c;
		}
		else if (number0) {
			number1=number0*11;
			number0=0,ccword=0;
		}
	
		col=min(255, pos-nl);
		int above=buf[nl1+col]; // text column context
		if (col<=2) frstchar=(col==2?min(c,96):0);
//	cm.set(spafdo|col<<8);
		cm.set(spafdo|spaces<<8);
		if (frstchar=='[' && c==32)	{if(buf(3)==']' || buf(4)==']' ) frstchar=96,xword0=0;}
		cm.set(frstchar<<11|c);
		cm.set(col<<8|frstchar);
		cm.set(spaces<<8|(words&255));
		
		cm.set(number0+word2);
		cm.set(number0+word1);
		cm.set(number1+c);
		cm.set(number0+number1);
		cm.set(word0+number1);

		cm.set(frstchar<<7);
//	cm.set(wordlen<<16|c);
		cm.set(wordlen1<<8|col);
		cm.set(c*64+spacecount/2);
		U32 h=wordcount*64+spacecount;
		cm.set((c<<13)+h);
//	cm.set(h); // 
		cm.set(h+spafdo*8);

		U32 d=c4&0xf0ff;
		cm.set(d<<9|frstchar);

		h=word0*271;
		cm.set(h+ccword);
		h=h+buf(1);
		cm.set(h);
		cm.set(word0);
		//cm.set(word1+(c==32));
		cm.set(h+word1);
		cm.set(word0+word1*31);
		cm.set(h+word1+word2*29);
		cm.set(text0&0xffffff);
		cm.set(text0&0xfffff);
          cm.set(word0+xword0*31);
		  cm.set(word0+xword1*31);
		  cm.set(word0+xword2*31);
		  cm.set(frstchar+xword2*31);
        
        cm.set(word0+cword0*31);
		cm.set(number0+cword0*31);
		cm.set(h+word2);
		cm.set(h+word3);
		cm.set(h+word4);
		cm.set(h+word5);
//		cm.set(buf(1)|buf(3)<<8|buf(5)<<16);
//		cm.set(buf(2)|buf(4)<<8|buf(6)<<16);
		cm.set(h+word1+word3);
		cm.set(h+word2+word3);
		if (f) {
			word5=word4*29;
			word4=word3*31;
			word3=word2*37;
			word2=word1*41;
			word1='.';
		}
		cm.set(col<<16|buf(1)<<8|above);
		cm.set(buf(1)<<8|above);
		cm.set(col<<8|buf(1));
		cm.set(col*(c==32));
		cm.set(col);
		
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
 
    cm.set(mask);
    cm.set((mask<<8)|buf(1));
    //cm.set((mask&0xff)|col<<8);
    cm.set((mask<<17)|(buf(2)<<8)|buf(3));
    cm.set((mask&0x1ff)|((f4&0x00fff0)<<9));
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

// Model 2-D data with fixed record length.  Also order 1-2 models
// that include the distance to the last match.

void recordModel(Mixer& m) {
  static int cpos1[256] , cpos2[256], cpos3[256], cpos4[256];
  static int wpos1[0x10000]; // buf(1..2) -> last position
  static int rlen=2, rlen1=3, rlen2=4;  // run length and 2 candidates
  static int rcount1=0, rcount2=0;  // candidate counts
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
    assert(rlen>0);
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

    //cp.set(w*16);
    //cp.set(d*32);
    //cp.set(c*64);
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
  static int wpos1[0x10000]; // buf(1..2) -> last position
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

// Model order 1-2 contexts with gaps.

void sparseModel(Mixer& m, int seenbefore, int howmany) {
  static ContextMap cm(MEM*2, 48);
  static int mask = 0;

  if (bpos==0) {

    cm.set( c4&0x00f0f0f0);
    cm.set((c4&0xf0f0f0f0)+1);
    cm.set((c4&0x00f8f8f8)+2);
    cm.set((c4&0xf8f8f8f8)+3);
    cm.set((c4&0x00e0e0e0)+4);
    cm.set((c4&0xe0e0e0e0)+5);
    cm.set((c4&0x00f0f0ff)+6);

    cm.set(seenbefore);
    cm.set(howmany);
    cm.set(c4&0x00ff00ff);
    cm.set(c4&0xff0000ff);
    cm.set(buf(1)|buf(5)<<8);
    cm.set(buf(1)|buf(6)<<8);
    cm.set(buf(3)|buf(6)<<8);
    cm.set(buf(4)|buf(8)<<8);
    
    for (int i=1; i<8; ++i) {
      cm.set((buf(i+1)<<8)|buf(i+2));
      cm.set((buf(i+1)<<8)|buf(i+3));
      cm.set(seenbefore|buf(i)<<8);
    }

    int fl = 0;
    if( c4&(0xff != 0 )){
           if( isalpha( c4&0xff ) ) fl = 1;
      else if( ispunct( c4&0xff ) ) fl = 2;
      else if( isspace( c4&0xff ) ) fl = 3;
      else if( c4&(0xff == 0xff) ) fl = 4;
      else if( c4&(0xff < 16) ) fl = 5;
      else if( c4&(0xff < 64) ) fl = 6;
      else fl = 7;
    }
    mask = (mask<<3)|fl;
    cm.set(mask);
    cm.set(mask<<8|buf(1));
    cm.set(mask<<17|buf(2)<<8|buf(3));
    cm.set((mask&0x1ff)|((c4&0xf0f0f0f0)<<9));
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

// Model for modelling distances between symbols

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

//////////////////////////// exeModel /////////////////////////

// Model x86 code.  The contexts are sparse containing only those
// bits relevant to parsing (2 prefixes, opcode, and mod and r/m fields
// of modR/M byte).

// Get context at buf(i) relevant to parsing 32-bit x86 code
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

// The context is a byte string history that occurs within a
// 1 or 2 byte context.

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

// Model using DMC.  The bitwise context is represented by a state graph,
// initilaized to a bytewise order 1 model as in 
// http://plg.uwaterloo.ca/~ftp/dmc/dmc.c but with the following difference:
// - It uses integer arithmetic.
// - The threshold for cloning a state increases as memory is used up.
// - Each state maintains both a 0,1 count and a bit history (as in a
//   context model).  The 0,1 count is best for stationary data, and the
//   bit history for nonstationary data.  The bit history is mapped to
//   a probability adaptively using a StateMap.  The two computed probabilities
//   are combined.
// - When memory is used up the state graph is reinitialized to a bytewise
//   order 1 context as in the original DMC.  However, the bit histories
//   are not cleared.

struct DMCNode {  // 12 bytes
  unsigned int nx[2];  // next pointers
  U8 state;  // bit history
  unsigned int c0:12, c1:12;  // counts * 256
};

void dmcModel(Mixer& m) {
  static int top=0, curr=0;  // allocated, current node
  static Array<DMCNode> t(MEM*2);  // state graph
  static StateMap sm;
  static int threshold=256;

  // clone next state
  if (top>0 && top<t.size()) {
    int next=t[curr].nx[y];
    int n=y?t[curr].c1:t[curr].c0;
    int nn=t[next].c0+t[next].c1;
    if (n>=threshold*2 && nn-n>=threshold*3) {
      int r=n*4096/nn;
      assert(r>=0 && r<=4096);
      t[next].c0 -= t[top].c0 = t[next].c0*r>>12;
      t[next].c1 -= t[top].c1 = t[next].c1*r>>12;
      t[top].nx[0]=t[next].nx[0];
      t[top].nx[1]=t[next].nx[1];
      t[top].state=t[next].state;
      t[curr].nx[y]=top;
      ++top;
      if (top==(t.size()*4)/8) { //        5/8, 4/8, 3/8
        threshold=512;
      } else if (top==(t.size()*6)/8) { // 6/8, 6/8, 7/8
        threshold=768;
      }
    }
  }

  // Initialize to a bytewise order 1 model at startup or when flushing memory
  if (top==t.size() && bpos==1) top=0;
  if (top==0) {
    assert(t.size()>=65536);
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

  // update count, state
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

// This combines all the context models with a Mixer.

U32 last_prediction = 2048;

int contextModel2() {
  static ContextMap cm(MEM*31, 9);
  static RunContextMap rcm7(MEM), rcm9(MEM), rcm10(MEM);
  static Mixer m(1200, 10800, 8, 128);
  static U32 cxt[16];  // order 0-11 contexts
  static Filetype filetype=EXE;

  m.update();
  m.add(256);

  // Test for special file types
  int ismatch=ilog(matchModel(m));  // Length of longest matching context

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

  order = order-3;
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

// A Predictor estimates the probability that the next bit of
// uncompressed data is 1.  Methods:
// p() returns P(1) as a 12 bit number (0-4095).
// update(y) trains the predictor with the actual bit (0 or 1).

class Predictor {
  int pr;  // next prediction
public:
  Predictor();
  int p() const {assert(pr>=0 && pr<4096); return pr;}
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
