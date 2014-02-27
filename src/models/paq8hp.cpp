// This is adapted from paq8hp12any

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

#include "paq8hp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <algorithm>
#define NDEBUG  // remove for debugging (turns on Array bound checks)
#include <assert.h>

#ifndef DEFAULT_OPTION
#define DEFAULT_OPTION 8
#endif

#ifndef NOASM
#define NOASM
#endif

namespace {

// 8, 16, 32 bit unsigned types (adjust as appropriate)
typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;

// min, max functions
#ifndef min
inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}
#endif

// Error handler: print message if any, and exit
void quit(const char* message=0) {
  throw message;
}

#if 0
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
#endif

typedef enum {DEFAULT, JPEG, EXE, BINTEXT, TEXT } Filetype; 

#define preprocFlag 1220

#define OPTION_UTF8							1
#define OPTION_USE_NGRAMS					2
#define OPTION_CAPITAL_CONVERSION			4
#define OPTION_WORD_SURROROUNDING_MODELING	8
#define OPTION_SPACE_AFTER_EOL				16
#define OPTION_EOL_CODING					32
#define OPTION_NORMAL_TEXT_FILTER			64
#define OPTION_USE_DICTIONARY				128
#define OPTION_RECORD_INTERLEAVING			256
#define OPTION_DNA_QUARTER_BYTE				512
#define OPTION_TRY_SHORTER_WORD				1024
#define OPTION_TO_LOWER_AFTER_PUNCTUATION	2048
#define OPTION_SPACELESS_WORDS				4096
#define OPTION_ADD_SYMBOLS_0_5				8192
#define OPTION_ADD_SYMBOLS_14_31			16384
#define OPTION_ADD_SYMBOLS_A_Z				32768
#define OPTION_ADD_SYMBOLS_MISC				65536
#define OPTION_SPACE_AFTER_CC_FLAG			131072
#define IF_OPTION(option) ((preprocFlag & option)!=0)

//////////////////////// Program Checker /////////////////////

// Track time and memory used
class ProgramChecker {
  int memused;  // bytes allocated by Array<T> now
  clock_t start_time;  // in ticks
public:
  int maxmem;   // most bytes allocated ever
  void alloc(int n) {  // report memory allocated, may be negative
    memused+=n;
    if (memused>maxmem) maxmem=memused;
  }
  ProgramChecker(): memused(0), maxmem(0) {
    start_time=clock();
    assert(sizeof(U8)==1);
    assert(sizeof(U16)==2);
    assert(sizeof(U32)==4);
    assert(sizeof(short)==2);
    assert(sizeof(int)==4);
  }
  void print() const {  // print time and memory used
    printf("Time %1.2f sec, used %d bytes of memory\n",
      double(clock()-start_time)/CLOCKS_PER_SEC, maxmem);
  }
} programChecker;

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
#ifndef NDEBUG
    if (i<0 || i>=n) fprintf(stderr, "%d out of bounds %d\n", i, n), quit();
#endif
    return data[i];
  }
  const T& operator[](int i) const {
#ifndef NDEBUG
    if (i<0 || i>=n) fprintf(stderr, "%d out of bounds %d\n", i, n), quit();
#endif
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
  if (savedata && saveptr) {
    memcpy(data, savedata, sizeof(T)*min(i, saven));
    programChecker.alloc(-ALIGN-n*sizeof(T));
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
  programChecker.alloc(sz);
  ptr = (char*)calloc(sz, 1);
  if (!ptr) quit("Out of memory");
  data = (ALIGN ? (T*)(ptr+ALIGN-(((long)ptr)&(ALIGN-1))) : (T*)ptr);
  assert((char*)data>=ptr && (char*)data<=ptr+ALIGN);
}

template<class T, int ALIGN> Array<T, ALIGN>::~Array() {
  programChecker.alloc(-ALIGN-n*sizeof(T));
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
  U32 table[64];
  int i;
public:
  Random() {
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
U32 b1=0, b2=0, b3=0, b4=0, b5=0, b6=0, b7=0, b8=0, tt=0, c4=0, x4=0, x5=0, w4=0, w5=0, f4=0; // Last 4 whole bytes, packed.  Last byte is bits 0-7.
int order, bpos=0, cxtfl=3, sm_shft=7, sm_add=65535+127, sm_add_y=0; // bits in c0 (0 to 7)
Buf buf;  // Rotating input queue set by Predictor

///////////////////////////// ilog //////////////////////////////

// ilog(x) = round(log2(x) * 16), 0 <= x < 64K
class Ilog {
  U8 t[65536];
public:
  int operator()(U16 x) const {return t[x];}
  Ilog();
} ilog;

// Compute lookup table by numerical integration of 1/x
Ilog::Ilog() {
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

#if 1 // change to #if 0 to generate this table at run time (4% slower)
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

// The code used to generate the above table at run time (4% slower).
// To print the table, uncomment the 4 lines of print statements below.
// In this code x,y = n0,n1 is the number of 0,1 bits represented by a state.
#else

class StateTable {
  Array<U8> ns;  // state*4 -> next state if 0, if 1, n0, n1
  enum {B=5, N=64}; // sizes of b, t
  static const int b[B];  // x -> max y, y -> max x
  static U8 t[N][N][2];  // x,y -> state number, number of states
  int num_states(int x, int y);  // compute t[x][y][1]
  void discount(int& x);  // set new value of x after 1 or y after 0
  void next_state(int& x, int& y, int b);  // new (x,y) after bit b
public:
  int operator()(int state, int sel) {return ns[state*4+sel];}
  StateTable();
} nex;

const int StateTable::b[B]={42,41,13,6,5};  // x -> max y, y -> max x
U8 StateTable::t[N][N][2];

int StateTable::num_states(int x, int y) {
  if (x<y) return num_states(y, x);
  if (x<0 || y<0 || x>=N || y>=N || y>=B || x>=b[y]) return 0;

  // States 0-30 are a history of the last 0-4 bits
  if (x+y<=4) {  // x+y choose x = (x+y)!/x!y!
    int r=1;
    for (int i=x+1; i<=x+y; ++i) r*=i;
    for (int i=2; i<=y; ++i) r/=i;
    return r;
  }

  // States 31-255 represent a 0,1 count and possibly the last bit
  // if the state is reachable by either a 0 or 1.
  else
    return 1+(y>0 && x+y<16);
}

// New value of count x if the opposite bit is observed
void StateTable::discount(int& x) {
  if (x>2) x=ilog(x)/6-1;
}

// compute next x,y (0 to N) given input b (0 or 1)
void StateTable::next_state(int& x, int& y, int b) {
  if (x<y)
    next_state(y, x, 1-b);
  else {
    if (b) {
      ++y;
      discount(x);
    }
    else {
      ++x;
      discount(y);
    }
    while (!t[x][y][1]) {
      if (y<2) --x;
      else {
        x=(x*(y-1)+(y/2))/y;
        --y;
      }
    }
  }
}

// Initialize next state table ns[state*4] -> next if 0, next if 1, x, y
StateTable::StateTable(): ns(1024) {

  // Assign states
  int state=0;
  for (int i=0; i<256; ++i) {
    for (int y=0; y<=i; ++y) {
      int x=i-y;
      int n=num_states(x, y);
      if (n) {
        t[x][y][0]=state;
        t[x][y][1]=n;
        state+=n;
      }
    }
  }

  // Print/generate next state table
  state=0;
  for (int i=0; i<N; ++i) {
    for (int y=0; y<=i; ++y) {
      int x=i-y;
      for (int k=0; k<t[x][y][1]; ++k) {
        int x0=x, y0=y, x1=x, y1=y;  // next x,y for input 0,1
        int ns0=0, ns1=0;
        if (state<15) {
          ++x0;
          ++y1;
          ns0=t[x0][y0][0]+state-t[x][y][0];
          ns1=t[x1][y1][0]+state-t[x][y][0];
          if (x>0) ns1+=t[x-1][y+1][1];
          ns[state*4]=ns0;
          ns[state*4+1]=ns1;
          ns[state*4+2]=x;
          ns[state*4+3]=y;
        }
        else if (t[x][y][1]) {
          next_state(x0, y0, 0);
          next_state(x1, y1, 1);
          ns[state*4]=ns0=t[x0][y0][0];
          ns[state*4+1]=ns1=t[x1][y1][0]+(t[x1][y1][1]>1);
          ns[state*4+2]=x;
          ns[state*4+3]=y;
        }
          // uncomment to print table above
//        printf("{%3d,%3d,%2d,%2d},", ns[state*4], ns[state*4+1], 
//          ns[state*4+2], ns[state*4+3]);
//        if (state%4==3) printf(" // %d-%d\n  ", state-3, state);
        assert(state>=0 && state<256);
        assert(t[x][y][1]>0);
        assert(t[x][y][0]<=state);
        assert(t[x][y][0]+t[x][y][1]>state);
        assert(t[x][y][1]<=6);
        assert(t[x0][y0][1]>0);
        assert(t[x1][y1][1]>0);
        assert(ns0-t[x0][y0][0]<t[x0][y0][1]);
        assert(ns0-t[x0][y0][0]>=0);
        assert(ns1-t[x1][y1][0]<t[x1][y1][1]);
        assert(ns1-t[x1][y1][0]>=0);
        ++state;
      }
    }
  }
//  printf("%d states\n", state); exit(0);  // uncomment to print table above
}

#endif

///////////////////////////// Squash //////////////////////////////

// return p = 1/(1 + exp(-d)), d scaled by 8 bits, p scaled by 12 bits
int squash(int d) {
  static const int t[33]={
    1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
    1546,2047,2549,2994,3348,3607,3785,3901,3975,4022,
    4050,4068,4079,4085,4089,4092,4093,4094};
  if (d>2047) return 4095;
  if (d<-2047) return 0;
  int w=d&127;
  d=(d>>7)+16;
  return (t[d]*(128-w)+t[(d+1)]*w+64) >> 7;
}

//////////////////////////// Stretch ///////////////////////////////

// Inverse of squash. d = ln(p/(1-p)), d scaled by 8 bits, p by 12 bits.
// d has range -2047 to 2047 representing -8 to 8.  p has range 0 to 4095.

class Stretch {
  short t[4096];
public:
  Stretch();
  int operator()(int p) const {
    assert(p>=0 && p<4096);
    return t[p];
  }
} stretch;

Stretch::Stretch() {
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
#ifdef NOASM  // no assembly language
int dot_product(short *t, short *w, int n) {
  int sum=0;
  n=(n+7)&-8;
  for (int i=0; i<n; i+=2)
    sum+=(t[i]*w[i]+t[i+1]*w[i+1]) >> 8;
  return sum;
}
#else  // The NASM version uses MMX and is about 8 times faster.
extern "C" int dot_product(short *t, short *w, int n);  // in NASM
#endif

// Train neural network weights w[n] given inputs t[n] and err.
// w[i] += t[i]*err, i=0..n-1.  t, w, err are signed 16 bits (+- 32K).
// err is scaled 16 bits (representing +- 1/2).  w[i] is clamped to +- 32K
// and rounded.  n is rounded up to a multiple of 8.
#ifdef NOASM
void train(short *t, short *w, int n, int err) {
  n=(n+7)&-8;
  for (int i=0; i<n; ++i) {
    int wt=w[i]+((((t[i]*err*2)>>16)+1)>>1);
    if (wt<-32768) wt=-32768;
    if (wt>32767) wt=32767;
    w[i]=wt;
  }
}
#else
extern "C" void train(short *t, short *w, int n, int err);  // in NASM
#endif

class Mixer {
  const int N, M, S;   // max inputs, max contexts, max context sets
  Array<short, 16> wx; // N*M weights
  Array<int> cxt;  // S contexts
  int ncxt;        // number of contexts (0 to S)
  int base;        // offset of next context
  Array<int> pr;   // last result (scaled 12 bits)
  Mixer* mp;       // points to a Mixer to combine results
public:
  Array<short, 16> tx; // N inputs from add()
  int nx;          // Number of inputs in tx, 0 to N
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

  void update2() {
    train(&tx[0], &wx[0], nx, ((y<<12)-base)*3/2);
    nx=0;
  }

  // Input x (call up to N times)
  void add(int x) {
    assert(nx<N);
    tx[nx++]=x;
  }

  void mul(int x) {
    int z=tx[nx];
    z=z*x/4;
    tx[nx++]=z;
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
      mp->update2();
      for (int i=0; i<ncxt; ++i) {
        int dp=dot_product(&tx[0], &wx[cxt[i]*N], nx);
	dp=(dp*9)>>9;
        pr[i] = squash(dp);
        mp->add(dp);
      }
      return mp->p();
    }
    else {  // S=1 context
	int z=dot_product(&tx[0], &wx[0], nx);
	base=squash( (z*15) >>13);
	return squash(z>>9);
    }
  }
  ~Mixer();
};

Mixer::~Mixer() {
  delete mp;
}


Mixer::Mixer(int n, int m, int s, int w):
    N((n+7)&-8), M(m), S(s), wx(N*M),
    cxt(S), ncxt(0), base(0), pr(S), mp(0), tx(N), nx(0) {
  assert(n>0 && N>0 && (N&7)==0 && M>0);
  int i;
  for (i=0; i<S; ++i)
    pr[i]=2048;
  for (i=0; i<N*M; ++i)
    wx[i]=w;
  if (S>1) mp=new Mixer(S, 1, 1, 0x7fff);
}

//////////////////////////// APM //////////////////////////////

// APM maps a probability and a context into a new probability
// that bit y will next be 1.  After each guess it updates
// its state to improve future guesses.  Methods:
//
// APM a(N) creates with N contexts, uses 66*N bytes memory.
// a.p(pr, cx, rate=8) returned adjusted probability in context cx (0 to
//   N-1).  rate determines the learning rate (smaller = faster, default 8).
//   Probabilities are scaled 16 bits (0-65535).

class APM {
  int index;     // last p, context
//const int N;   // number of contexts
  Array<U16> t;  // [N][33]:  p, context -> p
public:
  APM(int n);
  int p(int pr=2048, int cxt=0, int rate=8) {
    assert(pr>=0 && pr<4096 && cxt>=0 && cxt<N && rate>0 && rate<32);
    pr=stretch(pr);
    int g=(y<<16)+(y<<rate)-y*2;
    t[index]   += (g-t[index])   >> rate;
    t[index+1] += (g-t[index+1]) >> rate;
    const int w=pr&127;  // interpolation weight (33 points)
    index=((pr+2048)>>7)+cxt*33;
    return (t[index]*(128-w)+t[index+1]*w) >> 11;
  }
};

// maps p, cxt -> p initially
APM::APM(int n): index(0), t(n*33) {
    for (int j=0; j<33; ++j) t[j]=squash((j-16)*128)*16;
    for (int i=33; i<n*33; ++i) t[i]=t[i-33];
}

//////////////////////////// StateMap //////////////////////////

// A StateMap maps a nonstationary counter state to a probability.
// After each mapping, the mapping is adjusted to improve future
// predictions.  Methods:
//
// sm.p(cx) converts state cx (0-255) to a probability (0-4095).

// Counter state -> probability * 256
class StateMap {
protected:
  int cxt;  // context
  U16 t[256]; // 256 states -> probability * 64K
public:
  StateMap();
  int p(int cx) {
    assert(cx>=0 && cx<t.size());
    int q=t[cxt];
    t[cxt]=q + ( (sm_add_y - q) >> sm_shft);
    return t[cxt=cx] >> 4;
  }
};

StateMap::StateMap(): cxt(0) {
  for (int i=0; i<256; ++i) {
    int n0=nex(i,2);
    int n1=nex(i,3);
    if (n0==0) n1*=128;
    if (n1==0) n0*=128;
    t[i] = 65536*(n1+1)/(n0+n1+2);
  }
}

//////////////////////////// hash //////////////////////////////

// Hash 2-5 ints.
#if 0
inline U32 hash(U32 a, U32 b, U32 c=0xffffffff, U32 d=0xffffffff,
    U32 e=0xffffffff) {
  U32 h=a*200002979u+b*30005491u+c*50004239u+d*70004807u+e*110002499u;
  return h^h>>9^a>>2^b>>3^c>>4^d>>5^e>>6;
}
#else
inline U32 hash(U32 a, U32 b, U32 c=0xffffffff) {
  U32 h=a*110002499u+b*30005491u+c*50004239u; //+d*70004807u+e*110002499u;
  return h^h>>9^a>>3^b>>3^c>>4;
}
#endif

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
  enum {M=7};  // search limit
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
  U8 *p=&t[i*B-B];
  int j;
  for (j=0; j<M; ++j) {
    p+=B;
    if (p[2]==0) { *(U16*)p=chk; break; }
    if (*(U16*)p==chk) break;  // found
  }
  if (j==0) return p;  // front
#if 0
  static U8 tmp[B];  // element to move to front
  if (j==M) {
    --j;
    memset(tmp, 0, B);
    *(U16*)tmp=chk;
    if ( /*M>2&&*/ p[2]>p[2-B]) --j;
  }
  else memcpy(tmp, p, B);
  p = &t[i*B];
  memmove(p+B, p, j*B);
  memcpy(p, tmp, B);
#else
  if (j==M) {
    --j;
    if ( /*M>2&&*/ p[2]>p[-2]) --j;
  }
  else chk=*(int*)p;
  p = &t[i*4];
  memmove(p+4, p, j*4);
  *(int*)p=chk;
#endif
  return p;
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
  int st=stretch(p1);
 if (cxtfl) {
  m.add(st/4);
  int p0=4095-p1;
  m.add((p1-p0)*3/64);
  m.add(st*(n1-n0)*3/16);
  m.add(((p1&n0)-(p0&n1)) /16);
  m.add(((p0&n0)-(p1&n1))*7/64);
  return s>0;
 }
  m.add(st*9/32);
  m.add(st*(n1-n0)*3/16);
  int p0=4095-p1;
  m.add(((p1&n0)-(p0&n1)) /16);
  m.add(((p0&n0)-(p1&n1))*7/64);
  return s>0;
}

// A RunContextMap maps a context into the next byte and a repeat
// count up to M.  Size should be a power of 2.  Memory usage is 3M/4.
class RunContextMap {
  BH<4> t;
  U8 *cp;
  int mulc;
public:
  RunContextMap(int m, int c): t(m/4), mulc(c) {cp=t[0]+2;}
  void set(U32 cx) {  // update count
    if (cp[0]==0 || cp[1]!=b1) cp[0]=1, cp[1]=b1;
    else if (cp[0]<255) ++cp[0];
    cp=t[cx]+2;
  }
  int p() {  // predict next bit
    if ((cp[1]+256)>>(8-bpos)==c0)
      return (((cp[1]>>(7-bpos))&1)*2-1)*ilog(cp[0]+1)*mulc;
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
  int cxt, mulc;
  U16 *cp;
public:
  SmallStationaryContextMap(int m, int c): t(m/2), cxt(0), mulc(c) {
    assert((m/2&m/2-1)==0); // power of 2?
    for (int i=0; i<t.size(); ++i)
      t[i]=32768;
    cp=&t[0];
  }
  void set(U32 cx) {
    cxt=(cx*256)&(t.size()-256);
  }
  void mix(Mixer& m/*, int rate=7*/) {
#if 1
    if (pos<4000000)
    *cp += ((y<<16)-(*cp)+(1<<8)) >> 9;
    else
    *cp += ((y<<16)-(*cp)+(1<<9)) >> 10;
#else
    int q=*cp;
    if(y)
    q+=65790-q >> 8;
    else
    q+= -q >> 8;
    *cp=q;
#endif
    cp=&t[cxt+c0];
    m.add(stretch(*cp>>4)*mulc/32);
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
// as follows: <count:7,d:1> <b1> <b2> <b3> where <b1> is the last byte
// seen, possibly repeated, and <b2> and <b3> are the two bytes seen
// before the first <b1>.  <count:7,d:1> is a 7 bit count and a 1 bit
// flag.  If d=0 then <count> = 1..127 is the number of repeats of <b1>
// and no other bytes have been seen, and <b2><b3> are not used.
// If <d> = 1 then the history is <b3>, <b2>, and <count> - 2 repeats
// of <b1>.  In this case, <b3> is valid only if <count> >= 3 and
// <b2> is valid only if <count> >= 2.
//
// As an optimization, the last two hash elements of each byte (representing
// contexts with 2-7 bits) are not updated until a context is seen for
// a second time.  This is indicated by <count,d> = <1,0>.  After update,
// <count,d> is updated to <2,0> or <2,1>.

class ContextMap {
  const int C, Sz;  // max number of contexts
  class E {  // hash element, 64 bytes
    U16 chk[7];  // byte context checksums
    U8 last;     // last 2 accesses (0-6) in low, high nibble
  public:
    U8 bh[7][7]; // byte context, 3-bit context -> bit history state
      // bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
      // bh[][0] is also a replacement priority, 0 = empty
    U8* get(U16 chk, int i);  // Find element (0-6) matching checksum.
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
  int mix1(Mixer& m, int cc, int c1, int y1);
    // mix() with global context passed as arguments to improve speed.
public:
  ContextMap(int m, int c=1);  // m = memory in bytes, a power of 2, C = c
  void set(U32 cx);   // set next whole byte context
  int mix(Mixer& m) {return mix1(m, c0, b1, y);}
};

// Find or create hash element matching checksum ch
inline U8* ContextMap::E::get(U16 ch, int j) {
  ch+=j;
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
ContextMap::ContextMap(int m, int c): C(c), Sz((m>>6)-1), t(m>>6), cp(c), cp0(c),
    cxt(c), runp(c), cn(0) {
  assert(m>=64 && (m&m-1)==0);  // power of 2?
  assert(sizeof(E)==64);
  sm=new StateMap[C];
  for (int i=0; i<C; ++i) {
    cp0[i]=cp[i]=&t[0].bh[0][0];
    runp[i]=cp[i]+3;
  }
}

// Set the i'th context to cx
inline void ContextMap::set(U32 cx) {
  int i=cn++;
  assert(i>=0 && i<C);
  cx=cx*123456791+i;  // permute (don't hash) cx to spread the distribution
  cx=cx<<16|cx>>16;
  cxt[i]=cx*987654323+i;
}

// Update the model with bit y1, and predict next bit to mixer m.
// Context: cc=c0, bp=bpos, c1=buf(1), y1=y.
int ContextMap::mix1(Mixer& m, int cc, int c1, int y1) {

  // Update model with y
  int result=0;
  for (int i=0; i<cn; ++i) {
	U8 *cpi=cp[i];
    if (cpi) {
      assert(cpi>=&t[0].bh[0][0] && cpi<=&t[Sz].bh[6][6]);
      assert((long(cpi)&63)>=15);
      int ns=nex(*cpi, y1);
      if (ns>=204 && (rnd() << ((452-ns)>>3))) ns-=4; // probabilistic increment
      *cpi=ns;
    }

    // Update context pointers
    if (bpos>1 && runp[i][0]==0)
      cpi=0;
    else if (bpos==1||bpos==3||bpos==6)
      cpi=cp0[i]+1+(cc&1);
    else if (bpos==4||bpos==7)
      cpi=cp0[i]+3+(cc&3);
    else {
      cp0[i]=cpi=t[(cxt[i]+cc)&Sz].get(cxt[i]>>16,i);

      // Update pending bit histories for bits 2-7
      if (bpos==0) {
        if (cpi[3]==2) {
          const int c=cpi[4]+256;
          U8 *p=t[(cxt[i]+(c>>6))&Sz].get(cxt[i]>>16,i);
          p[0]=1+((c>>5)&1);
          p[p[0]        ]=1+((c>>4)&1);
          p[3+((c>>4)&3)]=1+((c>>3)&1);
          p=t[(cxt[i]+(c>>3))&Sz].get(cxt[i]>>16,i);
          p[0]=1+((c>>2)&1);
          p[p[0]        ]=1+((c>>1)&1);
          p[3+((c>>1)&3)]=1+(c&1);
          cpi[6]=0;
        }

	U8 c0 = runp[i][0];
        // Update run count of previous context
        if (c0==0)  // new context
          c0=2, runp[i][1]=c1;
        else if (runp[i][1]!=c1)  // different byte in context
          c0=1, runp[i][1]=c1;
        else if (c0<254)  // same byte in context
          c0+=2;
	runp[i][0] = c0;
        runp[i]=cpi+3;
      }
    }

    // predict from last byte in context
    int rc=runp[i][0];  // count*2, +1 if 2 different bytes seen
    if ((runp[i][1]+256)>>(8-bpos)==cc) {
      int b=((runp[i][1]>>(7-bpos))&1)*2-1;  // predicted bit + for 1, - for 0
      int c=ilog(rc+1);
	if (rc&1) c=(c*15)/4;
	     else c*=13;
      m.add(b*c);
    }
    else
      m.add(0);

    // predict from bit context
    result+=mix2(m, cpi ? *cpi : 0, sm[i]);
    cp[i]=cpi;
  }
  if (bpos==7) cn=0;
  return result;
}

//////////////////////////// Models //////////////////////////////

// All of the models below take a Mixer as a parameter and write
// predictions to it.

//////////////////////////// matchModel ///////////////////////////

// matchModel() finds the longest matching context and returns its length

int matchModel(Mixer& m) {
  const int MAXLEN=2047;  // longest allowed match + 1
  static Array<int> t(MEM);  // hash table of pointers to contexts
  static int h=0;  // hash of last 7 bytes
  static int ptr=0;  // points to next byte of match if any
  static int len=0;  // length of match, or 0 if no match
  static int result=0;

  if (!bpos) {
    h=(h*887*8+b1+1)&(t.size()-1);  // update context hash
    if (len) ++len, ++ptr;
    else {  // find match
      ptr=t[h];
      if (ptr && pos-ptr<buf.size())
        while (buf(len+1)==buf[ptr-len-1] && len<MAXLEN) ++len;
    }
    t[h]=pos;  // update hash table
    result=len;
//    if (result>0 && !(result&0xfff)) printf("pos=%d len=%d ptr=%d\n", pos, len, ptr);
  }

  // predict
  if (len>MAXLEN) len=MAXLEN;
  int sgn;
  if (len && b1==buf[ptr-1] && c0==(buf[ptr]+256)>>(8-bpos)) {
    if ((buf[ptr]>>(7-bpos))&1) sgn=8;
    else sgn=-8;
  }
  else sgn=len=0;
  m.add(sgn*ilog(len));
  m.add(sgn*8*min(len, 32));
  return result;
}

#if 0
//////////////////////////// picModel //////////////////////////

// Model a 1728 by 2376 2-color CCITT bitmap image, left to right scan,
// MSB first (216 bytes per row, 513216 bytes total).  Insert predictions
// into m.

void picModel(Mixer& m) {
  static U32 r0, r1, r2, r3;  // last 5 rows, bit 8 is over current pixel
  static Array<U8> t(0x10200);  // model: cxt -> state
  const int N=3;  // number of contexts
  static int cxt[N];  // contexts
  static StateMap sm[N];
  int i;

  // update the model
  for (i=0; i<N; ++i)
    t[cxt[i]]=nex(t[cxt[i]],y);

  // update the contexts (pixels surrounding the predicted one)
  r0+=r0+y;
  r1+=r1+((buf(215)>>(7-bpos))&1);
  r2+=r2+((buf(431)>>(7-bpos))&1);
  r3+=r3+((buf(647)>>(7-bpos))&1);
  cxt[0]=r0&0x7|r1>>4&0x38|r2>>3&0xc0;
  cxt[1]=0x100+(r0&1|r1>>4&0x3e|r2>>2&0x40|r3>>1&0x80);
  cxt[2]=0x200+(r0&0x3f^r1&0x3ffe^r2<<2&0x7f00^r3<<5&0xf800);

  // predict
  for (i=0; i<N; ++i)
    m.add(stretch(sm[i].p(t[cxt[i]])));
}
#endif

  static U32 col, frstchar=0, spafdo=0, spaces=0, spacecount=0, words=0, wordcount=0, fails=0, failz=0, failcount=0;

//////////////////////////// wordModel /////////////////////////

// Model English text (words and columns/end of line)

void wordModel(Mixer& m) {
  static U32 word0=0, word1=0, word2=0, word3=0, word4=0;  // hashes
  static ContextMap cm(MEM*64, 46);
  static int nl1=-3, nl=-2;  // previous, current newline position
  static U32 t1[256];
  static U16 t2[0x10000];

  // Update word hashes
  if (bpos==0) {
    U32 c=b1, f=0;

	if (spaces&0x80000000) --spacecount;
	if (words&0x80000000) --wordcount;
	spaces=spaces*2;
	words=words*2;

    if ( (c-'a') <= ('z'-'a') || c==8 || c==6 || (c>127&&b2!=12)) {
      ++words, ++wordcount;
      word0=word0*263*8+c;
    }
    else {
	if (c==32 || c==10) { ++spaces, ++spacecount; if (c==10) nl1=nl, nl=pos-1; }
	if (word0) {
	  word4=word3*43;
	  word3=word2*47;
	  word2=word1*53;
	  word1=word0*83;
	  word0=0;
	  if( c=='.'||c=='O'||c==('}'-'{'+'P') ) f=1, spafdo=0; else { ++spafdo; spafdo=min(63,spafdo); }
	}
    }
    
    U32 h=word0*271+c;
    cm.set(word0);
    cm.set(h+word1);
    cm.set(  word0*91+word1*89);
    cm.set(h+word1*79+word2*71);

    cm.set(h+word2);
    cm.set(h+word3);
    cm.set(h+word4);
    cm.set(h+word1*73+word3*61);
    cm.set(h+word2*67+word3*59);

	  if (f) {
	    word4=word3*31;
	    word3=word2*37;
	    word2=word1*41;
	    word1='.';
	  }

    cm.set(b3|b4<<8);
    cm.set(spafdo*8 * ((w4&3)==1) );

    col=min(31, pos-nl);
	if (col<=2) {
		if (col==2) frstchar=min(c,96); else frstchar=0;
	}
	if (frstchar=='[' && c==32)	{ if(b3==']' || b4==']' ) frstchar=96; }
    cm.set(frstchar<<11|c);

    int above=buf[nl1+col]; // text column context

    // Text column models
    cm.set(col<<16|c<<8|above);
    cm.set(col<<8|c);
    cm.set(col*(c==32));

    h = wordcount*64+spacecount;
    cm.set(spaces&0x7fff);
    cm.set(frstchar<<7);
    cm.set(spaces&0xff);
    cm.set(c*64+spacecount/2);
    cm.set((c<<13)+h);
    cm.set(        h);


    U32 d=c4&0xffff;
    h=w4<<6;
    cm.set(c+(h&0xffffff00));
    cm.set(c+(h&0x00ffff00));
    cm.set(c+(h&0x0000ff00));
    h<<=6;
    cm.set(d+(h&0xffff0000));
    cm.set(d+(h&0x00ff0000));
    h<<=6, f=c4&0xffffff;
    cm.set(f+(h&0xff000000));

    U16& r2=t2[f>>8];
    r2=r2<<8|c;
    U32& r1=t1[d>>8];
    r1=r1<<8|c;
    U32 t=c|t1[c]<<8;
    cm.set(t&0xffff);
    cm.set(t&0xffffff);
    cm.set(t);
    cm.set(t&0xff00);
    t=d|t2[d]<<16;
    cm.set(t&0xffffff);
    cm.set(t);

    cm.set(x4&0x00ff00ff);
    cm.set(x4&0xff0000ff);
    cm.set(x4&0x00ffff00);
    cm.set(c4&0xff00ff00);
    cm.set(c +b5*256+(1<<17));
    cm.set(c +b6*256+(2<<17));
    cm.set(b4+b8*256+(4<<17));

    cm.set(d);
    cm.set(w4&15);
    cm.set(f4);
    cm.set((w4&63)*128+(5<<17));
    cm.set(d<<9|frstchar);
    cm.set((f4&0xffff)<<11|frstchar);
  }
  cm.mix(m);
}

//////////////////////////// recordModel ///////////////////////

// Model 2-D data with fixed record length.  Also order 1-2 models
// that include the distance to the last match.

void recordModel(Mixer& m) {
  static int cpos1[256]; //, cpos2[256], cpos3[256], cpos4[256]; //buf(1)->last 3 pos
  static int wpos1[0x10000]; // buf(1..2) -> last position
///  static int rlen=2, rlen1=3, rlen2=4;  // run length and 2 candidates
///  static int rcount1=0, rcount2=0;  // candidate counts
  static ContextMap cm(32768/4, 2), cn(32768/2, 5), co(32768, 4), cp(32768*2, 3), cq(32768*4, 3);

  // Find record length
  if (!bpos) {
    int c=b1, w=(b2<<8)+c, d=w&0xf0ff, e=c4&0xffffff;
#if 0
    int r=pos-cpos1[c];
    if (r>1 && r==cpos1[c]-cpos2[c]
        && r==cpos2[c]-cpos3[c] && r==cpos3[c]-cpos4[c]
        && (r>15 || (c==buf(r*5+1)) && c==buf(r*6+1))) {
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
///    cm.set(rlen|buf(rlen)<<10|buf(rlen*2)<<18);
    cn.set(w    );
    cn.set(d<<8 );
    cn.set(c<<16);
    cn.set((f4&0xffff)<<3);
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
///    cpos4[c]=cpos3[c];
///    cpos3[c]=cpos2[c];
///    cpos2[c]=cpos1[c];
    cpos1[c]=pos;
    wpos1[w]=pos;
  }
  co.mix(m);
  cp.mix(m);
    cxtfl=0;
  cm.mix(m);
  cn.mix(m);
  cq.mix(m);
    cxtfl=3;
}

//////////////////////////// sparseModel ///////////////////////

// Model order 1-2 contexts with gaps.

void sparseModel(Mixer& m) {
  static ContextMap cn(MEM*2, 5);
  static SmallStationaryContextMap scm1(0x20000,17), scm2(0x20000,12), scm3(0x20000,12),
				   scm4(0x20000,13), scm5(0x10000,12), scm6(0x20000,12),
				   scm7(0x2000 ,12), scm8(0x8000 ,13), scm9(0x1000 ,12), scma(0x10000,16);

  if (bpos==0) {
    cn.set(words&0x1ffff);
    cn.set((f4&0x000fffff)*7);
    cn.set((x4&0xf8f8f8f8)+3);
    cn.set((tt&0x00000fff)*9);
    cn.set((x4&0x80f0f0ff)+6);
      scm1.set(b1);
      scm2.set(b2);
      scm3.set(b3);
      scm4.set(b4);
      scm5.set(words&127);
      scm6.set((words&12)*16+(w4&12)*4+(b1>>4));
      scm7.set(w4&15);
      scm8.set(spafdo*((w4&3)==1));
      scm9.set(col*(b1==32));
      scma.set(frstchar);
  }
  cn.mix(m);
  scm1.mix(m);
  scm2.mix(m);
  scm3.mix(m);
  scm4.mix(m);
  scm5.mix(m);
  scm6.mix(m);
  scm7.mix(m);
  scm8.mix(m);
  scm9.mix(m);
  scma.mix(m);
}

#if 0
//////////////////////////// bmpModel /////////////////////////////////

// Model a 24-bit color uncompressed .bmp or .tif file.  Return
// width in pixels if an image file is detected, else 0.

// 32-bit little endian number at buf(i)..buf(i-3)
inline U32 i4(int i) {
  assert(i>3);
  return buf(i)+256*buf(i-1)+65536*buf(i-2)+16777216*buf(i-3);
}

// 16-bit
inline int i2(int i) {
  assert(i>1);
  return buf(i)+256*buf(i-1);
}

// Square buf(i)
inline int sqrbuf(int i) {
  assert(i>0);
  return buf(i)*buf(i);
}

int bmpModel(Mixer& m) {
  static int w=0;  // width of image in bytes (pixels * 3)
  static int eoi=0;     // end of image
  static U32 tiff=0;  // offset of tif header
  const int SC=0x20000;
  static SmallStationaryContextMap scm1(SC), scm2(SC),
    scm3(SC), scm4(SC), scm5(SC), scm6(SC*2);
  static ContextMap cm(MEM*4, 8);

  // Detect .bmp file header (24 bit color, not compressed)
  if (!bpos && buf(54)=='B' && buf(53)=='M'
      && i4(44)==54 && i4(40)==40 && i4(24)==0) {
    w=(i4(36)+3&-4)*3;  // image width
    const int height=i4(32);
    eoi=pos;
    if (w<0x30000 && height<0x10000) {
      eoi=pos+w*height;  // image size in bytes
      printf("BMP %dx%d ", w/3, height);
    }
    else
      eoi=pos;
  }

  // Detect .tif file header (24 bit color, not compressed).
  // Parsing is crude, won't work with weird formats.
  if (!bpos) {
    if (c4==0x49492a00) tiff=pos;  // Intel format only
    if (pos-tiff==4 && c4!=0x08000000) tiff=0; // 8=normal offset to directory
    if (tiff && pos-tiff==200) {  // most of directory should be read by now
      int dirsize=i2(pos-tiff-4);  // number of 12-byte directory entries
      w=0;
      int bpp=0, compression=0, width=0, height=0;
      for (int i=tiff+6; i<pos-12 && --dirsize>0; i+=12) {
        int tag=i2(pos-i);  // 256=width, 257==height, 259: 1=no compression
          // 277=3 samples/pixel
        int tagfmt=i2(pos-i-2);  // 3=short, 4=long
        int taglen=i4(pos-i-4);  // number of elements in tagval
        int tagval=i4(pos-i-8);  // 1 long, 1-2 short, or points to array
        if ((tagfmt==3||tagfmt==4) && taglen==1) {
          if (tag==256) width=tagval;
          if (tag==257) height=tagval;
          if (tag==259) compression=tagval; // 1 = no compression
          if (tag==277) bpp=tagval;  // should be 3
        }
      }
      if (width>0 && height>0 && width*height>50 && compression==1
          && (bpp==1||bpp==3))
        eoi=tiff+width*height*bpp, w=width*bpp;
      if (eoi>pos)
        printf("TIFF %dx%dx%d ", width, height, bpp);
      else
        tiff=w=0;
    }
  }
  if (pos>eoi) return w=0;

  // Select nearby pixels as context
  if (!bpos) {
    assert(w>3);
    int color=pos%3;
    int mean=buf(3)+buf(w-3)+buf(w)+buf(w+3);
    const int var=sqrbuf(3)+sqrbuf(w-3)+sqrbuf(w)+sqrbuf(w+3)-mean*mean/4>>2;
    mean>>=2;
    const int logvar=ilog(var);
    int i=0;
    cm.set(hash(++i, buf(3)>>2, buf(w)>>2, color));
    cm.set(hash(++i, buf(3)>>2, buf(1)>>2, color));
    cm.set(hash(++i, buf(3)>>2, buf(2)>>2, color));
    cm.set(hash(++i, buf(w)>>2, buf(1)>>2, color));
    cm.set(hash(++i, buf(w)>>2, buf(1)>>2, color));
    cm.set(hash(++i, buf(3)+buf(w)>>1, color));
    cm.set(hash(++i, buf(3)+buf(w)>>3, buf(1)>>5, buf(2)>>5, color));
    cm.set(hash(++i, mean, logvar>>5, color));
    scm1.set(buf(3)+buf(w)>>1);
    scm2.set(buf(3)+buf(w)-buf(w+3)>>1);
    scm3.set(buf(3)*2-buf(6)>>1);
    scm4.set(buf(w)*2-buf(w*2)>>1);
    scm5.set(buf(3)+buf(w)-buf(w-3)>>1);
    scm6.set(mean>>1|logvar<<1&0x180);
  }

  // Predict next bit
  scm1.mix(m);
  scm2.mix(m);
  scm3.mix(m);
  scm4.mix(m);
  scm5.mix(m);
  scm6.mix(m);
  cm.mix(m);
  return w;
}

//////////////////////////// jpegModel /////////////////////////

// Model JPEG. Return 1 if a JPEG file is detected or else 0.
// Only the baseline and 8 bit extended Huffman coded DCT modes are
// supported.  The model partially decodes the JPEG image to provide
// context for the Huffman coded symbols.

// Print a JPEG segment at buf[p...] for debugging
void dump(const char* msg, int p) {
  printf("%s:", msg);
  int len=buf[p+2]*256+buf[p+3];
  for (int i=0; i<len+2; ++i)
    printf(" %02X", buf[p+i]);
  printf("\n");
}

// Detect invalid JPEG data.  The proper response is to silently
// fall back to a non-JPEG model.
#define jassert(x) if (!(x)) { \
/*  printf("JPEG error at %d, line %d: %s\n", pos, __LINE__, #x); */ \
  jpeg=0; \
  return next_jpeg;}

struct HUF {U32 min, max; int val;}; // Huffman decode tables
  // huf[Tc][Th][m] is the minimum, maximum+1, and pointer to codes for
  // coefficient type Tc (0=DC, 1=AC), table Th (0-3), length m+1 (m=0-15)

int jpegModel(Mixer& m) {

  // State of parser
  enum {SOF0=0xc0, SOF1, SOF2, SOF3, DHT, RST0=0xd0, SOI=0xd8, EOI, SOS, DQT,
    DNL, DRI, APP0=0xe0, COM=0xfe, FF};  // Second byte of 2 byte codes
  static int jpeg=0;  // 1 if JPEG is header detected, 2 if image data
  static int next_jpeg=0;  // updated with jpeg on next byte boundary
  static int app;  // Bytes remaining to skip in APPx or COM field
  static int sof=0, sos=0, data=0;  // pointers to buf
  static Array<int> ht(8);  // pointers to Huffman table headers
  static int htsize=0;  // number of pointers in ht

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
  static int linesize=0; // width of image in MCU
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
  static U32 huff1=0, huff2=0, huff3=0, huff4=0;  // hashes of last codes
  static int rs1, rs2, rs3, rs4;  // last 4 RS codes
  static int ssum=0, ssum1=0, ssum2=0, ssum3=0, ssum4=0;
    // sum of S in RS codes in block and last 4 values

  // Be sure to quit on a byte boundary
  if (!bpos) next_jpeg=jpeg>1;
  if (bpos && !jpeg) return next_jpeg;
  if (!bpos && app>0) --app;
  if (app>0) return next_jpeg;
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

    // Detect JPEG (SOI, APPx)
    if (!jpeg && buf(4)==FF && buf(3)==SOI && buf(2)==FF && buf(1)>>4==0xe) {
      jpeg=1;
      app=sos=sof=htsize=data=mcusize=linesize=0;
      huffcode=huffbits=huffsize=mcupos=cpos=0, rs=-1;
      memset(&huf[0], 0, huf.size()*sizeof(HUF));
      memset(&pred[0], 0, pred.size()*sizeof(int));
    }

    // Detect end of JPEG when data contains a marker other than RSTx
    // or byte stuff (00).
    if (jpeg && data && buf(2)==FF && buf(1) && (buf(1)&0xf8)!=RST0) {
      jassert(buf(1)==EOI);
      jpeg=0;
    }
    if (!jpeg) return next_jpeg;

    // Detect APPx or COM field
    if (!data && !app && buf(4)==FF && (buf(3)>>4==0xe || buf(3)==COM))
      app=buf(2)*256+buf(1)+2;

    // Save pointers to sof, ht, sos, data,
    if (buf(5)==FF && buf(4)==SOS) {
      int len=buf(3)*256+buf(2);
      if (len==6+2*buf(1) && buf(1) && buf(1)<=4)  // buf(1) is Ns
        sos=pos-5, data=sos+len+2, jpeg=2;
    }
    if (buf(4)==FF && buf(3)==DHT && htsize<8) ht[htsize++]=pos-4;
    if (buf(4)==FF && buf(3)==SOF0) sof=pos-4;

    // Restart
    if (buf(2)==FF && (buf(1)&0xf8)==RST0) {
      huffcode=huffbits=huffsize=mcupos=0, rs=-1;
      memset(&pred[0], 0, pred.size()*sizeof(int));
    }
  }

  {
	int i,j;
    // Build Huffman tables
    // huf[Tc][Th][m] = min, max+1 codes of length m, pointer to byte values
    if (pos==data && bpos==1) {
      jassert(htsize>0);
      for (i=0; i<htsize; ++i) {
        int p=ht[i]+4;  // pointer to current table after length field
        int end=p+buf[p-2]*256+buf[p-1]-2;  // end of Huffman table
        int count=0;  // sanity check
        while (p<end && end<pos && end<p+2100 && ++count<10) {
          int tc=buf[p]>>4, th=buf[p]&15;
          if (tc>=2 || th>=4) break;
          jassert(tc>=0 && tc<2 && th>=0 && th<4);
          HUF* h=&huf[tc*64+th*16]; // [tc][th][0]; 
          int val=p+17;  // pointer to values
          int hval=tc*1024+th*256;  // pointer to RS values in hbuf
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

      // Build Huffman table selection table (indexed by mcupos).
      // Get image width.
      if (!sof && sos) return next_jpeg;
      int ns=buf[sos+4];
      int nf=buf[sof+9];
      jassert(ns<=4 && nf<=4);
      mcusize=0;  // blocks per MCU
      int hmax=0;  // MCU horizontal dimension
      for (i=0; i<ns; ++i) {
        for (j=0; j<nf; ++j) {
          if (buf[sos+2*i+5]==buf[sof+3*j+10]) { // Cs == C ?
            int hv=buf[sof+3*j+11];  // packed dimensions H x V
            if (hv>>4>hmax) hmax=hv>>4;
            hv=(hv&15)*(hv>>4);  // number of blocks in component C
            jassert(hv>=1 && hv+mcusize<=10);
            while (hv) {
              jassert(mcusize<10);
              hufsel[0][mcusize]=buf[sos+2*i+6]>>4&15;
              hufsel[1][mcusize]=buf[sos+2*i+6]&15;
              jassert (hufsel[0][mcusize]<4 && hufsel[1][mcusize]<4);
              color[mcusize]=i;
              --hv;
              ++mcusize;
            }
          }
        }
      }
      jassert(hmax>=1 && hmax<=10);
      width=buf[sof+7]*256+buf[sof+8];  // in pixels
      int height=buf[sof+5]*256+buf[sof+6];
      printf("JPEG %dx%d ", width, height);
      width=(width-1)/(hmax*8)+1;  // in MCU
      jassert(width>0);
      mcusize*=64;  // coefficients per MCU
      row=column=0;
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
          huff4=huff3;
          huff3=huff2;
          huff2=huff1;
          huff1=hash(huffcode, huffbits);
          rs4=rs3;
          rs3=rs2;
          rs2=rs1;
          rs1=rs;
          int x=0;  // decoded extra bits
          if (mcupos&63) {  // AC
            if (rs==0) { // EOB
              mcupos=mcupos+63&-64;
              jassert(mcupos>=0 && mcupos<=mcusize && mcupos<=640);
              while (cpos&63) cbuf[cpos++]=0;
            }
            else {  // rs = r zeros + s extra bits for the next nonzero value
                    // If first extra bit is 0 then value is negative.
              jassert((rs&15)<=10);
              const int r=rs>>4;
              const int s=rs&15;
              jassert(mcupos>>6==mcupos+r>>6);
              mcupos+=r+1;
              x=huffcode&(1<<s)-1;
              if (s && !(x>>s-1)) x-=(1<<s)-1;
              for (int i=r; i>=1; --i) cbuf[cpos++]=i<<4|s;
              cbuf[cpos++]=s<<4|huffcode<<2>>s&3|12;
              ssum+=s;
            }
          }
          else {  // DC: rs = 0S, s<12
            jassert(rs<12);
            ++mcupos;
            x=huffcode&(1<<rs)-1;
            if (rs && !(x>>rs-1)) x-=(1<<rs)-1;
            jassert(mcupos>=0 && mcupos>>6<10);
            const int comp=color[mcupos>>6];
            jassert(comp>=0 && comp<4);
            dc=pred[comp]+=x;
            jassert((cpos&63)==0);
            cbuf[cpos++]=dc+1023>>3;
            ssum4=ssum3;
            ssum3=ssum2;
            ssum2=ssum1;
            ssum1=ssum;
            ssum=rs;
          }
          jassert(mcupos>=0 && mcupos<=mcusize);
          if (mcupos>=mcusize) {
            mcupos=0;
            if (++column==width) column=0, ++row;
          }
          huffcode=huffsize=huffbits=0, rs=-1;
        }
      }
    }
  }

  // Estimate next bit probability
  if (!jpeg || !data) return next_jpeg;

  // Context model
  const int N=19;  // size of t, number of contexts
  static BH<9> t(MEM);  // context hash -> bit history
    // As a cache optimization, the context does not include the last 1-2
    // bits of huffcode if the length (huffbits) is not a multiple of 3.
    // The 7 mapped values are for context+{"", 0, 00, 01, 1, 10, 11}.
  static Array<U32> cxt(N);  // context hashes
  static Array<U8*> cp(N);  // context pointers
  static StateMap sm[N];
  static Mixer m1(32, 800, 4);
  static APM a1(1024), a2(0x10000);
  const static U8 zzu[64]={  // zigzag coef -> u,v
    0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,
    3,2,1,0,1,2,3,4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};
  const static U8 zzv[64]={
    0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,
    4,5,6,7,7,6,5,4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};


  // Update model
  if (cp[N-1]) {
    for (int i=0; i<N; ++i)
      *cp[i]=nex(*cp[i],y);
  }
  m1.update();

  // Update context
  const int comp=color[mcupos>>6];
  const int coef=(mcupos&63)|comp<<6;
  const int hc=huffcode|1<<huffbits;
  static int hbcount=2;
  if (++hbcount>2 || huffbits==0) hbcount=0;
  jassert(coef>=0 && coef<256);
  const int zu=zzu[mcupos&63], zv=zzv[mcupos&63];
  if (hbcount==0) {
    const int mpos=mcupos>>4|!(mcupos&-64)<<7;
    int n=0;
    cxt[0]=hash(++n, hc, mcupos>>2, min(3, mcupos&63));
    cxt[1]=hash(++n, hc, mpos>>4, cbuf[cpos-mcusize]);
    cxt[2]=hash(++n, hc, mpos>>4, cbuf[cpos-width*mcusize]);
    cxt[3]=hash(++n, hc, ilog(ssum3), coef);
    cxt[4]=hash(++n, hc, coef, column>>3);
    cxt[5]=hash(++n, hc, coef, column>>1);
    cxt[6]=hash(++n, hc, rs1, mpos);
    cxt[7]=hash(++n, hc, rs1, rs2);
    cxt[8]=hash(++n, hc, rs1, rs2, rs3);
    cxt[9]=hash(++n, hc, ssum>>4, mcupos);
    cxt[10]=hash(++n, hc, mpos, cbuf[cpos-1]);
    cxt[11]=hash(++n, hc, dc);
    cxt[12]=hash(++n, hc, rs1, coef);
    cxt[13]=hash(++n, hc, rs1, rs2, coef);
    cxt[14]=hash(++n, hc, mcupos>>3, ssum3>>3);
    cxt[15]=hash(++n, hc, huff1);
    cxt[16]=hash(++n, hc, coef, huff1);
    cxt[17]=hash(++n, hc, zu, comp);
    cxt[18]=hash(++n, hc, zv, comp);
  }

  // Predict next bit
  m1.add(128);
  assert(hbcount<=2);
  for (int i=0; i<N; ++i) {
    if (hbcount==0) cp[i]=t[cxt[i]]+1;
    else if (hbcount==1) cp[i]+=1+(huffcode&1)*3;
    else cp[i]+=1+(huffcode&1);
    int sp=stretch(sm[i].p(*cp[i]));
    m1.add(sp);
  }
  m1.set(0, 1);
  m1.set(coef, 64);
  m1.set(mcupos, 640);
  int pr=m1.p();
  pr=a1.p(pr, hc&1023);
  pr=a2.p(pr, hc&255|coef<<8);
  m.add(stretch(pr));
  return 1;
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
  static ContextMap cm(MEM*2, N);
  if (!bpos) {
    for (int i=0; i<N; ++i)
      cm.set(execxt(i, buf(1)*(i>4)));
  }
  cm.mix(m);
}
//#endif

//////////////////////////// indirectModel /////////////////////

// The context is a byte string history that occurs within a
// 1 or 2 byte context.

void indirectModel(Mixer& m) {
  static ContextMap cm(MEM*8, 5);
  static U32 t1[256];
  static U16 t2[0x10000];

  if (!bpos) {
    U32 d=c4&0xffff, c=d&255;
    U32& r1=t1[d>>8];
    r1=r1<<8|c;
    U16& r2=t2[c4>>8&0xffff];
    r2=r2<<8|c;
    U32 t=c|t1[c]<<8;
    cm.set(t&0xffff);
    cm.set(t&0xffffff);
    cm.set(t);
    t=d|t2[d]<<16;
    cm.set(t&0xffffff);
    cm.set(t);
  }
  cm.mix(m);
}
#endif

int primes[]={ 0, 257,251,241,239,233,229,227,223,211,199,197,193,191 };
static U32 WRT_mpw[16]= { 3, 3, 3, 2, 2, 2, 1, 1,  1, 1, 1, 1, 1, 0, 0, 0 }, tri[4]={0,4,3,7}, trj[4]={0,6,6,12};
static U32 WRT_mtt[16]= { 0, 0, 1, 2, 3, 4, 5, 5,  6, 6, 6, 6, 6, 7, 7, 7 };

//////////////////////////// contextModel //////////////////////

// file types (order is important: the last will be sorted by filetype detection as the first)

// This combines all the context models with a Mixer.

int contextModel2() {
  static ContextMap cm(MEM*32, 7);
  static RunContextMap rcm7(MEM/4,14), rcm9(MEM/4,18), rcm10(MEM/2,20);
  static Mixer m(456, 128*(16+14+14+12+14+16), 6, 512);
  static U32 cxt[16];  // order 0-11 contexts
  //static Filetype filetype=DEFAULT;
  static int size=0;  // bytes remaining in block
//  static const char* typenames[4]={"", "jpeg ", "exe ", "text "};

  // Parse filetype and size
  if (bpos==0) {
    --size;
    //if (size==-1) filetype=(Filetype)b1;
    if (size==-5) {
      size=c4;
//      if (filetype<=3) printf("(%s%d)", typenames[filetype], size);
/////      if (filetype==EXE) size+=8;
    }
  }

  m.update();
  m.add(64);

  // Test for special file types
/////  int isjpeg=jpegModel(m);  // 1 if JPEG is detected, else 0
 // int ismatch=matchModel(m);  // Length of longest matching context
/////  int isbmp=bmpModel(m);  // Image width (bytes) if BMP or TIFF detected, or 0

  //if (ismatch>1024) {   // Model long matches directly
  //  m.set(0, 8);
  //  return m.p();
  //}
/////  else if (isjpeg) {
/////    m.set(1, 8);
/////    m.set(c0, 256);
/////    m.set(buf(1), 256);
/////    return m.p();
/////  }
/////  else if (isbmp>0) {
/////    static int col=0;
/////    if (++col>=24) col=0;
/////    m.set(2, 8);
/////    m.set(col, 24);
/////    m.set(buf(isbmp)+buf(3)>>4, 32);
/////    m.set(c0, 256);
/////    return m.p();
/////  }

  // Normal model
  if (bpos==0) {
    int i=0, f2=buf(2);

    if(f2=='.'||f2=='O'||f2=='M'||f2=='!'||f2==')'||f2==('}'-'{'+'P')) {
      if (b1!=(unsigned int)f2 && buf(3)!=f2 ) i=13, x4=x4*256+f2;
    }

    for (; i>0; --i)  // update order 0-11 context hashes
      cxt[i]=cxt[i-1]*primes[i];
    
    for (i=13; i>0; --i)  // update order 0-11 context hashes
      cxt[i]=cxt[i-1]*primes[i]+b1;

    cm.set(cxt[3]);
    cm.set(cxt[4]);
    cm.set(cxt[5]);
    cm.set(cxt[6]);
    cm.set(cxt[8]);
    cm.set(cxt[13]);
    cm.set(0);

    rcm7.set(cxt[7]);
    rcm9.set(cxt[9]);
    rcm10.set(cxt[11]);

	x4=x4*256+b1;
  }
  rcm7.mix(m);
  rcm9.mix(m);
  rcm10.mix(m);
    int qq=m.nx;
  order=cm.mix(m)-1;
  if(order<0) order=0;
    int zz=(m.nx-qq)/7;

  m.nx=qq+zz*3;
    for (qq=zz*2;qq!=0;--qq) m.mul(5);
    for (qq=zz; qq!=0; --qq) m.mul(6);
    for (qq=zz; qq!=0; --qq) m.mul(9);

  if (level>=4) {
    wordModel(m);
    sparseModel(m);
    recordModel(m);
    //indirectModel(m);
/////    if (filetype==EXE) exeModel(m);
/////    if (fsize==513216) picModel(m);
  }

//m.tx[420]=0;
//m.tx[425]=0;
//m.tx[430]=0;
		U32 c1=b1, c2=b2, c;
		if (c1==9 || c1==10 || c1==32) c1=16;
		if (c2==9 || c2==10 || c2==32) c2=16;

		m.set(256*order + (w4&240) + (c2>>4), 256*7);

		c=(words>>1)&63;
		m.set((w4&3)*64+c+order*256, 256*7);

		c=(w4&255)+256*bpos;
		m.set(c, 256*8);

		if(bpos)
		{
		c=c0<<(8-bpos);  if(bpos==1)c+=b3/2;
		c=(min(bpos,5))*256+(tt&63)+(c&192);
		}
		else c=(words&12)*16+(tt&63);
		m.set(c, 1536);

		c=bpos; c2=(c0<<(8-bpos)) | (c1>>bpos);
		m.set(order*256 + c + (c2&248), 256*7);

		c=c*256+((c0<<(8-bpos))&255);
		c1 = (words<<bpos) & 255;
		m.set(c+(c1>>bpos), 2048);

  return m.p();
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
  static APM a1(256), a2(0x8000), a3(0x8000), a4(0x20000), a5(0x10000), a6(0x10000);

  // Update global context: pos, bpos, c0, c4, buf
  c0+=c0+y;
  if (c0>=256) {
    buf[pos++]=c0;
    c0-=256;
	if (pos<=1024*1024) {
		if (pos==1024*1024) sm_shft=9, sm_add=65535+511;
		if (pos== 512*1024) sm_shft=8, sm_add=65535+255;
		sm_add_y = sm_add & (-y);
	}
	int i=WRT_mpw[c0>>4];
	w4=w4*4+i;
	if (b1==12) i=2;
	w5=w5*4+i;
	b8=b7, b7=b6, b6=b5, b5=b4,
	b4=b3; b3=b2; b2=b1; b1=c0;
	if(c0=='.' || c0=='O' || c0=='M' || c0=='!' || c0==')' || c0==('}'-'{'+'P')) {
					w5=(w5<<8)|0x3ff, x5=(x5<<8)+c0, f4=(f4&0xfffffff0)+2;
					if(c0!='!' && c0!='O') w4|=12;
					if(c0!='!') b2='.', tt=(tt&0xfffffff8)+1;
                                }
    c4=(c4<<8)+c0;
	x5=(x5<<8)+c0;
	if (c0==32) --c0;
	f4=f4*16+(c0>>4);
	tt=tt*8+WRT_mtt[c0>>4];
    c0=1;
  }
  bpos=(bpos+1)&7;

        if (fails&0x00000080) --failcount;
        fails=fails*2;
        failz=failz*2;
        if (y) pr^=4095;
        if (pr>=1820) ++fails, ++failcount;
        if (pr>= 848) ++failz;

  // Filter the context model with APMs
  pr=contextModel2();

  int rate=6 + (pos>14*256*1024) + (pos>28*512*1024);
  int pt, pu=(a1.p(pr, c0, 3)+7*pr+4)>>3, pv, pz=failcount+1;
        pz+=tri[(fails>>5)&3];
        pz+=trj[(fails>>3)&3];
        pz+=trj[(fails>>1)&3];
        if (fails&1) pz+=8;
        pz=pz/2;

  pu=a4.p(pu,   (c0*2)^(hash(b1, (x5>>8)&255, (x5>>16)&0x80ff)&0x1ffff), rate);
  pv=a2.p(pr,   (c0*8)^(hash(29,failz&2047)&0x7fff), rate+1);
  pv=a5.p(pv,          hash(c0,w5&0xfffff)&0xffff, rate);
  pt=a3.p(pr,  (c0*32)^(hash(19,     x5&0x80ffff)&0x7fff), rate);
  pz=a6.p(pu,   (c0*4)^(hash(min(9,pz),x5&0x80ff)&0xffff), rate);

  if (fails&255)   pr =(pt*6+pu  +pv*11+pz*14 +16)>>5;
  else		   pr =(pt*4+pu*5+pv*12+pz*11 +16)>>5;
}

Predictor paq8;

void PAQ8Init() {
  buf.setsize(MEM*8);
}

int PAQ8Predict() {
  return paq8.p();
}

void PAQ8Perceive(int bit) {
  y = bit;
  paq8.update();
}

}  // namespace

PAQ8HP::PAQ8HP() {
  PAQ8Init();
}

float PAQ8HP::Predict() {
  int p = PAQ8Predict();
  return (1.0 + p) / 4097;
}

void PAQ8HP::Perceive(int bit) {
  PAQ8Perceive(bit);
}
