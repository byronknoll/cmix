// This SSE code was written by Eugene Shelwien.
// https://encode.ru/threads/2515-mod_ppmd

#include "sse.h"

#include <math.h>

namespace SSE_sh {

typedef unsigned short word;
typedef unsigned int uint;
typedef unsigned char byte;

struct SSEi_updstr {
  int P;
  int sw;
  word *C1;
  word *T;
};

template<int SSEQuant=7, int SCALElog=15, int InitFlag=0>
struct SSEi {

  static const int SCALE = 1<<SCALElog;
  static const int mSCALE = SCALE-1;

  word P[SSEQuant];

  void Init( int Wi=0 ) {
    int i,p1;

    int SCw = (SCALE-Wi)/(SSEQuant-1);
    int INC = Wi/2 + 8192;
    for( i=0,p1=INC; i<SSEQuant; i++,p1+=SCw ) P[i] = p1 ;
  }

  int SSE_Pred( int iP, SSEi_updstr& X ) {

    int sseFreq = ((SSEQuant-1)*iP)>>SCALElog;

    X.sw = ((SSEQuant-1)*iP)&mSCALE;
    X.C1 = &P[sseFreq+0];
    int f = (((SCALE-X.sw)*X.C1[0]+X.sw*X.C1[1])>>SCALElog) - 8192;

    if( f<=0 ) f=1;
    if( f>=SCALE ) f=mSCALE;

    X.P = f;

    return f;
  }

  void SSE_Update( byte c, int wr0, SSEi_updstr& X ) {

    X.P = X.P*(SCALE-wr0)>>SCALElog;
    if( c==0 ) X.P+=wr0;

    int dC = (X.C1[0]-X.C1[1]);
    int sw_dC= ((X.sw*dC+mSCALE)>>SCALElog);
    X.C1[0] = X.P + sw_dC +8192;
    X.C1[1] = X.P - (dC-sw_dC) +8192;
  }

};


enum {
  SCALElog = 15,
  SCALE = 1<<SCALElog,
  hSCALE = SCALE/2,
  mSCALE = SCALE-1
};

#ifndef M_LOG2E
#define M_LOG2E 1.44269504088896340736
#endif

double log2( double a ) { return M_LOG2E*log(a); }

double exp2( double a ) { return exp( a/M_LOG2E ); }

double st( double p ) { return log2((1-p)/p); }

double sq( double p ) { return 1.0/(1.0+exp2(p)); }

const double st_coef = (hSCALE-1)/log2(SCALE-1);
const double sq_coef = 1.0 / st_coef;

uint st_i( uint p ) {
  p = st(double(p)/SCALE) * st_coef + hSCALE;
  return p;
}

uint sq_i( uint p ) {
  p = sq( double(int(p)-hSCALE) * sq_coef ) * SCALE;
  return p;
}

double st_d( double p ) {
  p = st(p/SCALE);
  return p;
}

double sq_d( double p ) {
  p = sq(p) * SCALE;
  return p;
}

word t_st[SCALE];
word t_sq[SCALE];

uint Init_ST_SQ( void ) {
  uint i,s,x,y;

  for( i=1; i<SCALE; i++ ) t_sq[i] = sq_i(i);

  x=0; t_st[x]=0;
  for( i=1; i<SCALE; i++ ) {
    s = st_i(i);
    t_st[i] = s;

    if( s!=t_st[x] ) {
      y = i-1;
      t_sq[ t_st[x] ] = (x+y+1)/2;
      x = i;
    }
  }

  for( i=1; i<SCALE; i++ ) {
    s = t_st[i];
  }
  return 0;
}

uint unused_0 = Init_ST_SQ();


uint Extrap( int p1, int C ) {
  p1 = (((p1-hSCALE)*C) >> 13) + hSCALE;
  if( p1<1 ) p1=1;
  if( p1>mSCALE ) p1=mSCALE;
  return p1;
}

uint WExtrap( int _p1, int C ) {
  double p1 = st_d(_p1-hSCALE);
  p1 = (p1*C)/8192;
  _p1 = sq_d(p1);
  if( _p1<0 ) _p1=0;
  if( _p1>SCALE ) _p1=SCALE;
  return _p1+hSCALE;
}

struct Mixer {
  int w;

  void Init( int w0 ) {
    w = w0 + hSCALE;
  }

  int rdiv( int x, int a, int d ) {
    return x>=0 ? (x+a)>>d : -((-x+a)>>d);
  }


  int Mixup( int w, int s1, int s0 ) {
    int x = s1 + rdiv( (w-hSCALE)*(s0-s1), 1<<(SCALElog-1), SCALElog );
    x = (x>0) ? (x<SCALE) ? x : SCALE-1 : 1;
    return x;
  }

  void Update( int y, int p0,int p1, int wq, int pm ) {
    int py = SCALE - (y<<SCALElog);
    int e = (py-pm);
    int d = rdiv( e*(p0-p1), 1<<(SCALElog-1), SCALElog );
    d = rdiv( d*wq, 1<<(SCALElog-1), SCALElog );
    w += d;
  }
};

static const byte M_mx1mask0[256]={ 0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41,42,42,43,43,44,44,45,45,46,46,47,47,47,47,48,48,48,48,49,49,49,49,50,50,50,50,51,51,51,51,52,52,52,52,53,53,53,53,54,54,54,54,55,55,55,55,56,56,56,56,57,57,57,57,58,58,58,58,59,59,59,59,60,60,60,60,61,61,61,61,62,62,62,62,63,63,63,63,63,63,63,63,64,64,64,64,64,64,64,64,65,65,65,65,65,65,65,65,66,66,66,66,66,66,66,66,67,67,67,67,67,67,67,67,68,68,68,68,68,68,68,68,69,69,69,69,69,69,69,69,70,70,70,70,70,70,70,70,71,71,71,71,71,71,71,71,72,72,72,72,72,72,72,72,73,73,73,73,73,73,73,73,74,74,74,74,74,74,74,74,75,75,75,75,75,75,75,75,76,76,76,76,76,76,76,76,77,77,77,77,77,77,77,77,78,78,78,78,78,78,78,78 };



static const byte M_sm7mask0[256]={ 0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254 };


static const int M_f0C = (10240+0) * (1);
static const int M_f1C = (7935+0) * (1);
static const int M_f2C = (9592+0) * (1);
static const int M_sm6wrA = (0+0) * (1);
static const int M_sm6wrB = (106+0) * (1);
static const int M_sm6mw = (0+0) * (1);
static const int M_sm6C1 = (8092+0) * (1);
static const int M_x1W0 = (7648+1) * (1);
static const int M_x1wr = (6202+0) * (1);
static const int M_f3C = (8200+0) * (1);
static const int M_f4C = (7677+0) * (1);
static const int M_sm7wrB = (127+0) * (1);
static const int M_sm7mw = (8192+0) * (1);
static const int M_sm7C1 = (8202+0) * (1);
static const int M_x2W0 = (2560+1) * (1);
static const int M_x2wr = (8320+0) * (1);

static const int M_mix1_Volume = 1* 4* (1<<8)* (1<<3)* 79;
static const int M_mix2_Volume = 1* 3* (1<<1)* (1<<8)* 256;
static const int M_sm6x_Volume = 1* 3* (1<<7)* (1<<8)* 256;
static const int M_sm7x_Volume = 1* 3* (1<<5)* (1<<8)* 255;

struct M_T {

  SSEi_updstr su6,su7;
  int sm6x,mix1, sm7x,mix2;
  uint mix1_s0,mix1_s1,mix1_p, mix1_p1,mix1_p2;
  uint mix2_s0,mix2_s1,mix2_p;
  uint M_j, M_pc, M_ffl;


  SSEi<7> s6[M_sm6x_Volume];
  Mixer x1[M_mix1_Volume];
  SSEi<7> s7[M_sm7x_Volume];
  Mixer x2[M_mix2_Volume];

  void M_Init( void ) {

    uint i;

    for( i=0; i<M_sm6x_Volume; i++ ) s6[i].Init(M_sm6mw);
    for( i=0; i<M_mix1_Volume; i++ ) x1[i].Init(M_x1W0);

    for( i=0; i<M_sm7x_Volume; i++ ) s7[i].Init(M_sm7mw);
    for( i=0; i<M_mix2_Volume; i++ ) x2[i].Init(M_x2W0);

    M_j=1; M_pc=0; M_ffl=0;

  }
  void M_Quit( void ) {
  }
};

struct M_T1 : M_T {

uint M_Estimate( uint p ) {
  uint p0,p1,s0,s1,s2;
  uint p2,s4,s5;
  uint j=M_j, pc=M_pc, ffl=M_ffl, prq=p>>11;

  sm7x = 0;
  sm7x = sm7x*3 + ((prq>0+(1-1))+(prq>14+(1-1)));
  sm7x = (sm7x<<5) + ((ffl)&31);
  sm7x = (sm7x<<8) + ((pc)&255);
  sm7x = (sm7x*255) + M_sm7mask0[j];

  mix2 = 0;
  mix2 = mix2*3 + ((prq>0+(1-1))+(prq>14+(1-1)));
  mix2 = (mix2<<1) + ((ffl)&1);
  mix2 = (mix2<<8) + ((pc)&255);
  mix2 = (mix2*256) + (j);

  sm6x = 0;
  sm6x = sm6x*3 + ((prq>0+(1-1))+(prq>14+(1-1)));
  sm6x = (sm6x<<7) + ((ffl)&127);
  sm6x = (sm6x<<8) + ((pc)&255);
  sm6x = (sm6x*256) + (j);

  mix1 = 0;
  mix1 = mix1*4 + ((prq>0+(1-1))+(prq>7+(1-1))+(prq>14+(1-1)));
  mix1 = (mix1<<8) + ((ffl)&255);
  mix1 = (mix1<<3) + (((pc)>>5)&7);
  mix1 = (mix1*79) + M_mx1mask0[j];

  p0 = p;
  p1 = s6[sm6x].SSE_Pred( t_sq[Extrap(t_st[p0],M_f0C)], su6 );
  s0 = t_st[p0]; s0 = Extrap(s0,M_f1C);
  s1 = t_st[p1]; s1 = Extrap(s1,M_f2C);
  mix1_s0=s0; mix1_s1=s1;
  s2 = x1[mix1].Mixup( x1[mix1].w, mix1_s0, mix1_s1 ); s2 = Extrap(s2,M_sm6C1);
  mix1_p = t_sq[s2];

  p2 = s7[sm7x].SSE_Pred( t_sq[Extrap(t_st[p0],M_f3C)], su7 );
  s4 = t_st[p2]; s4 = Extrap(s4,M_f4C);
  mix2_s0=s2; mix2_s1=s4;
  s5 = x2[mix2].Mixup( x2[mix2].w, mix2_s0, mix2_s1 ); s5 = Extrap(s5,M_sm7C1);
  mix2_p = t_sq[s5];
  mix2_s0=s2;
  mix2_s1=s4;

  return mix2_p;
}

void M_Update( uint bit ) {
  s6[sm6x].SSE_Update( bit, M_sm6wrA*128+M_sm6wrB, su6 );
  x1[mix1].Update( bit, mix1_s0,mix1_s1, M_x1wr, mix1_p );

  s7[sm7x].SSE_Update( bit, M_sm6wrA*128+M_sm7wrB, su7 );
  x2[mix2].Update( bit, mix2_s0,mix2_s1, M_x2wr, mix2_p );

  M_j += M_j+bit;

  if( M_j>=256 ) {
    M_ffl= byte(M_ffl*2+(M_pc>=0x40));
    M_pc = byte(M_j);
    M_j = 1;
  }
}
};

}

using SSE_sh::M_T1;

SSE::SSE() : sse_(new M_T1()) {
  sse_->M_Init();
}

SSE::~SSE() {
  delete sse_;
}

float SSE::Predict(float input) {
  int discrete = 1 + (1 - input) * 32766;
  int estimate = sse_->M_Estimate(discrete);
  return 1 - ((estimate - 1) / 32766.0);
}

void SSE::Perceive(int bit) {
  sse_->M_Update(bit);
}
