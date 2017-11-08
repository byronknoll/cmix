// This preprocessor is adapted from paq8l, paq8hp12any and paq8px_v101.

#include <vector>
#include <cstdlib>
#include <string.h>

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

#include "textfilter.cpp"
#include "preprocessor.h"

namespace preprocessor {

typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;

inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}

bool IsAscii(int byte) {
  if (byte >= 9 && byte <= 13) return true;
  if (byte >= 32 && byte <= 126) return true;
  return false;
}

#define bswap(x) \
+   ((((x) & 0xff000000) >> 24) | \
+    (((x) & 0x00ff0000) >>  8) | \
+    (((x) & 0x0000ff00) <<  8) | \
+    (((x) & 0x000000ff) << 24))

#define IMG_DET_NOHDR(type,start_pos,width,height) return detd=(width)*(height),info=(width),fseek(in, start+(start_pos), SEEK_SET),(type)

#define IMG_DET(type,start_pos,header_len,width,height) return dett=(type),deth=(header_len),detd=(width)*(height),info=(width),fseek(in, start+(start_pos), SEEK_SET),HDR

int info;

void Pretrain(Predictor* p, FILE* dictionary) {
  fseek(dictionary, 0L, SEEK_END);
  unsigned int len = ftell(dictionary);
  fseek(dictionary, 0L, SEEK_SET);

  std::vector<unsigned char> header;
  header.push_back(DEFAULT);
  header.push_back(len>>24);
  header.push_back(len>>16);
  header.push_back(len>>8);
  header.push_back(len);

  for (unsigned int i = 0; i < header.size(); ++i) {
    for (int j = 7; j >= 0; --j) {
      p->Pretrain((header[i]>>j)&1);
    }
  }

  unsigned int percent = 1 + (len / 100);
  for (unsigned int i = 0; i < len; ++i) {
    unsigned char c = getc(dictionary);
    if (c == '\n') c = ' ';
    if (i % percent == 0) {
      printf("\rpretraining: %d%%", i / percent);
      fflush(stdout);
    }
    for (int j = 7; j >= 0; --j) {
      p->Pretrain((c>>j)&1);
    }
  }
  printf("\r                 \r");
  fflush(stdout);
}

Filetype detect(FILE* in, int n, Filetype type) {
  U32 buf1=0, buf0=0;
  long start=ftell(in);

  // For EXE detection
  std::vector<int> abspos(256, 0),
    relpos(256, 0);
  int e8e9count=0;
  int e8e9pos=0;
  int e8e9last=0;

  // For JPEG detection
  int soi=0, sof=0, sos=0, app=0;

  // For TEXT detection
  int ascii_start = -1;
  int ascii_run = 0;
  int space_count = 0;

  // For BMP detection
  uint64_t bmp=0;
  int imgbpp=0,bmpx=0,bmpy=0,bmpof=0;
  static int deth=0,detd=0;  // detected header/data size in bytes
  static Filetype dett;  // detected block type
  if (deth) return fseek(in, start+deth, SEEK_SET),deth=0,dett;
  else if (detd) return fseek(in, start+detd, SEEK_SET),detd=0,DEFAULT;
  // For TGA detection
  uint64_t tga=0;
  int tgaid=0, tgaw=0, tgah=0;

  for (int i=0; i<n; ++i) {
    int c=getc(in);
    if (c==EOF) return (Filetype)(-1);
    buf1=buf1<<8|buf0>>24;
    buf0=buf0<<8|c;

    if (!soi && i>=3 && (buf0&0xffffff00)==0xffd8ff00 && ((buf0&0xFE)==0xC0 || (U8)buf0==0xC4 || ((U8)buf0>=0xDB && (U8)buf0<=0xFE) )) soi=i, app=i+2, sos=sof=0;
    if (soi) {
      if (app==i && (buf0>>24)==0xff &&
         ((buf0>>16)&0xff)>0xc1 && ((buf0>>16)&0xff)<0xff) app=i+(buf0&0xffff)+2;
      if (app<i && (buf1&0xff)==0xff && (buf0&0xfe0000ff)==0xc0000008) sof=i;
      if (sof && sof>soi && i-sof<0x1000 && (buf0&0xffff)==0xffda) {
        sos=i;
        if (type!=JPEG) return fseek(in, start+soi-3, SEEK_SET), JPEG;
      }
      if (i-soi>0x40000 && !sos) soi=0;
    }
    if (type==JPEG && sos && i>sos && (buf0&0xff00)==0xff00
        && (buf0&0xff)!=0 && (buf0&0xf8)!=0xd0) return DEFAULT;

    if (((buf1&0xfe)==0xe8 || (buf1&0xfff0)==0x0f80) && ((buf0+1)&0xfe)==0) {
      int r=buf0>>24;
      int a=((buf0>>24)+i)&0xff;
      int rdist=i-relpos[r];
      int adist=i-abspos[a];
      if (adist<rdist && adist<0x800 && abspos[a]>5) {
        e8e9last=i;
        ++e8e9count;
        if (e8e9pos==0 || e8e9pos>abspos[a]) e8e9pos=abspos[a];
      }
      else e8e9count=0;
      if (type!=EXE && e8e9count>=4 && e8e9pos>5)
        return fseek(in, start+e8e9pos-5, SEEK_SET), EXE;
      abspos[a]=i;
      relpos[r]=i;
    }
    if (type==EXE && i-e8e9last>0x4000)
      return fseek(in, start+e8e9last, SEEK_SET), DEFAULT;

    // Detect TEXT
    if (type == DEFAULT) {
      if (IsAscii(c)) {
        if (ascii_start == -1) {
          ascii_start = i;
          ascii_run = 0;
          space_count = 0;
        }
        if (c == ' ') ++space_count;
        ++ascii_run;
        if (ascii_run > 500) {
          if (space_count < 5) {
            ascii_start = -1;
          } else {
            return fseek(in, start + ascii_start, SEEK_SET), TEXT;
          }
        }
      } else {
        ascii_start = -1;
      }
    } else if (type == TEXT) {
      if (IsAscii(c)) {
        ascii_run -= 2;
        if (ascii_run < 0) ascii_run = 0;
      } else {
        ascii_run += 3;
        if (ascii_run > 300) {
          return fseek(in, ftell(in) - 100, SEEK_SET), DEFAULT;
        }
      }
    }

    // Detect .bmp image
    if ((buf0&0xffff)==16973) imgbpp=bmpx=bmpy=bmpof=0,bmp=i;
    if (bmp) {
      const int p=int(i-bmp);
      if (p==12) bmpof=bswap(buf0);
      else if (p==16 && buf0!=0x28000000) bmp=0;
      else if (p==20) bmpx=bswap(buf0),bmp=((bmpx==0||bmpx>0x40000)?0:bmp);
      else if (p==24) bmpy=abs((int)bswap(buf0)),bmp=((bmpy==0||bmpy>0x20000)?0:bmp);
      else if (p==27) imgbpp=c,bmp=((imgbpp!=1 && imgbpp!=4 && imgbpp!=8 && imgbpp!=24)?0:bmp);
      else if (p==31) {
        if (imgbpp!=0 && buf0==0 && bmpx>1) {
          if (imgbpp==24) {
            int width = ((bmpx*3)+3)&-4;
            IMG_DET(IMAGE24, (bmp-1), bmpof, width, bmpy);
          }
        }
        bmp=0;
      }
    }
    // detect .tga image
    if ((buf1&0xFFFFFF)==0x000200 && buf0==0x00000000) tga=i, tgaid=buf1>>24;
    if (tga){
      if ((i-tga)==8){
        tgaw = bswap(buf0)&0xffff;
        tgah = bswap(buf0)>>16;
        tga*=(!buf1 && tgaw && tgaw<0x4000 && tgah && tgah<0x4000);
      }
      else if ((i-tga)==10){
        if (tga && (buf0&0xFFDF)==0x1800)
          IMG_DET(IMAGE24, (tga-7), 18+tgaid, tgaw*3, tgah);
      }
    }
    
    // Detect .tiff image
    if (buf1==0x49492a00 && n>i+(int)bswap(buf0)) {
      long savedpos=ftell(in);
      fseek(in, start+i+bswap(buf0)-7, SEEK_SET);

      // read directory
      int dirsize=getc(in);
      int tifx=0,tify=0,tifz=0,tifzb=0,tifc=0,tifofs=0,tifofval=0,b[12];
      if (getc(in)==0) {
        for (int i=0; i<dirsize; i++) {
          for (int j=0; j<12; j++) b[j]=getc(in);
          if (b[11]==EOF) break;
          int tag=b[0]+(b[1]<<8);
          int tagfmt=b[2]+(b[3]<<8);
          int taglen=b[4]+(b[5]<<8)+(b[6]<<16)+(b[7]<<24);
          int tagval=b[8]+(b[9]<<8)+(b[10]<<16)+(b[11]<<24);
          if (tagfmt==3||tagfmt==4) {
            if (tag==256) tifx=tagval;
            else if (tag==257) tify=tagval;
            else if (tag==258) tifzb=taglen==1?tagval:8; // bits per component
            else if (tag==259) tifc=tagval; // 1 = no compression
            else if (tag==273 && tagfmt==4) tifofs=tagval,tifofval=(taglen<=1);
            else if (tag==277) tifz=tagval; // components per pixel
          }
        }
      }
      if (tifx && tify && tifzb && (tifz==1 || tifz==3) && (tifc==1) && (tifofs && tifofs+i<n)) {
        if (!tifofval) {
          fseek(in, start+i+tifofs-7, SEEK_SET);
          for (int j=0; j<4; j++) b[j]=getc(in);
          tifofs=b[0]+(b[1]<<8)+(b[2]<<16)+(b[3]<<24);
        }
        if (tifofs && tifofs<(1<<18) && tifofs+i<n) {
          if (tifz==1 && tifzb==1) IMG_DET_NOHDR(IMAGE1,(i-7)+tifofs,((tifx-1)>>3)+1,tify);
          else if (tifz==1 && tifzb==8) IMG_DET_NOHDR(IMAGE8, (i-7)+tifofs, tifx, tify);
          else if (tifz==3 && tifzb==8) IMG_DET_NOHDR(IMAGE24, (i-7)+tifofs, tifx*3, tify);
        }
      }
      fseek(in, savedpos, SEEK_SET);
    }
  }
  return type;
}

void encode_default(FILE* in, FILE* out, int len) {
  while (len--) putc(getc(in), out);
}

int decode_default(FILE* in) {
  return getc(in);
}

void encode_jpeg(FILE* in, FILE* out, int len) {
  while (len--) putc(getc(in), out);
}

int decode_jpeg(FILE* in) {
  return getc(in);
}

void encode_bmp(FILE* in, FILE* out, int len, int width) {
  fprintf(out, "%c%c%c%c", width>>24, width>>16, width>>8, width);
  int r,g,b;
  for (int i=0; i<len/width; i++) {
    for (int j=0; j<width/3; j++) {
      b=getc(in), g=getc(in), r=getc(in);
      putc(g, out);
      putc(g-r, out);
      putc(g-b, out);
    }
    for (int j=0; j<width%3; j++) putc(getc(in), out);
  }
}

int decode_bmp(FILE *in, int &reset) {
  static int width = 0;
  if (width == 0 || reset) {
    width=getc(in)<<24;
    width|=getc(in)<<16;
    width|=getc(in)<<8;
    width|=getc(in);
    reset=0;
  }

  static int r,g,b;
  static int state1 = 0, state2 = 0;

  if (state1 < width/3) {
    if (state2 == 0) {
      b=getc(in), g=getc(in), r=getc(in);
      ++state2;
      return (b-r)&255;
    } else if (state2 == 1) {
      ++state2;
      return b;
    } else if (state2 == 2) {
      ++state1;
      if (width%3 == 0) state1 = 0;
      state2 = 0;
      return (b-g)&255;
    }
  } else {
    ++state2;
    if (state2 == width%3) {
      state1 = 0;
      state2 = 0;
    }
    return getc(in);
  }
  return -1;
}

void encode_exe(FILE* in, FILE* out, int len, int begin) {
  const int BLOCK=0x10000;
  std::vector<U8> blk(BLOCK);
  fprintf(out, "%c%c%c%c", len>>24, len>>16, len>>8, len);
  fprintf(out, "%c%c%c%c", begin>>24, begin>>16, begin>>8, begin); 

  for (int offset=0; offset<len; offset+=BLOCK) {
    int size=min(len-offset, BLOCK);
    int bytesRead=fread(&blk[0], 1, size, in);
    if (bytesRead!=size) abort();
    for (int i=bytesRead-1; i>=5; --i) {
      if ((blk[i-4]==0xe8 || blk[i-4]==0xe9 || (blk[i-5]==0x0f && (blk[i-4]&0xf0)==0x80))
         && (blk[i]==0||blk[i]==0xff)) {
        int a=(blk[i-3]|blk[i-2]<<8|blk[i-1]<<16|blk[i]<<24)+offset+begin+i+1;
        a<<=7;
        a>>=7;
        blk[i]=a>>24;
        blk[i-1]=a^176;
        blk[i-2]=(a>>8)^176;
        blk[i-3]=(a>>16)^176;
      }
    }
    fwrite(&blk[0], 1, bytesRead, out);
  }
}

int decode_exe(FILE* in) {
  const int BLOCK=0x10000;
  static int offset=0, q=0;
  static int size=0;
  static int begin=0;
  static U8 c[6];

  while (offset==size && q==0) {
    offset=0;
    size=getc(in)<<24;
    size|=getc(in)<<16;
    size|=getc(in)<<8;
    size|=getc(in);
    begin=getc(in)<<24;
    begin|=getc(in)<<16;
    begin|=getc(in)<<8;
    begin|=getc(in);
  }

  while (offset<size && q<6) {
    memmove(c+1, c, 5);
    c[0]=getc(in);
    ++q;
    ++offset;
  }

  if (q==6 && (c[0]==0x00 || c[0]==0xFF) && (c[4]==0xE8 || c[4]==0xE9 || (c[5]==0x0F && (c[4]&0xF0)==0x80))
   && (((offset-1)^(offset-6))&-BLOCK)==0 && offset<=size) {
    int a=((c[1]^176)|(c[2]^176)<<8|(c[3]^176)<<16|c[0]<<24)-offset-begin;
    a<<=7;
    a>>=7;
    c[3]=a;
    c[2]=a>>8;
    c[1]=a>>16;
    c[0]=a>>24;
  }

  return c[--q];
}

void encode_text(FILE* in, FILE* out, int len, string temp_path,
    FILE* dictionary) {
  string path1 = temp_path + "2";
  FILE* temp_input = fopen(path1.c_str(), "wb+");
  if (!temp_input) abort();

  for (int i = 0; i < len; ++i) {
    putc(getc(in), temp_input);
  }
  rewind(temp_input);

  string path2 = temp_path + "3";
  FILE* temp_output = fopen(path2.c_str(), "wb+");
  if (!temp_output) abort();

  WRT wrt;
  wrt.defaultSettings(0, NULL);
  wrt.WRT_start_encoding(temp_input, temp_output, len, false, dictionary);

  int size = ftell(temp_output);
  if (size > len - 50) {
    for (int i = 0; i < 4; ++i) {
      putc(0, out);
    }
    rewind(temp_input);
    for (int i = 0; i < len; ++i) {
      putc(getc(temp_input), out);
    }
  } else {
    rewind(temp_output);
    for (int i = 0; i < size; ++i) {
      putc(getc(temp_output), out);
    }
  }

  fclose(temp_input);
  fclose(temp_output);
  remove(path1.c_str());
  remove(path2.c_str());
}

FILE* wrt_temp;
WRT* wrt_decoder = NULL;
bool wrt_enabled = true;

void reset_text_decoder(FILE* in, string path) {
  if (wrt_temp) {
    fclose(wrt_temp);
    remove(path.c_str());
  }
  wrt_temp = fopen(path.c_str(), "wb+");
  if (!wrt_temp) abort();

  unsigned int size = 0;
  for (int i = 4; i != 0; --i) {
    int c = getc(in);
    size = size * 256 + c;
    putc(c, wrt_temp);
  }
  if (size == 0) {
    wrt_enabled = false;
    return;
  }
  wrt_enabled = true;
  size -= 8;

  for (unsigned int i = 0; i < size; ++i) {
    putc(getc(in), wrt_temp);
  }
  rewind(wrt_temp);

  if (wrt_decoder != NULL) delete wrt_decoder;
  wrt_decoder = new WRT();
  wrt_decoder->defaultSettings(0, NULL);
  wrt_decoder->WRT_prepare_decoding();
}

int decode_text(FILE* in, FILE* dictionary) {
  if (!wrt_enabled) return getc(in);
  return wrt_decoder->WRT_decode_char(wrt_temp, NULL, 0, dictionary);
}

void Encode(FILE* in, FILE* out, int n, string temp_path, FILE* dictionary) {
  Filetype type=DEFAULT;
  long begin=ftell(in);

  long start = begin;
  int remainder = n;
  int text_bytes = 0;
  while (remainder > 0) {
    Filetype nextType=detect(in, remainder, type);
    long end=ftell(in);
    int len=int(end-begin);
    if (type == TEXT) text_bytes += len;
    remainder-=len;
    type=nextType;
    begin=end;
  }
  fseek(in, start, SEEK_SET);
  type = DEFAULT;
  begin = start;

  double text_fraction = text_bytes;
  text_fraction /= n;
  if (text_fraction > 0.95) {
    fprintf(out, "%c%c%c%c%c", TEXT, n>>24, n>>16, n>>8, n);
    encode_text(in, out, n, temp_path, dictionary);
    return;
  }

  while (n>0) {
    Filetype nextType=detect(in, n, type);
    long end=ftell(in);
    fseek(in, begin, SEEK_SET);
    int len=int(end-begin);
    if (len>0) {
      fprintf(out, "%c%c%c%c%c", type, len>>24, len>>16, len>>8, len);
      switch(type) {
        case JPEG:    encode_jpeg(in, out, len); break;
        case IMAGE24: encode_bmp(in, out, len, info); break;
        case EXE:     encode_exe(in, out, len, begin); break;
        case TEXT:    encode_text(in, out, len, temp_path, dictionary); break;
        default: {
          if (HasInfo(type))
            fprintf(out, "%c%c%c%c", info>>24, info>>16, info>>8, info); // write info
          encode_default(in, out, len); break;
        }
      }
    }
    n-=len;
    type=nextType;
    begin=end;
  }
}

void NoPreprocess(FILE* in, FILE* out, int n) {
  fprintf(out, "%c%c%c%c%c", DEFAULT, n>>24, n>>16, n>>8, n);
  encode_default(in, out, n);
}

int decode2(FILE* in, string temp_path, FILE* dictionary) {
  static Filetype type=DEFAULT;
  static int len=0, reset=0;
  while (len==0) {
    int c = getc(in);
    if (c == EOF) return -1;
    reset=1;
    type=(Filetype)c;
    len=getc(in)<<24;
    len|=getc(in)<<16;
    len|=getc(in)<<8;
    len|=getc(in);
    if (len<0) len=1;
    if (type == TEXT) reset_text_decoder(in, temp_path);
  }
  --len;
  switch (type) {
    case JPEG:    return decode_jpeg(in);
    case IMAGE24: return decode_bmp(in, reset);
    case EXE:     return decode_exe(in);
    case TEXT:    return decode_text(in, dictionary);
    default: {
      if (reset && HasInfo(type)){
        for (int i=info=0;i<4;i++) info=(info<<8)|getc(in); //read info
        reset=0;
      }
      return decode_default(in);
    }
  }
}

void Decode(FILE* in, FILE* out, string temp_path, FILE* dictionary) {
  string path = temp_path + "2";
  while (true) {
    int result = decode2(in, path, dictionary);
    if (result == -1) {
      if (wrt_temp) {
        fclose(wrt_temp);
        remove(path.c_str());
      }
      return;
    }
    putc(result, out);
  }
}

}
