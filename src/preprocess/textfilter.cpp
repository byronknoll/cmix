// TextFilter 3.0 for PAQ (based on WRT 4.6) by P.Skibinski, 02.03.2006, inikep@o2.pl

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <map>
#include <string>
#ifdef WIN32
	#include <io.h>
	#include <windows.h>
#else
	#include <sys/types.h>
	#include <dirent.h>
#endif

#define USE_EOLC 0

#define CHAR_FIRSTUPPER		64 	// for encode lower word with first capital letter
#define CHAR_UPPERWORD		7	// for encode upper word
#define CHAR_LOWERWORD		6	// for encode lower word with a few capital letter
#define CHAR_PUNCTUATION	8	// for punctuation marks modeling
#define CHAR_NOSPACE		8   // the same as CHAR_PUNCTUATION
#define CHAR_CR_LF			14
#define CHAR_ESCAPE			12	// for encode reserved chars (CHAR_ESCAPE,CHAR_FIRSTUPPER,...)
#define NGRAM_FIRST			'A'
#define NGRAM_LAST			'Z'
#define BINARY_FIRST		128
#define BINARY_LAST			255

#define AUTO_SWITCH			8	// param for !OPTION_NORMAL_TEXT_FILTER
#define WORD_MIN_SIZE		1
#define FUNCTION_CHECK_ERRORS
#define WRT_HEADER "WRT4"

//#define min(a,b) (((a)>(b))?(b):(a))
#ifndef SHORTEN_CODE
#define TOLOWER(c) ((c>='A' && c<='Z')?(c+32):((upperSet[0][c]>0)?lowerSetRev[0][upperSet[0][c]]:c))
#define TOUPPER(c) ((c>='a' && c<='z')?(c-32):((lowerSet[0][c]>0)?upperSetRev[0][lowerSet[0][c]]:c))
#define TOUPPER_SET(c) ((c>='a' && c<='z')?(c-32):((lowerSet[usedSet][c]>0)?upperSetRev[usedSet][lowerSet[usedSet][c]]:c))
#else
#define TOLOWER(c)	((c>='A' && c<='Z')?(c+32):c)
#define TOUPPER(c)	((c>='a' && c<='z')?(c-32):c)
#define TOUPPER_SET	TOUPPER
#endif

#define OPTION(option) ((wrt.preprocFlag & option)!=0)
#define TURN_OFF(option) ;//{if (preprocFlag & option) preprocFlag-=option;}
#define TURN_ON(option)	 ;//{preprocFlag|=option;}
#define RESET_OPTIONS 	 ;//preprocFlag=0

#define COND_BIN_FILTER(c) (((c<32)?(c!=10 && c!=13):(0)) || (c>=BINARY_FIRST)) 

#define PRINT_CHARS(data) ;//printf data
#define PRINT_CODEWORDS(data) ;//printf data
#define PRINT_DICT(data) ;//printf data
//#define LOG_ARITHMETIC_ENCODER

#define HASH_TABLE_SIZE		(1<<21)

using namespace std;

int word_hash[HASH_TABLE_SIZE];

typedef unsigned int  uint;
typedef unsigned char uc;

// filesize() function

static uint flen( FILE* f )
{
	fseek( f, 0, SEEK_END );
	uint len = ftell(f);
	fseek( f, 0, SEEK_SET );
	return len;
}

class WRT
{
public:

bool fileCorrupted=false;

WRT() : restartEnc(false), WRT_verbose(false), preprocType(PAQ), dict(NULL), dictlen(NULL), dictmem(NULL), langCount(0), lastShortDict(-1) { };
~WRT() { WRT_deinitialize(); freeNames(); }

enum EPreprocessType { LZ77, BWT, PPM, PAQ };
enum EWordType { LOWERWORD, FIRSTUPPER, UPPERWORD };
enum EEOLType { UNDEFINED, CRLF, LF };
enum EUpperType { UFALSE, UTRUE, FORCE };
enum ESpaceType { NONE, SPACE, EOL };


bool restartEnc,initOrder,forceNormalTextFilter,forceWordSurroroundModeling,forceEOLcoding;
int tryShorterBound,preprocessing,s_size,WRTd_c,WRTd_qstart,WRTd_qend,WRTd_type;
int fftell,fftelld,originalFileLen,autoSwitch,WRTd_binCount;
int bufferedChar,lastEOL,EOLcount,lastChar,llast,llbckp;
bool swapCase,WRT_verbose,WRTd_upper;
unsigned char WRTd_s[1024];
unsigned char WRTd_queue[128];
EUpperType upperWord;
EEOLType EOLType;
ESpaceType spaceBefore;
EPreprocessType preprocType;

#ifdef POWERED_BY_PAQ
#define DECODE_GETC(c,file)\
{\
	if (fftelld<originalFileLen) \
	{ \
		c=WRTd_filter->read(); \
		fftelld++; \
	} \
	else \
		c=EOF; \
}
#else
#define DECODE_GETC(c,file)\
{\
	if (fftelld<originalFileLen) \
	{ \
		c=getc(file); \
		fftelld++; \
	} \
	else \
		c=EOF; \
}
#endif

#define ENCODE_PUTC(c,file)\
{ \
	putc(c,file); \
}

#define MAX_FREQ_ORDER1		2520
#define ORDER1_STEP	4
	
int mZero[MAX_FREQ_ORDER1];
int mOne[MAX_FREQ_ORDER1];

#define	UPDATE_ORDER1(prev,value)	UpdateOrder1(prev,value,ORDER1_STEP)
#define	ENCODE_ORDER1(prev,value)	EncodeOrder1(prev,value)
#define	DECODE_ORDER1(prev)		DecodeOrder1(prev)
#define	INIT_ORDER1			InitOrder1(MAX_FREQ_ORDER1)

#define DICTNAME_EXT "dic"
#define DICTNAME "english"
#define SHORT_DICTNAME "english_short"
#define WRT_DICT_DIR "./"

#define HASH_DOUBLE_MULT	29
#define HASH_MULT		23

int sizeDict;
unsigned char** dict;
unsigned char* dictlen;
unsigned char* dictmem;

int ngram_hash[256][256];

#define CHARSET_COUNT		6

int lowerSet[CHARSET_COUNT][256];
int upperSet[CHARSET_COUNT][256];
int lowerSetRev[CHARSET_COUNT][256];
int upperSetRev[CHARSET_COUNT][256];
int freeUpper[CHARSET_COUNT],freeLower[CHARSET_COUNT];
int usedSet;

int reservedSet[256]; 
int addSymbols[256]; 
int sym2codeword[256]; 
int codeword2sym[256]; 
int value[256];

int dictionary,dict1size,dict2size,dict3size,dict4size,dict1plus2plus3,dict1plus2;
int bound4,bound3,dict123size,dict12size;



// convert upper string to lower
inline void toLower(unsigned char* s,int s_size)
{
	for (int i=0; i<s_size; i++)
		s[i]=TOLOWER(s[i]);
}


// convert lower string to upper
inline void toUpper(unsigned char* s,int s_size)
{
	for (int i=0; i<s_size; i++)
		s[i]=TOUPPER(s[i]); 
}

#ifndef SHORTEN_CODE
#define ORIGINAL_CHARSET(c)\
{\
	if (usedSet>0)\
	{\
		if (lowerSet[0][c]>0)\
			c=lowerSetRev[usedSet][lowerSet[0][c]];\
		else\
		if (upperSet[0][c]>0)\
			c=upperSetRev[usedSet][upperSet[0][c]];\
	}\
}
#else
#define ORIGINAL_CHARSET(c) c
#endif

// make hash from string
inline unsigned int stringHash(const unsigned char *ptr, int len)
{
	unsigned int hash;
	for (hash = 0; len>0; len--, ptr++)
		hash = HASH_MULT * hash + *ptr;
 
	return hash&(HASH_TABLE_SIZE-1);
}


// check if word "s" does exist in the dictionary using hash "h" 
inline int checkHashExactly(const unsigned char* s,int s_size,int h)
{
	int i;

	i=word_hash[h];
	if (i>0)
	{
		if (dictlen[i]!=s_size || memcmp(dict[i],s,s_size)!=0)
		{
			i=word_hash[(h+s_size*HASH_DOUBLE_MULT)&(HASH_TABLE_SIZE-1)];
			if (i>0)
			{
				if (dictlen[i]!=s_size || memcmp(dict[i],s,s_size)!=0)
				{
					i=word_hash[(h+s_size*HASH_DOUBLE_MULT*HASH_DOUBLE_MULT)&(HASH_TABLE_SIZE-1)];
					if (i>0)
					{
						if (dictlen[i]!=s_size || memcmp(dict[i],s,s_size)!=0)
							i=-1;
					}
					else
						i=-1;
				}
			}
			else
				i=-1;
		}
	}
	else
		i=-1;

	if (i>dictionary)
		i=-1;

	return i;
}

// check if word "s" (prefix of original word) does exist in the dictionary using hash "h" 
inline int checkHash(const unsigned char* s,int s_size,int h)
{
	int i;

	i=word_hash[h];
	if (i>0)
	{
		if (dictlen[i]>s_size || memcmp(dict[i],s,s_size)!=0)
		{
			i=word_hash[(h+s_size*HASH_DOUBLE_MULT)&(HASH_TABLE_SIZE-1)];
			if (i>0)
			{
				if (dictlen[i]>s_size || memcmp(dict[i],s,s_size)!=0)
				{
					i=word_hash[(h+s_size*HASH_DOUBLE_MULT*HASH_DOUBLE_MULT)&(HASH_TABLE_SIZE-1)];
					if (i>0)
					{
						if (dictlen[i]>s_size || memcmp(dict[i],s,s_size)!=0)
							i=-1;
					}
					else
						i=-1;
				}
			}
			else
				i=-1;
		}
	}
	else
		i=-1;

	if (i>dictionary)
		i=-1;

	return i;
}


// check if word "s" or prefix of word "s" does exist in the dictionary using hash "h" 
inline int findShorterWord(const unsigned char* s,int s_size)
{
	int ret, i, best;
	unsigned int hash;

	hash = 0;
	for (i=0; i<WORD_MIN_SIZE+tryShorterBound; i++)
		hash = HASH_MULT * hash + s[i];
 
	best=-1;
	for (; i<s_size; i++)
	{
		ret=checkHash(s,i,hash&(HASH_TABLE_SIZE-1));	
		if (ret>=0)
			best=ret;
		hash = HASH_MULT * hash + s[i];
	}

	return best;
}

inline int findShorterWordRev(const unsigned char* s,int s_size)
{
	int ret, i;

	for (i=s_size-1; i>=WORD_MIN_SIZE+tryShorterBound; i--)
	{
		ret=checkHash(s+s_size-i,i,stringHash(s+s_size-i,i));	
		if (ret>=0)
			return ret;
	}

	return -1;
}

// encode word (should be lower case) using n-gram array (when word doesn't exist in the dictionary)
inline void encodeAsText(unsigned char* s,int s_size,FILE* fileout)
{
	int i,ngram;

	if (spaceBefore!=NONE)
	{
		if (spaceBefore==SPACE)
			ENCODE_PUTC(' ',fileout);
		spaceBefore=NONE;
	}

	if (usedSet>0)
	{
		for (i=0; i<s_size; i++)
		{
			ORIGINAL_CHARSET(s[i]);
			ENCODE_PUTC(s[i],fileout);
		}
	}
	else
	{
		ngram=0;
		for (i=0; i<s_size; )
		{
			if (IF_OPTION(OPTION_USE_NGRAMS))
				ngram=ngram_hash[s[i]][s[i+1]];

			if (ngram>0 && ngram<dict1size)	///// && preprocType!=LZ77)
			{
				encodeCodeWord(ngram,fileout);
				i+=2;
			}
			else
			{
				ENCODE_PUTC(s[i],fileout);
				i++;
			}
		}
	}
}

inline void encodeCodeWord(int i,FILE* fileout)
{
	int first,second,third;

	first=i-1;

	if (first>=80*49) //bound3)
	{
		first-=80*49; //bound3;

		third=first/dict12size;		
		first=first%dict12size;
		second=first/dict1size;		
		first=first%dict1size;

		ENCODE_PUTC(sym2codeword[dict1plus2+third],fileout);
		PRINT_CODEWORDS(("1st=%d(%d) ",sym2codeword[dict1plus2+third],third));

		ENCODE_PUTC(sym2codeword[dict1size+second],fileout);
		PRINT_CODEWORDS(("2nd=%d(%d) ",sym2codeword[dict1size+second],second));

		ENCODE_PUTC(sym2codeword[first],fileout);
		PRINT_CODEWORDS(("3rd=%d(%d) ",sym2codeword[first],first));
	}
	else
		if (first>=dict1size)
		{
			first-=dict1size;

			second=first/dict1size;		
			first=first%dict1size;

			ENCODE_PUTC(sym2codeword[dict1size+second],fileout);
			PRINT_CODEWORDS(("1st=%d ",sym2codeword[dict1size+second]));
	
			ENCODE_PUTC(sym2codeword[first],fileout);
			PRINT_CODEWORDS(("2nd=%d ",sym2codeword[first]));
		}
		else
		{
			ENCODE_PUTC(sym2codeword[first],fileout);
			PRINT_CODEWORDS(("1st=%d ",sym2codeword[first]));
		}

		PRINT_CODEWORDS((" no=%d %s\n", no-1,dict[no]));
}

// encode word "s" using dictionary
inline void encodeWord(FILE* fileout,unsigned char* s,int s_size,EWordType wordType)
{
	int i,j,d,e;
	int size=0;
	int flagToEncode=-1;

	if (s_size<1)
	{
		if (spaceBefore!=NONE)
		{
			if (spaceBefore==SPACE)
				ENCODE_PUTC(' ',fileout);
			spaceBefore=NONE;
		}
		return;
	}

	s[s_size]=0;

	if (wordType!=LOWERWORD)
	{
		if (IF_OPTION(OPTION_CAPITAL_CONVERSION))
		{
		if (wordType==FIRSTUPPER)
		{
			flagToEncode=CHAR_FIRSTUPPER;
			s[0]=TOLOWER(s[0]);
		}
		else // wordType==UPPERWORD
		{
			flagToEncode=CHAR_UPPERWORD;
			toLower(s,s_size);
		}
		}
		else
			wordType=LOWERWORD;
	}
	

	if (IF_OPTION(OPTION_USE_DICTIONARY) && s_size>=WORD_MIN_SIZE)
	{
		i=checkHashExactly(s,s_size,stringHash(s,s_size));
		PRINT_CODEWORDS(("checkHashExactly i=%d %d=%s\n",i,s_size,s));

		if (i<0 && IF_OPTION(OPTION_TRY_SHORTER_WORD))
		{
			// try to find shorter version of word in dictionary
			i=findShorterWord(s,s_size);
			j=findShorterWordRev(s,s_size);
			PRINT_CODEWORDS(("findShorterWord i=%d\n",i));

			d=e=0;
			if (i>=0) d=dictlen[i]-(i>80)-(i>3920)-1;
			if (j>=0) e=dictlen[j]-(j>80)-(j>3920)-1;
			if (d>=e) { if (d> 0) size=dictlen[i]; }
			else
				if (!IF_OPTION(OPTION_SPACELESS_WORDS))
				{
					i=j;
					PRINT_CODEWORDS(("findShorterWordRev i=%d\n",i));
					if (e> 0)
					{
						if (wordType!=LOWERWORD)
						{
							ENCODE_PUTC(flagToEncode,fileout);
							if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
								ENCODE_PUTC(' ',fileout);
							wordType=LOWERWORD;
						}
						
						s[s_size-dictlen[i]]=0;
						encodeAsText(s,s_size-dictlen[i],fileout);
                                                  ENCODE_PUTC(CHAR_NOSPACE,fileout);
						s+=dictlen[i];
						s_size-=dictlen[i];
					}
				}
		}
	}
	else
		i=-1;


	if (i>=0)
	{
		if (wordType!=LOWERWORD)
		{
			ENCODE_PUTC(flagToEncode,fileout);
			if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
				ENCODE_PUTC(' ',fileout);
		}

		if (spaceBefore==NONE)
		{
			if (IF_OPTION(OPTION_SPACELESS_WORDS))
				ENCODE_PUTC(CHAR_NOSPACE,fileout);
		}
		else
			spaceBefore=NONE;

		//////if (preprocType==LZ77)
		//////	encodeCodeWord_LZ(i,fileout);
		//////else
			encodeCodeWord(i,fileout);

		if (size>0)
			encodeAsText(s+size,s_size-size,fileout);
	}
	else
	{
		if (wordType!=LOWERWORD)
		{
			if (spaceBefore!=NONE)
			{
				if (spaceBefore==SPACE)
					ENCODE_PUTC(' ',fileout);
				spaceBefore=NONE;
			}

			ENCODE_PUTC(flagToEncode,fileout);
			if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
				ENCODE_PUTC(' ',fileout);
		}
 
		encodeAsText(s,s_size,fileout);
	}
}


// decode word using dictionary
#define DECODE_WORD(dictNo,i)\
{\
		switch (dictNo)\
		{\
			case 4:\
				i+=bound4;\
				break;\
			case 3:\
				i+=80*49; /*bound3;*/\
				break;\
			case 2:\
				i+=dict1size;\
				break;\
		}\
\
		if (i>=0 && i<sizeDict)\
		{\
			PRINT_CODEWORDS(("i=%d ",i)); \
			i++;\
			s_size=dictlen[i];\
			memcpy(s,dict[i],s_size+1);\
			PRINT_CODEWORDS(("%s\n",dict[i])); \
		}\
		else\
		{\
			printf("File is corrupted!\n");\
			fileCorrupted=true;\
		}\
}



inline int decodeCodeWord(FILE* file, unsigned char* s,int& c)
{
	int i,dictNo,s_size=0;

	if (codeword2sym[c]<dict1size)
	{
		i=codeword2sym[c];
		dictNo=1;
		DECODE_WORD(dictNo, i);
		return s_size;
	}
		i=dict1size*(codeword2sym[c]-dict1size);

	DECODE_GETC(c,file);
	PRINT_CODEWORDS(("DC1 c=%d i=%d\n",c,i));

	if (codeword2sym[c]<dict1size)
	{
		i+=codeword2sym[c];
		dictNo=2;
		DECODE_WORD(dictNo, i);
		return s_size;
	}
	{
		i=(i-dict12size)*dict2size;
		PRINT_CODEWORDS(("DC2b c=%d\n",codeword2sym[c]-dict1size));
		i+=dict1size*(codeword2sym[c]-dict1size);
	}

	DECODE_GETC(c,file);
	PRINT_CODEWORDS(("DC2 c=%d i=%d\n",c,i));

	{
		PRINT_CODEWORDS(("DC3b c=%d\n",codeword2sym[c]));
		i+=codeword2sym[c];
		dictNo=3;
		DECODE_WORD(dictNo, i);
		return s_size;
	}
}

unsigned char* loadDictionary(FILE* file,unsigned char* mem,int word_count)
{
	unsigned char* word;
	int i,j,c,collision,bound;

	collision=0;
	bound=sizeDict+word_count;


	while (!feof(file))
	{
		word=mem;
		do
		{
			c=getc(file);
#ifndef SHORTEN_CODE
			if (usedSet==CHARSET_COUNT-1)
			{
				if (lowerSet[0][c]>0)
					c=lowerSetRev[CHARSET_COUNT-1][lowerSet[0][c]];
			}
#endif
			word[0]=c;
			word++;
		}
		while (c>32);

		if (c==EOF)
			break;
		if (c=='\r') 
			c=getc(file);
		
		word[-1]=0;
		i=word-mem-1;
		
		dictlen[sizeDict]=i;
		dict[sizeDict]=mem;

		j=stringHash(mem,i);
		mem+=(i/4+1)*4;
		

		if (word_hash[j]!=0)
		{
			if (dictlen[sizeDict]!=dictlen[word_hash[j]] || memcmp(dict[sizeDict],dict[word_hash[j]],dictlen[sizeDict])!=0)
			{
				c=(j+i*HASH_DOUBLE_MULT)&(HASH_TABLE_SIZE-1);
				if (word_hash[c]!=0)
				{
					if (dictlen[sizeDict]!=dictlen[word_hash[c]] || memcmp(dict[sizeDict],dict[word_hash[c]],dictlen[sizeDict])!=0)
					{
						c=(j+i*HASH_DOUBLE_MULT*HASH_DOUBLE_MULT)&(HASH_TABLE_SIZE-1);
						if (word_hash[c]!=0)
						{
							collision++;
						}
						else
						{
								word_hash[c]=sizeDict++;
						}
					}
				}
				else
				{
							word_hash[c]=sizeDict++;
				}
			}
		}
		else
		{
				word_hash[j]=sizeDict++;
		}
		
		if (sizeDict>dictionary || sizeDict>=bound)
		{
			sizeDict--;
			break;
		}
	}

	return mem;
}

int loadCharset(FILE* file,int& freeChar,int* charset,int* charsetRev,bool *joinCharsets=NULL)
{
	int c,res,mult;
	
	res=0;
	c=getc(file); 
	mult=100;
#ifndef SHORTEN_CODE
	while (c>32)
	{
		{
			if (joinCharsets)
				joinCharsets[c]=true;

			charsetRev[freeChar]=c;
			charset[c]=freeChar++;

		}

		res+=mult*value[c];
		mult--;

		c=getc(file);
	}

#endif
	if (c==13)
		c=getc(file); // skip CR+LF or LF
	return res;
}

void initializeCodeWords()
{
	int c,charsUsed;

	for (c=0; c<256; c++)
	{
		addSymbols[c]=0;
		codeword2sym[c]=0;
		sym2codeword[c]=0;
		reservedSet[c]=0;
	}

	for (c=BINARY_FIRST; c<=BINARY_LAST; c++)
		addSymbols[c]=1;


	if (IF_OPTION(OPTION_ADD_SYMBOLS_MISC))
	{
		addSymbols[35]=1;
		addSymbols[38]=1;
		addSymbols[60]=1;
		addSymbols[62]=1;
		addSymbols[64]=1;
		addSymbols[94]=1;
		addSymbols[96]=1;
		addSymbols[123]=1;
		addSymbols[124]=1;
		addSymbols[125]=1;
		addSymbols[126]=1;
	}

	if (IF_OPTION(OPTION_ADD_SYMBOLS_0_5))
		for (c=0; c<=5; c++)
			addSymbols[c]=1;

	if (IF_OPTION(OPTION_ADD_SYMBOLS_14_31))
		for (c=14; c<=31; c++)
			addSymbols[c]=1;

	for (c=0; c<256; c++)
	{
		if (IF_OPTION(OPTION_USE_DICTIONARY) && addSymbols[c]
			&& (lowerSet[usedSet][c]!=0))
			addSymbols[c]=0;

		if (IF_OPTION(OPTION_USE_DICTIONARY) && addSymbols[c]
			&& (upperSet[usedSet][c]!=0) && !IF_OPTION(OPTION_CAPITAL_CONVERSION))
			addSymbols[c]=0;

		if ( (IF_OPTION(OPTION_USE_DICTIONARY) && addSymbols[c])
			|| c==CHAR_ESCAPE || c==CHAR_LOWERWORD || c==CHAR_FIRSTUPPER || c==CHAR_UPPERWORD || c==CHAR_CR_LF
			|| c==CHAR_NOSPACE
			|| (IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING) && c==CHAR_PUNCTUATION))
			reservedSet[c]=1;
		else
			reservedSet[c]=0;
	}

	if (IF_OPTION(OPTION_ADD_SYMBOLS_A_Z) && IF_OPTION(OPTION_CAPITAL_CONVERSION))
		for (c=NGRAM_FIRST; c<=NGRAM_LAST; c++)
			addSymbols[c]=1;


	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
		charsUsed=0;
		for (c=0; c<256; c++)
		{
			if (addSymbols[c])
			{
				codeword2sym[c]=charsUsed;
				sym2codeword[charsUsed]=c;
				charsUsed++;

			}
		}

					dict1size=80;
					dict2size=32;
					dict3size=16;
					dict4size=0;

		bound4=dict1size*dict2size*dict3size+dict1size*dict2size+dict1size;
		bound3=dict1size*dict2size+dict1size;
		dict123size=dict1size*dict2size*dict3size;
		dict12size=dict1size*dict2size;
		dict1plus2=dict1size+dict2size;
		dict1plus2plus3=dict1size+dict2size+dict3size;
		dictionary=dict123size + dict1size*(dict2size+dict3size+1);
		
		dict=(unsigned char**)calloc(sizeof(unsigned char*)*(dictionary+1),1);
		dictlen=(unsigned char*)calloc(sizeof(unsigned char)*(dictionary+1),1);

		PRINT_DICT(("usedSet=%d preprocType=%d %d %d %d %d(%d) charsUsed=%d sizeDict=%d\n",usedSet,preprocType,dict1size,dict2size,dict3size,dict4size,dictionary,charsUsed,sizeDict));
	}
}


// read dictionary from files to arrays
bool initialize(unsigned char* dictName,unsigned char* shortDictName,bool encoding, FILE* english_dictionary)
{
	PRINT_DICT(("dictName=%s shortDictName=%s\n",dictName,shortDictName));

	int i,j,c,set[CHARSET_COUNT],fileLen;
	FILE* file,*file2;
	unsigned char* mem;

	WRT_deinitialize();
	sizeDict=0;

	memset(&word_hash[0],0,HASH_TABLE_SIZE*sizeof(word_hash[0]));
	memset(lowerSet,0,sizeof(lowerSet));
	memset(upperSet,0,sizeof(upperSet));
	memset(lowerSetRev,0,sizeof(lowerSetRev));
	memset(upperSetRev,0,sizeof(upperSetRev));
	memset(ngram_hash,0,sizeof(ngram_hash));

	if (dictName==NULL && shortDictName==NULL)
	{
		initializeCodeWords();
		return true;
	}

	if (dictName==NULL && shortDictName!=NULL)
	{
		dictName=shortDictName;
		shortDictName=NULL;
	}


	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
    file = english_dictionary;
		if (file==NULL)
		{
			printf("Can't open dictionary %s\n",dictName);
			return false;
		}

		fileLen=flen(file);
		if (fscanf(file,"%d",&dict123size) < 0) abort();

		do { c=getc(file); } while (c>=32); if (c==13) c=getc(file); // skip CR+LF or LF

		for (i=0; i<CHARSET_COUNT; i++)
		{
			freeUpper[i]=1;
			set[i]=loadCharset(file,freeUpper[i],upperSet[i],upperSetRev[i]);
			freeLower[i]=1;
			set[i]+=loadCharset(file,freeLower[i],lowerSet[i],lowerSetRev[i]);
		}

		if (encoding)
		{
			j=0;
			for (i=1; i<CHARSET_COUNT-1; i++)
			{
				if (set[i]>set[j])
					j=i;
			}

			usedSet=j;

			if (set[usedSet]==0 && freeLower[CHARSET_COUNT-1]>1)
				usedSet=CHARSET_COUNT-1;
		}
	
		if (freeUpper[0]!=freeLower[0])
		{
			fclose(file);
			return false;
		}

		if (freeUpper[usedSet]!=freeLower[usedSet] || freeLower[usedSet]!=freeLower[0])
		{
			fclose(file);
			return false;
		}

		dictmem=(unsigned char*)calloc(fileLen*2,1);
		mem=dictmem;
		sizeDict=1;

		if (!dictmem)
		{
			initializeCodeWords();
			return true;
		}

		if (shortDictName)
		{
			file2=fopen((const char*)shortDictName,"rb");
			if (file2==NULL)
			{
				printf("Can't open dictionary %s\n",shortDictName);
				return false;
			}
				
			do c=getc(file2); while (c>=32); if (c==13) c=getc(file2); 
			for (i=0; i<CHARSET_COUNT; i++)
			{
				loadCharset(file2,freeUpper[i],upperSet[i],upperSetRev[i]);
				loadCharset(file2,freeLower[i],lowerSet[i],lowerSetRev[i]);
			}

			if (freeUpper[0]!=freeLower[0])
			{
				fclose(file2);
				fclose(file);
				return false;
			}

			if (freeUpper[usedSet]!=freeLower[usedSet] || freeLower[usedSet]!=freeLower[0])
			{
				fclose(file2);
				fclose(file);
				return false;
			}

			initializeCodeWords();

			if (dict==NULL || dictlen==NULL)
				return false;

			mem=loadDictionary(file2,mem,dictionary);
			fclose(file2);

		}
		else
		{
			initializeCodeWords();

			if (dict==NULL || dictlen==NULL)
				return false;
		}

		mem=loadDictionary(file,mem,dictionary);

		if (encoding && usedSet==CHARSET_COUNT-1)
		{
			memset(lowerSet,0,sizeof(lowerSet));
			memset(upperSet,0,sizeof(upperSet));
			memset(lowerSetRev,0,sizeof(lowerSetRev));
			memset(upperSetRev,0,sizeof(upperSetRev));
		}
	}
	else
	{
		initializeCodeWords();
	}

	return true;
}


void WRT_deinitialize()
{
	if (dict)
	{
		free(dict);
		dict=NULL;
	}
	if (dictlen)
	{
		free(dictlen);
		dictlen=NULL;
	}
	if (dictmem)
	{
		free(dictmem);
		dictmem=NULL;
	}

	sizeDict=0;
}

#define MAX_RECORD_LEN		1024
#define MAX_DICT_NUMBER		255
#define SAMPLE_WORDS_COUNT	250
#define SAMPLE_WORDS_COUNT_MAX	(SAMPLE_WORDS_COUNT*CHARSET_COUNT)

int lang[MAX_DICT_NUMBER];
int langCount,langSum;
unsigned char* langName[MAX_DICT_NUMBER];
int longDictLen,shortDictLen,longDict,shortDict,lastShortDict;
bool joinCharsets[256];

std::multimap<std::string,int> map;
std::multimap<std::string,int>::iterator it;

inline void checkWord(unsigned char* s,int s_size)
{
	if (s_size<WORD_MIN_SIZE)
		return;

	std::string str;
	str.append((char*)s,s_size);

	it=map.find(str);
	if (it==map.end())
		return;

	do
	{
		lang[it->second/SAMPLE_WORDS_COUNT_MAX]++;
		langSum++;

		it++;
	}
	while (it!=map.end() && it->first==str);

	return;
}


#define SET_PPF 1220

int WRT_detectFileType(FILE* file, int part_length, int parts, int& recordLen)
{
	dictionary=1<<30;
	memset(&lang[0],0,MAX_DICT_NUMBER*sizeof(lang[0]));
	memset(value,0,sizeof(value));
	longDictLen=1;
	longDict=shortDictLen=0;
	shortDict=lastShortDict=-1;
	return preprocFlag;
}

void WRT_set_options(char c,char c2) {}

void WRT_get_options(int& c,int& c2)
{
	if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
		TURN_OFF(OPTION_CAPITAL_CONVERSION);

	if (!IF_OPTION(OPTION_USE_DICTIONARY))
		TURN_OFF(OPTION_USE_NGRAMS);

	if ((IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING)) && (IF_OPTION(OPTION_SPACELESS_WORDS)))
	{
		TURN_OFF(OPTION_WORD_SURROROUNDING_MODELING);
	}

	c=c2=0;
	if (IF_OPTION(OPTION_USE_NGRAMS))
		c=c+128;
	if (IF_OPTION(OPTION_USE_DICTIONARY))
		c=c+64;
#if USE_EOLC
	if (IF_OPTION(OPTION_EOL_CODING))
		c=c+32;
#endif
	if (IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING))
		c=c+16;
	if (IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
		c=c+8;
	if (IF_OPTION(OPTION_SPACE_AFTER_EOL))
		c=c+4;

	c+=preprocType;

	if (IF_OPTION(OPTION_CAPITAL_CONVERSION))
		c2=c2+32;
	if (IF_OPTION(OPTION_UTF8))
		c2=c2+16;
}

int defaultSettings(int argc, char* argv[])
{
	forceNormalTextFilter=false;
	forceWordSurroroundModeling=false;

	RESET_OPTIONS;

#if USE_EOLC
	forceEOLcoding=false;
	TURN_ON(OPTION_EOL_CODING);
#endif
	TURN_ON(OPTION_TRY_SHORTER_WORD);
	TURN_ON(OPTION_USE_DICTIONARY);
	TURN_ON(OPTION_NORMAL_TEXT_FILTER);
	TURN_ON(OPTION_CAPITAL_CONVERSION);
	TURN_ON(OPTION_UTF8);
                        tryShorterBound=6;
			TURN_ON(OPTION_WORD_SURROROUNDING_MODELING);
	int optCount=0;
	return optCount;
}


inline bool addWord(std::string s,int& sizeFullDict)
{
	std::pair<std::string,int> pair(s,sizeFullDict++);

	map.insert(pair);

	return true;
}

bool readDicts(const char* pattern,char* dictPath,int dictPathLen, FILE* english_dictionary)
{
	FILE* file;
	bool nonlatin;
	int c,i,sizeDict,sizeFullDict=0;

	memset(joinCharsets,0,sizeof(joinCharsets));

	map.clear();

	langSum=0;
	langCount=0;
	sizeDict=0;
	dictPath[dictPathLen]=0;
	strcat(dictPath,pattern);
#ifdef WIN32
    struct _finddata_t c_file;
    long hFile=_findfirst( dictPath, &c_file );

	if (hFile != -1L) 
	do
	{
		char* filename=c_file.name;
#else
		{
		const char* filename="english.dic";
#endif

		dictPath[dictPathLen]=0;
		strcat(dictPath,filename);

    file = english_dictionary;
    rewind(file);
		if (file==NULL)
			return false; //continue;

		i=strlen(filename);

		langName[langCount]=(unsigned char*)malloc(i+1);
		memcpy(langName[langCount],(const char*)filename,i+1);

		memset(lowerSet,0,sizeof(lowerSet));
		memset(lowerSetRev,0,sizeof(lowerSetRev));

		do c=getc(file); while (c>=32); if (c==13) c=getc(file); 

		for (i=0; i<CHARSET_COUNT; i++)
		{
			do  c=getc(file); while (c>=32); if (c==13) c=getc(file); 
			freeLower[i]=1;
			loadCharset(file,freeLower[i],lowerSet[i],lowerSetRev[i],joinCharsets);
		}

		sizeFullDict=sizeDict;
		std::string s;

		while (!feof(file))
		{
			s.erase();
			nonlatin=false;
			while (true)
			{
				c=getc(file);
				if (c<32)
					break;
				s.append(1,(char)c);
				if (lowerSet[0][c]>0)
					nonlatin=true;
			}
			if (c==EOF)
				break;

			if (c==13)
				c=getc(file); // skip CR+LF or LF

			if (addWord(s,sizeFullDict))
			{
				sizeDict++;
				if (sizeDict%SAMPLE_WORDS_COUNT==0)
					break;
			}
			else
				continue;


			if (nonlatin)
			{
				std::string t;

				for (i=1; i<CHARSET_COUNT-1; i++)
				{
					if (freeLower[i]>1)
					{
						t.erase();
	
						for (c=0; c<(signed int)s.size(); c++)
						{
							unsigned char uc=s[c];
							if (lowerSet[0][uc]>0)
								t.append(1,(char)lowerSetRev[i][lowerSet[0][uc]]);
							else
								t.append(1,(char)uc);
						}

						if (addWord(t,sizeFullDict))
						{
							if (sizeDict%SAMPLE_WORDS_COUNT==0)
								break; 
						}
						else
							continue;
					}
				}
			} // end if (nonlatin)

			if (sizeDict%SAMPLE_WORDS_COUNT==0)
				break;
		}

		if (sizeDict%SAMPLE_WORDS_COUNT_MAX!=0)
			sizeDict=((sizeDict/SAMPLE_WORDS_COUNT_MAX)+1)*SAMPLE_WORDS_COUNT_MAX;

		langCount++;

	}
#ifdef WIN32
	while (_findnext( hFile, &c_file ) == 0);

	_findclose( hFile );

#else
#endif

	return true;
}


void freeNames()
{
	for (int i=0; i<langCount; i++)
		free(langName[i]);
}

int getSourcePath(char* buf, int buf_size)
{
#ifdef WIN32
	int pos;

	pos=GetModuleFileName(NULL,buf,buf_size);

	if (pos>0)
	{	
		for (int i=pos-1; i>=0; i--)
			if (buf[i]=='\\')
			{
				buf[i+1]=0;
				pos=i+1;
				break;
			}
	}
	else
		buf[0]=0;

	return pos;
#else
	buf[0]=0;
	return 0;
#endif
}

int WRT_getFileType(FILE* file,int& recordLen,FILE* english_dictionary)
{
	int dictPathLen;
	unsigned char dictPath[256];

	if (map.size()==0)
	{
		getSourcePath((char*)dictPath,sizeof(dictPath));
		strcat((char*)dictPath,WRT_DICT_DIR);
		dictPathLen=strlen((char*)dictPath);

		readDicts(DICTNAME "*" DICTNAME_EXT,(char*)dictPath,dictPathLen,english_dictionary);
	}

	return WRT_detectFileType(file,10240,5,recordLen);
}

#define SWAP_CASE(c) \
{\
\
	if (swapCase)\
	{\
		if (c!=' ' && c!='\r' && c!='\n' && c!='\'' && c!='"' && c!='\t')\
			swapCase=false;\
\
		if (c>='a' && c<='z')\
		{\
			c-=32; /* toupper(c); */ \
		}\
		else\
		if (c>='A' && c<='Z')\
		{\
			c+=32; /* tolower(c); */ \
		}\
	}\
	else\
	if (c=='.' || c=='!' || c=='?')\
		swapCase=true;\
}



#define ENCODE_GETC(c,file)\
{\
	if (llbckp!=0) \
	{ \
		if (llbckp==32*32) \
		{ \
			c=32; \
			llbckp=32; \
		} \
		else \
		{ \
			c=llbckp; \
			llbckp=0; \
		} \
	} \
	else \
	{ \
		c=getc(file); \
 \
		if (IF_OPTION(OPTION_SPACE_AFTER_EOL) && llast==10) \
		{ \
			if (c==32) \
			{ \
				c=getc(file); \
 \
				if (c==32) \
					llbckp=32*32; \
			} \
			else \
			{ \
				llbckp=c; \
				c=32; \
			} \
		} \
 \
		llast=c; \
	} \
 \
 	fftell++; \
 \
	if (c>127) \
	{ \
		if (IF_OPTION(OPTION_UTF8)) \
		{ \
			if (c!=194 && c!=195) \
			{ \
				TURN_OFF(OPTION_UTF8); \
				restartEnc=true; \
				return; \
			} \
			else \
			{ \
				int c2=fgetc(file); \
				if (c2<128 || c2>191) \
				{ \
					TURN_OFF(OPTION_UTF8); \
					restartEnc=true; \
					return; \
				} \
				else \
					c=c2+(c-194)*64; \
			} \
		} \
	} \
 \
	if (IF_OPTION(OPTION_TO_LOWER_AFTER_PUNCTUATION))\
		SWAP_CASE(c);\
}

#define DECODE_QUEUE(c)\
{\
	if (IF_OPTION(OPTION_SPACE_AFTER_EOL) && llast==10) \
	{ \
		if ((c)!=32) \
		{ \
			WRTd_queue[WRTd_qend++]=32;\
			WRTd_queue[WRTd_qend++]=c;\
		} \
 \
	} \
 	else \
	{ \
		WRTd_queue[WRTd_qend++]=c;\
	} \
}

#define DECODE_PUTC(c)\
{\
	if (IF_OPTION(OPTION_TO_LOWER_AFTER_PUNCTUATION)) \
		SWAP_CASE(c);\
 \
 \
	if (c>127) \
	{ \
		if (IF_OPTION(OPTION_UTF8)) \
		{ \
			DECODE_QUEUE((c >> 6) | 0xc0); \
			c = (c & 0x3f) | 0x80; \
		} \
	} \
 \
	DECODE_QUEUE(c); \
 \
	llast=c; \
 	fftell++; \
}

// preprocess the file
void WRT_encode(FILE* file,FILE* fileout,int fileLen) //,int fileType)
{
	unsigned char s[1024];
	EWordType wordType;
	int last_c,c,next_c;
	int binCount=0;

	preprocessing=0;
	s_size=0;
	last_c=0;
	lastEOL=0;
	EOLType=UNDEFINED;
	wordType=LOWERWORD;
	spaceBefore=NONE;
	initOrder=true;
	lastEOL=0;


	if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && !IF_OPTION(OPTION_USE_DICTIONARY))
	{
		autoSwitch=1<<(31-1); // MaxSignedInt
		preprocessing=autoSwitch;
	}
	else
		if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
			autoSwitch=AUTO_SWITCH*4;
		else
			autoSwitch=AUTO_SWITCH;


	ENCODE_GETC(c,file);
	fftell=0;

	while (!feof(file))
	{
		if (restartEnc)
			return;

		if (fileCorrupted)
			return;

		PRINT_CHARS(("c=%d (%c)\n",c,c));

		if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && preprocessing>0)
		{
			if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(c))
			{
				binCount++;
				preprocessing=autoSwitch;
				PRINT_CHARS(("preprocessing=%d c=%c(%d)\n",preprocessing,c,c));
			}
			else
			{
				preprocessing--;
				PRINT_CHARS(("preprocessing=%d c=%c(%d)\n",preprocessing,c,c));
				if (preprocessing==0)
				{
					initOrder=true;
					if (binCount*100/(fftell+5000)>25)
					{
						autoSwitch=AUTO_SWITCH*16;
						preprocessing=autoSwitch;
					}
				}
			}

			ENCODE_PUTC(c,fileout);
			ENCODE_GETC(c,file);
			continue;
		}

		if (reservedSet[c])
		{
			PRINT_CHARS(("reservedSet[c] c=%d (%c)\n",c,c));

			encodeWord(fileout,s,s_size,wordType);
			s_size=0;
			ENCODE_PUTC(CHAR_ESCAPE,fileout);
			ENCODE_PUTC(c,fileout);

			if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(c))
				preprocessing=autoSwitch;

			last_c=c;
			ENCODE_GETC(c,file);
			continue;
		}


		if ((c>='a' && c<='z') || lowerSet[usedSet][c]>0)
		{
			if (usedSet>0 && lowerSet[usedSet][c]>0)
				c=lowerSetRev[0][lowerSet[usedSet][c]];

			PRINT_CHARS(("a-z c=%d (%c)\n",c,c));

			if (s_size==0)
			{
				wordType=LOWERWORD;
			}
			else
			{
				if (wordType==UPPERWORD)
				{
					encodeWord(fileout,s,s_size,wordType);
					if (IF_OPTION(OPTION_CAPITAL_CONVERSION))
						ENCODE_PUTC(CHAR_LOWERWORD,fileout);

					wordType=LOWERWORD;
					s_size=1;
					s[0]=c;
					last_c=c;
					ENCODE_GETC(c,file);
					continue;
				}
			}

			s[s_size++]=c;
			if (s_size>=(signed int)(sizeof(s)-1))
			{
				encodeWord(fileout,s,s_size,wordType);
				s_size=0;
			}
			last_c=c;
			ENCODE_GETC(c,file);
			continue;
		}


		if ((c<='Z' && c>='A') || upperSet[usedSet][c]>0)
		{
			if (usedSet>0 && upperSet[usedSet][c]>0)
				c=upperSetRev[0][upperSet[usedSet][c]];

			PRINT_CHARS(("A-Z c=%d (%c)\n",c,c));


			if (s_size==0)
			{
				wordType=FIRSTUPPER;
			}
			else
			{
				if (wordType==FIRSTUPPER)
				{
					if (s_size==1)
						wordType=UPPERWORD;
					else
					{
						encodeWord(fileout,s,s_size,wordType);

						wordType=FIRSTUPPER;
						s_size=1;
						s[0]=c;
						last_c=c;
						ENCODE_GETC(c,file);
						continue;
					}
				}
				else if (wordType==LOWERWORD)
				{
					encodeWord(fileout,s,s_size,wordType);

					wordType=FIRSTUPPER;
					s_size=1;
					s[0]=c;
					last_c=c;
					ENCODE_GETC(c,file);
					continue;
				}
			}

			s[s_size++]=c;
			if (s_size>=(signed int)(sizeof(s)-1))
			{
				encodeWord(fileout,s,s_size,wordType);
				s_size=0;
			}

			last_c=c;
			ENCODE_GETC(c,file);
			continue;
		}

		encodeWord(fileout,s,s_size,wordType);

		s_size=0;

		PRINT_CHARS(("other c=%d\n",c));

		ENCODE_GETC(next_c,file);

#if USE_EOLC
		if (IF_OPTION(OPTION_EOL_CODING) && (c==32 || c==10 || (c==13 && next_c==10)))
		{
			if (initOrder)
			{
				initOrder=false;
				INIT_ORDER1;
			}


			if (c==13)
			{
				if (EOLType==CRLF || EOLType==UNDEFINED)
				{
					c=next_c;
					ENCODE_GETC(next_c,file);
					lastEOL++;

						last_c=ContextEncode(last_c,c,next_c,fftell-lastEOL+(next_c<0?1:0));


					if (EOLType==UNDEFINED && last_c!=c)
					{
						EOLType=CRLF;
						EncodeEOLformat(EOLType);
					}

					lastEOL=fftell;

					if (last_c==10)
					{
						ENCODE_PUTC(CHAR_CR_LF,fileout);
					}
					else
					{
						if ((IF_OPTION(OPTION_SPACELESS_WORDS)) && last_c==' ')
							spaceBefore=SPACE;
						else
							ENCODE_PUTC(last_c,fileout);
					}
				}
				else
					ENCODE_PUTC(c,fileout);
			}
			else
			{
				if (c==10 && EOLType==CRLF)
				{ ENCODE_PUTC(c,fileout); }
				else
				{
					if (IF_OPTION(OPTION_EOL_CODING))
					{		
							last_c=ContextEncode(last_c,c,next_c,fftell-lastEOL+(next_c<0?1:0));
					}
					else
						last_c=c;

					if (EOLType==UNDEFINED && last_c!=c)
					{
						EOLType=LF;
						EncodeEOLformat(EOLType);
					}

					if ((IF_OPTION(OPTION_SPACELESS_WORDS)) && last_c==' ')
						spaceBefore=SPACE;
					else
						ENCODE_PUTC(last_c,fileout);
				}

				if (c==10)
				{
					lastEOL=fftell;
				}
			}
		}
		else
#endif
		{
			if (c!=EOF)
			{
				if ((IF_OPTION(OPTION_WORD_SURROROUNDING_MODELING)) && c>' ')
				{
					if ((last_c>='a' && last_c<='z') || (last_c>='A' && last_c<='Z'))
					{
						if (fftell<fileLen)
						{
							ENCODE_PUTC(' ',fileout);
							ENCODE_PUTC(CHAR_PUNCTUATION,fileout);
							ENCODE_PUTC(c,fileout);
						}
						else
							ENCODE_PUTC(c,fileout);
					}
					else 
						if (next_c>='a' && next_c<='z')
						{
							ENCODE_PUTC(c,fileout);
							ENCODE_PUTC(CHAR_LOWERWORD,fileout);
							ENCODE_PUTC(' ',fileout);
						}
						else
							ENCODE_PUTC(c,fileout);
				}
				else 
					if ((IF_OPTION(OPTION_SPACELESS_WORDS)) && c==' ')
						spaceBefore=SPACE;
					else
						ENCODE_PUTC(c,fileout);
			}
		}

		last_c=c;
		c=next_c;
	}

	encodeWord(fileout,s,s_size,wordType);
	s_size=0;
}


inline void hook_putc(int c)
{
	if (bufferedChar<0 && c==' ') 
	{
		bufferedChar=c;
		return;
	}

	if (bufferedChar>=0)
	{	
#if USE_EOLC
		if (IF_OPTION(OPTION_EOL_CODING))
		{
			if (initOrder)
			{
				initOrder=false;
				INIT_ORDER1;

			}

				bufferedChar=ContextDecode(lastChar,c,fftell-lastEOL);

			if (bufferedChar==10)
			{
				if (EOLType==UNDEFINED)
					EOLType=DecodeEOLformat();
				if (EOLType==CRLF)
				{
					lastChar=13;
					DECODE_PUTC(lastChar);
				}

				lastEOL=fftell;
			}
		}
#endif
		DECODE_PUTC(bufferedChar);

		if (c==' ')
		{
			lastChar=bufferedChar;
			bufferedChar=c;
			return;
		}

		bufferedChar=-1;
	}
	
	if (c==10)
		lastEOL=fftell;


	lastChar=c;


	if (c==EOF)
		return;

	DECODE_PUTC(c);
	return;
}


int readEOLstream(FILE* file)
{
	unsigned int EOLstream_len=0;

	fseek(file, -4, SEEK_END );

	for (int i=0; i<4; i++)
	    EOLstream_len=EOLstream_len*256+fgetc(file);

#if USE_EOLC
	if (IF_OPTION(OPTION_EOL_CODING) && EOLstream_len>0 && EOLstream_len<fileLen)
	{
		fseek(file, -((int)EOLstream_len)-4, SEEK_END );

		EOLcount=0;
		INIT_ORDER1;
		coder.StartDecode(file,EOLstream_len);
	}
#endif
	return EOLstream_len;
}



void writeEOLstream(FILE* fileout)
{
	unsigned int EOLstream_len;

	EOLstream_len=0;

#if USE_EOLC
	if (IF_OPTION(OPTION_EOL_CODING) && EOLcount>0)
	{
		EOLstream_len=coder.FinishEncode(fileout);
	}
#endif

	fprintf(fileout, "%c%c%c%c", EOLstream_len>>24, EOLstream_len>>16, EOLstream_len>>8, EOLstream_len);
}

inline void WRT_decode(FILE* file)
{
		PRINT_CHARS(("c=%d (%c)\n",WRTd_c,WRTd_c));

		if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && preprocessing>0)
		{
			DECODE_PUTC(WRTd_c);

			if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(WRTd_c))
			{
				preprocessing=autoSwitch;
				WRTd_binCount++;
			}
			else
			{
				preprocessing--;
				PRINT_CHARS(("preprocessing=%d c=%c(%d)\n",preprocessing,WRTd_c,WRTd_c));
				if (preprocessing==0)
				{
					initOrder=true;
					if (WRTd_binCount*100/(fftell+5000)>25)
					{
						autoSwitch=AUTO_SWITCH*16;
						preprocessing=autoSwitch;
					}
				}
			}

			DECODE_GETC(WRTd_c,file);
			return;
		}


		if (addSymbols[WRTd_c] && IF_OPTION(OPTION_USE_DICTIONARY))
		{
			PRINT_CHARS(("addSymbols[c] && IF_OPTION(OPTION_USE_DICTIONARY) upperWord=%d\n",upperWord));
	
			if (spaceBefore==SPACE)
			{
				if (upperWord==FORCE)
					upperWord=UTRUE;
				else
					upperWord=UFALSE;
			}

				s_size=decodeCodeWord(file,WRTd_s,WRTd_c);

			if (WRTd_upper)
			{
				WRTd_upper=false;
				WRTd_s[0]=TOUPPER(WRTd_s[0]);
			}

			if (upperWord!=UFALSE)
				toUpper(WRTd_s,s_size);

			if (IF_OPTION(OPTION_SPACELESS_WORDS))
			{
				if (spaceBefore==SPACE)
					hook_putc(' '); 
				else
					spaceBefore=SPACE;
			}

			int i;
			PRINT_CHARS(("word="));
			for (i=0; i<s_size; i++)
				PRINT_CHARS(("%c",WRTd_s[i]));
			PRINT_CHARS((" upperWord=%d\n",upperWord));


			for (i=0; i<s_size; i++)
			{
				ORIGINAL_CHARSET(WRTd_s[i]);
				hook_putc(WRTd_s[i]);
			}

			DECODE_GETC(WRTd_c,file);
			return;
		}


		if (reservedSet[WRTd_c])
		{
			PRINT_CHARS(("reservedSet[%d] OPTION_SPACELESS_WORDS=%d\n",WRTd_c,IF_OPTION(OPTION_SPACELESS_WORDS)));

			if (WRTd_c==CHAR_ESCAPE)
			{
				WRTd_upper=false;
				upperWord=UFALSE;

				DECODE_GETC(WRTd_c,file);
				PRINT_CHARS(("c==CHAR_ESCAPE, next=%x\n",WRTd_c));
				hook_putc(WRTd_c);

				if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && COND_BIN_FILTER(WRTd_c))
					preprocessing=autoSwitch;

				DECODE_GETC(WRTd_c,file);
				return;
			}

			if (WRTd_c==CHAR_NOSPACE)
			{
				PRINT_CHARS(("c==CHAR_NOSPACE\n"));

				if (upperWord==FORCE)
					upperWord=UTRUE;

				DECODE_GETC(WRTd_c,file);
				spaceBefore=NONE;
				return;
			}

			if (WRTd_c==CHAR_CR_LF)
			{
				PRINT_CHARS(("c==CHAR_CR_LF\n"));

				hook_putc('\r');
				WRTd_c='\n';
			}

			if (WRTd_c==CHAR_FIRSTUPPER)
			{
				PRINT_CHARS(("c==CHAR_FIRSTUPPER\n"));

				if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
				{			
					DECODE_GETC(WRTd_c,file); // skip space
				}
				WRTd_upper=true;
				upperWord=UFALSE;
				DECODE_GETC(WRTd_c,file);
				PRINT_CHARS(("c==CHAR_FIRSTUPPER WRTd_c=%d\n",WRTd_c));
				return;
			}
			
			if (WRTd_c==CHAR_UPPERWORD)
			{
				PRINT_CHARS(("c==CHAR_UPPERWORD\n"));

				if (IF_OPTION(OPTION_SPACE_AFTER_CC_FLAG))
				{
					DECODE_GETC(WRTd_c,file); // skip space
				}
				upperWord=FORCE;
				DECODE_GETC(WRTd_c,file);
				return;
			}
			
			if (WRTd_c==CHAR_LOWERWORD)
			{
				PRINT_CHARS(("c==CHAR_LOWERWORD\n"));

				WRTd_upper=false;
				upperWord=UFALSE;
				DECODE_GETC(WRTd_c,file);
				if (WRTd_c==' ') // skip space
				{
					DECODE_GETC(WRTd_c,file);
				}
				return;
			}
		}

 		PRINT_CHARS(("other c=%d (%d %d)\n",WRTd_c,fftell+((bufferedChar>=0)?1:0),originalFileLen));

		if (upperWord!=UFALSE)
		{
			if (upperWord==FORCE)
				upperWord=UTRUE;

			if ((WRTd_c>='a' && WRTd_c<='z') || lowerSet[usedSet][WRTd_c]>0)
				WRTd_c=TOUPPER_SET(WRTd_c);
			else
				upperWord=UFALSE;
		}
		else
		if (WRTd_upper==true)
		{
			WRTd_upper=false;
			WRTd_c=TOUPPER_SET(WRTd_c);
		}

		hook_putc(WRTd_c);

		DECODE_GETC(WRTd_c,file);
		return;
}

void WRT_start_encoding(FILE* file,FILE* fileout,unsigned int fileLen,bool type_detected,FILE* english_dictionary)
{
  unsigned int start_pos = ftell(fileout);
	int i,c,c2,recordLen=0,dictPathLen;
	unsigned char s[256];
	unsigned char t[256];
	unsigned char dictPath[256];
	s[0]=0;
	t[0]=0;

	llbckp=0;
	swapCase=false;
	usedSet=0;

	getSourcePath((char*)dictPath,sizeof(dictPath));
	strcat((char*)dictPath,WRT_DICT_DIR);
	dictPathLen=strlen((char*)dictPath);

	if (!type_detected)
		WRT_getFileType(file,recordLen,english_dictionary);

	if (IF_OPTION(OPTION_USE_DICTIONARY) && longDictLen>0)
		memcpy(s,langName[longDict],strlen((const char*)langName[longDict])+1);
	else
		longDictLen=0;
	
	if (IF_OPTION(OPTION_USE_DICTIONARY) && shortDictLen>0)
		memcpy(t,langName[shortDict],strlen((const char*)langName[shortDict])+1);
	else
		shortDictLen=0;

	if (dictPathLen>0)
	{
		dictPath[dictPathLen]=0;
		strcat((char*)dictPath,(char*)s);
		strcpy((char*)s,(char*)dictPath);

		dictPath[dictPathLen]=0;
		strcat((char*)dictPath,(char*)t);
		strcpy((char*)t,(char*)dictPath);
	}


restart:

	int pos=ftell(fileout);
	fprintf(fileout, "%c%c%c%c", 0,0,0,0);

	WRT_get_options(c,c2); // before initialize
	putc(c,fileout);
	putc(c2,fileout);

	WRT_deinitialize();
	if (!initialize((longDictLen<=0)?NULL:s,(shortDictLen<=0)?NULL:t,true,english_dictionary))
		return;

	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
		putc(usedSet,fileout);
		putc(shortDictLen,fileout);
		c=strlen(SHORT_DICTNAME);
		for (i=c; i<shortDictLen+c; i++)
			fputc(langName[shortDict][i],fileout);
		
		putc(longDictLen,fileout);
		c=strlen(DICTNAME);
		for (i=c; i<longDictLen+c; i++)
			fputc(langName[longDict][i],fileout);
	}

  {
#if USE_EOLC
		    EOLcount=0;
		    if (IF_OPTION(OPTION_EOL_CODING)) {
			INIT_ORDER1;
			coder.StartEncode();
		    }
#endif	

			WRT_encode(file,fileout,fileLen); 
			if (restartEnc)
			{
				restartEnc=false;
				fseek(fileout, pos, SEEK_SET );
				fseek(file, 0, SEEK_SET );
				llbckp=0;
				swapCase=false;
				goto restart;
			}

			unsigned int fileLen=ftell(fileout) + 4;
      unsigned int store = fileLen - start_pos;
			fseek(fileout, pos, SEEK_SET );
			fprintf(fileout, "%c%c%c%c", store>>24, store>>16, store>>8, store);
			fseek(fileout, fileLen-4, SEEK_SET );
		}
}


void WRT_start_decoding(FILE* file,FILE* fileout,int header,FILE* english_dictionary)
{
	int i,c,c2,dictPathLen;
	unsigned char s[256];
	unsigned char t[256];
	unsigned char dictPath[256];
	unsigned int fileLen;
	s[0]=0;
	t[0]=0;

	llbckp=0;
	swapCase=false;
	usedSet=0;
	WRTd_binCount=0;

	for (i=0, fileLen=0; i<4; i++)
	    fileLen=fileLen*256+fgetc(file);

	i=0;
	c=getc(file);
	c2=getc(file);
	
	preprocType=(EPreprocessType)(c%4); // { LZ77, BWT, PPM, PAQ };
	
	defaultSettings(0,NULL); // after setting preprocType 
	
	WRT_set_options(c,c2);

	if (IF_OPTION(OPTION_USE_DICTIONARY))
	{
		usedSet=getc(file);
		shortDictLen=getc(file);
		c=strlen(SHORT_DICTNAME);
		memcpy(t,SHORT_DICTNAME,c);
		for (i=c; i<shortDictLen+c; i++)
			t[i]=getc(file);
		c=strlen(DICTNAME_EXT);
		memcpy(t+i,DICTNAME_EXT,c);
		i+=c;
		t[i]=0;
		
		longDictLen=getc(file);
		c=strlen(DICTNAME);
		memcpy(s,DICTNAME,c);
		for (i=c; i<longDictLen+c; i++)
			s[i]=getc(file);
		c=strlen(DICTNAME_EXT);
		memcpy(s+i,DICTNAME_EXT,c);
		i+=c;
		s[i]=0;
		
		i=longDictLen+1+shortDictLen+1+(IF_OPTION(OPTION_USE_DICTIONARY)?1:0); // usedSet
	}
	else
	{
		longDictLen=0;
		shortDictLen=0;
	}
	header+=4;
	i+=2+header; // WRT4
	
	getSourcePath((char*)dictPath,sizeof(dictPath));
	strcat((char*)dictPath,WRT_DICT_DIR);
	dictPathLen=strlen((char*)dictPath);

	if (dictPathLen>0)
	{
		dictPath[dictPathLen]=0;
		strcat((char*)dictPath,(char*)s);
		strcpy((char*)s,(char*)dictPath);

		dictPath[dictPathLen]=0;
		strcat((char*)dictPath,(char*)t);
		strcpy((char*)t,(char*)dictPath);
	}

	WRT_deinitialize();

	if (!initialize((longDictLen<=0)?NULL:s,(shortDictLen<=0)?NULL:t,false,english_dictionary))
		return;

		{
#if USE_EOLC
			int EOLlen= readEOLstream(file)+4; // +sizeof(EOLlen)
#else
			int EOLlen= 4;
#endif

			fseek(file, i, SEEK_SET ); // skip "WRTx" header
			EOLlen+=i;	// header + fileLen

#ifdef POWERED_BY_PAQ
			WRTd_filter->reads+=i;
#endif

			originalFileLen=fileLen-EOLlen;
			bufferedChar=-1;
			lastChar=0;
			fftell=0;
			fftelld=0;
			WRTd_upper=false;
			upperWord=UFALSE;
			preprocessing=0;
			s_size=0;
			initOrder=true;
			lastEOL=-1;
			EOLType=UNDEFINED;
			
			
			if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER) && !IF_OPTION(OPTION_USE_DICTIONARY))
			{
				autoSwitch=1<<(31-1); // MaxSignedInt
				preprocessing=autoSwitch;
			}
			else
				if (!IF_OPTION(OPTION_NORMAL_TEXT_FILTER))
					autoSwitch=AUTO_SWITCH*4;
				else
					autoSwitch=AUTO_SWITCH;
				
				if (IF_OPTION(OPTION_SPACELESS_WORDS))
					spaceBefore=SPACE;
				else
					spaceBefore=NONE;
				
				
				DECODE_GETC(WRTd_c,file);
				PRINT_CHARS(("WRT_start_decoding WRTd_c=%d ftell=%d\n",WRTd_c,ftell(file)));
		} 
}

void WRT_prepare_decoding()
{
	WRTd_type=0;
}

int WRT_decode_char(FILE* file,FILE* fileout,int header,FILE* english_dictionary)
{
	switch (WRTd_type)
	{
		default:
		case 0:
			WRT_start_decoding(file,fileout,header,english_dictionary);
			WRTd_qstart=WRTd_qend=0;
			WRTd_type=1;
		case 1:
			if (WRTd_c!=EOF)
			{
				while (WRTd_qstart>=WRTd_qend && WRTd_c!=EOF)
				{
					WRTd_qstart=WRTd_qend=0;
					WRT_decode(file);
					if (fileCorrupted)
						WRTd_type=2;
				}

				if (WRTd_qstart<WRTd_qend)
					return WRTd_queue[WRTd_qstart++];
			}
			hook_putc(EOF);
			WRTd_type=2;
		case 2:
			if (WRTd_qstart<WRTd_qend)
				return WRTd_queue[WRTd_qstart++];
			else
				return -1;
	}
}

}; // end class 
