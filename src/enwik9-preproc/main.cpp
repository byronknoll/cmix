// This code comes from STARLIT: https://github.com/amargaritov/starlit
// It was created by Artemiy Margaritov and Alexander Rhatushnyak.
// A standalone version of the preprocessor was released by Eugene Shelwien:
// https://encode.su/threads/3635-Hutter-Prize-Entry-quot-STARLIT-quot-Open-For-Comments?p=70028&viewfull=1#post70028

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define ifstream ifstream1
#define ofstream ofstream1

#include <cstdint>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>

#undef ifstream
#undef ofstream

namespace std {

  struct ifstream : ifstream1 {
    ifstream( string s ) : ifstream1(s,ios::binary) {}
    ifstream( string s, ios_base::openmode x ) : ifstream1(s,x|ios::binary) {}
  };

  struct ofstream : ofstream1 {
    ofstream( string s ) : ofstream1(s,ios::binary) {}
    ofstream( string s, ios_base::openmode x ) : ofstream1(s,x|ios::binary) {}
  };

}

#include "misc.h"
#include "phda9_preprocess.h"
#include "article_reorder.h"

int main( int argc, char** argv ) {

  if( argv[1][0]!='d' ) {
    split4Comp( argv[2] );
    // change the order of articles in the input 
    reorder();

    // apply phda9 preprocessor
    phda9_prepr();

    // merge all input parts after preprocessing
    cat(".main_phda9prepr", ".intro", "un1");
    cat("un1", ".coda", ".ready4cmix");

  } else {
    split4Decomp( argv[2] );
    // apply phda9 preprocessor
    phda9_resto();

    // change the order of articles in the input 
    sort();

    // merge all input parts after preprocessing
    cat(".intro_decomp", ".main_decomp_restored_sorted", "un1_d");
    cat("un1_d", ".coda_decomp", "enwik9_uncompressed");
  }
}
