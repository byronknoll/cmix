
#define COMP_INTRO_END_LINE 29
#define COMP_MAIN_END_LINE  13146932
#define COMP_CODA_END_LINE  13147025

#define DECOMP_MAIN_END_LINE  13146905
#define DECOMP_INTRO_END_LINE 13146934 
#define DECOMP_CODA_END_LINE  13147027

void split4Comp(char const *enwik9_filename) {
  std::ifstream ifile(enwik9_filename); 
  std::ofstream ofile1(".intro"); 
  std::ofstream ofile2(".main"); 
  std::ofstream ofile3(".coda"); 

  int line_count = 0;

  std::string s;
  while (std::getline(ifile, s)) {
    if (line_count < COMP_INTRO_END_LINE) {
      ofile1 << s << std::endl;
    } else if (line_count < COMP_MAIN_END_LINE) {
        ofile2 << s << std::endl;
      } else if (line_count < COMP_CODA_END_LINE) {
          ofile3 << s << std::endl;
          } else {
            ofile3 << s; 
          }
    line_count++;
  }

  ofile1.close();
  ofile2.close();
  ofile3.close();
}

void split4Decomp( const char* inpnam ) {
  std::ifstream ifile(inpnam); 
  std::ofstream ofile1(".intro_decomp"); 
  std::ofstream ofile2(".main_decomp"); 
  std::ofstream ofile3(".coda_decomp"); 

  int line_count = 0;

  std::string s;
  while (std::getline(ifile, s)) {
    if (line_count < 13146905) {
      ofile2 << s << std::endl;
    } else if (line_count < 13146934) {
        ofile1 << s << std::endl;
      } else if (line_count < 13147027) {
          ofile3 << s << std::endl;
          } else {
            ofile3 << s; 
          }
    line_count++;
  }

  ofile1.close();
  ofile2.close();
  ofile3.close();
}
