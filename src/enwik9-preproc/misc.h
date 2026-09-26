
#define COMP_INTRO_END_LINE 29
#define COMP_MAIN_END_LINE  13146932
#define COMP_CODA_END_LINE  13147025

#define DECOMP_MAIN_END_LINE  13146905
#define DECOMP_INTRO_END_LINE 13146934 
#define DECOMP_CODA_END_LINE  13147027

void split4Comp(char const *enwik9_filename) {
  FILE* ifile = fopen(enwik9_filename, "rb");
  if (!ifile) return;
  FILE* ofile1 = fopen(".intro", "wb");
  FILE* ofile2 = fopen(".main", "wb");
  FILE* ofile3 = fopen(".coda", "wb");  
  if (!ofile1 || !ofile2 || !ofile3) {
    if (ofile1) fclose(ofile1);
    if (ofile2) fclose(ofile2);
    if (ofile3) fclose(ofile3);
    fclose(ifile);
    return;
  }
  int line_count = 0;
  
  constexpr size_t BUF_SIZE = 65536;
  char buf[BUF_SIZE];
  size_t n;
  FILE* cur_file = (line_count < COMP_INTRO_END_LINE) ? ofile1 :
                   (line_count < COMP_MAIN_END_LINE)  ? ofile2 : ofile3;

  while ((n = fread(buf, 1, BUF_SIZE, ifile)) > 0) {
    size_t start = 0;
    for (size_t i = 0; i < n; ++i) {
      if (buf[i] == '\n') {
        line_count++;
        FILE* next_file = (line_count < COMP_INTRO_END_LINE) ? ofile1 :
                          (line_count < COMP_MAIN_END_LINE)  ? ofile2 : ofile3;
        if (next_file != cur_file) {
          fwrite(buf + start, 1, i + 1 - start, cur_file);
          start = i + 1;
          cur_file = next_file;
        }
      }
    }
    if (start < n) {
      fwrite(buf + start, 1, n - start, cur_file);
    }
  }
  fclose(ifile);
  fclose(ofile1);
  fclose(ofile2);
  fclose(ofile3);
}

void split4Decomp( const char* inpnam ) {
  FILE* ifile = fopen(inpnam, "rb");
  if (!ifile) return;
  FILE* ofile1 = fopen(".intro_decomp", "wb");
  FILE* ofile2 = fopen(".main_decomp", "wb");
  FILE* ofile3 = fopen(".coda_decomp", "wb");  
  if (!ofile1 || !ofile2 || !ofile3) {
    if (ofile1) fclose(ofile1);
    if (ofile2) fclose(ofile2);
    if (ofile3) fclose(ofile3);
    fclose(ifile);
    return;
  }
  int line_count = 0;
  
  constexpr size_t BUF_SIZE = 65536;
  char buf[BUF_SIZE];
  size_t n;
  FILE* cur_file = (line_count < 13146906) ? ofile2 :
                   (line_count < 13146935) ? ofile1 : ofile3;

  while ((n = fread(buf, 1, BUF_SIZE, ifile)) > 0) {
    size_t start = 0;
    for (size_t i = 0; i < n; ++i) {
      if (buf[i] == '\n') {
        line_count++;
        FILE* next_file = (line_count < 13146906) ? ofile2 :
                          (line_count < 13146935) ? ofile1 : ofile3;
        if (next_file != cur_file) {
          fwrite(buf + start, 1, i + 1 - start, cur_file);
          start = i + 1;
          cur_file = next_file;
        }
      }
    }
    if (start < n) {
      fwrite(buf + start, 1, n - start, cur_file);
    }
  }
  fclose(ifile);
  fclose(ofile1);
  fclose(ofile2);
  fclose(ofile3);
}
