
#define P3_INPUT_SIZE 999988944
#define P4_INPUT_SIZE 985154324 // out3 size
#define TMP1A_SIZE 937793245
#define TMP1B_SIZE 16193156
#define TMP2A_SIZE 929737295
#define TMP2B_SIZE 24249088
#define OUT5_SIZE  980080701
#define R3_INPUT_SIZE 984215015 // size of out4 = out5d 
#define R3_OUTPUT_ARR_SIZE 1100000000 // set to maximum possible
#define R4_INPUT_SIZE      999049635 // size of out3d 
#define R4_OUTPUT_ARR_SIZE 1100000000 // set to maximum possible


int prepr1(char const* argv[]) {
  char s[8192], *ps;
  FILE *sf, *t1, *t2;
  int c, i = 0, mi = 0, lnu = 0, j = 0, f = 0, b1 = 0, lc = 0;
  sf = fopen(argv[0], "rb");
  t1 = fopen(argv[1], "wb");
  t2 = fopen(argv[2], "wb");

  do {
    fgets(s, 8192, sf);
    ++lnu;
    if (feof(sf)) {
      break;
    }
    j = strlen(s);
    for (i = 0; i < j; ++i) {
      if (*(int*)&s[i] == 'xet<') {
        b1 = lnu;
      }
    }

    if (f == 0) {
      if (s[0] == '[' && s[1] == '[') { // if starts with [[
        ps = s + (s[2] == ':' ? 1 : 0); // ps = s + 1 if s starts with [[:
        for (c = 2; c < j; ++c) {
          if ((s[c] < 'a' && s[c] != '-') || s[c] > 'z') {
            break;
          }
        }
        if (c < j && s[c] == ':') {
          if (memcmp(&ps[2], "http:", 5) == 0) { // if [[http:"...
            goto TOFILE1;
          }
          if (memcmp(&ps[2], "user:", 5) == 0) {
            goto TOFILE1;
          }
          if (memcmp(&ps[2], "media:", 6) == 0) {
            goto TOFILE1;
          }
          if (memcmp(&ps[3], "mage:", 5) == 0) {
            goto TOFILE1;
          }
          if (memcmp(&ps[3], "ategory:", 8) == 0) {
            goto TOFILE1;
          }
          if (memcmp(&ps[2], "fr:Wikipédia:Aide]]", 20) == 0) {
            goto TOFILE1;
          }
          if (memcmp(&ps[2], "de:Boogie Down Produ", 20) == 0) {
            goto TOFILE1;
          }
          if (memcmp(&ps[2], "da:Wikipedia:Hvordan", 20) == 0) {
            goto TOFILE1;
          }
          if (memcmp(&ps[2], "sv:Indiska musikinstrument", 26) == 0) {
            goto TOFILE1;
          }
          if (lnu - b1 < 4) {
            goto TOFILE1;
          }
          f = 1; // f = 1 means that we found where text section begins and waiting for it end  <text> ... (we are here)  </text>
          lc = 0;
          goto SEARCHEND;
        }
      }
    TOFILE1:
      fputs(s, t1);
      continue;
    }  // if (f==0)

  SEARCHEND: //search for end of text section
    ++lc;
    for (i = 0; i < j; i++) {
      if (*(int*)&s[i] == 'et/<') { //found end of text section
        f = 0;
      }
    }
    fputs(s, t2);
    ++mi;
  } while (!feof(sf));

  fclose(t2);
  fclose(t1);
  fclose(sf);

  return 0;
}


int prepr2(char const *argv[]) {
  char s[8192];
  FILE *sf, *t1, *t2;
  int i, line_length, f = 0, lastID = 0;
  sf = fopen(argv[0], "rb");
  t1 = fopen(argv[1], "wb");
  t2 = fopen(argv[2], "wb");

  int num_iteration = 0;

  do {
    fgets(s, 8192, sf);
    line_length= strlen(s);

    if (f == 2) {
      int curID = atoi(&s[8]);
      if (*(int *)&s[4] != '>di<') {
        if (*(int *)&s[4] != 'aeh<') {
          f = 0;
          goto normal;
        }
      }
      fprintf(t2, ">%d%c", curID - lastID, 10);
      lastID = curID;
      f = 1;
      continue;
    }

normal:
    if (f) {
      if (*(int *)&s[6] == 'mit<') { // process timestamp
        int year = atoi(&s[17]);
        int month = atoi(&s[22]);
        int day = atoi(&s[25]);
        int hour = atoi(&s[28]);
        int minute = atoi(&s[31]);
        int second = atoi(&s[34]);
        fprintf(t2, "timestamp>%d%d:%d%c", year - 2002, month * 31 + day - 32,
                hour * 3600 + minute * 60 + second, 10);
        continue;
      }
      char *p = strchr(s, '>');
      if (p) {
        p = strchr(p + 1, '<');
        if (p) {
          p[0] = 10, p[1] = 0;
        }
      }
      if (f == 3) { // expect contributor
        if (*(int *)&s[6] == 'oc/<') {
          fputs(&s[7], t2);
        } else { // expect something of 9 symbols
          fputs(&s[9], t2);
        }
      } else { // f != 3
        if (*(int *)&s[4] == 'ver<' || *(int *)&s[4] == 'ser<') {
          fputs(&s[5], t2);
        } else {
          fputs(&s[7], t2);
        }
        if (*(int *)&s[6] == 'noc<') {
          f = 3;
        }
      }
    } else {
      fputs(s, t1);
    }

    for (i = 0; i < line_length; i++) {
      if (*(int *)&s[i] == 'it/<' && *(int *)&s[i + 4] == '>elt' && 
          *(int *)s == '    ') {
        f = 2; // read line in for of <title>XXX</title>
      }
    }
    for (i = 0; i < line_length; i++) {
      if (*(int *)&s[i] == 'oc/<' && *(int *)&s[i + 4] == 'irtn') {
        f = 0;
      }
    }
    num_iteration++;
  } while (!feof(sf));

  fclose(t2);
  fclose(t1);
  fclose(sf);

  return 0;
}


int prepr3(char const *argv[]) {
  FILE *sf, *t2;
  int j, k;
  unsigned char *p1 = (unsigned char *)malloc(P3_INPUT_SIZE), *p0 = p1,
                *p3 = (unsigned char *)malloc(P3_INPUT_SIZE), *p4 = p3;

  sf = fopen(argv[0], "rb");
  if (sf == 0) {
  }
  fread(p1, P3_INPUT_SIZE, 1, sf);

  do {
    j = *p1++; //read symbol from input 
    *p4++ = j; //write symbol to output

    if (j == '&') {
      k = *(int *)p1; //read following 4 symbols from input
      if (k == ';pma') { // &amp; -> &&
        *p4++ = '&', p1 += 4;
      } else if (k == 'touq' && *(p1 + 4) == ';') { //&quot; -> &"
        *p4++ = '"', p1 += 5;
      } else {
        k = k * 256 + ' ';
        if (k == ';tl ') {
          *p4++ = '<', p1 += 3;
        } else if (k == ';tg ') {
          *p4++ = '>', p1 += 3;
        }
      }
    }
  } while (p1 != p0 + P3_INPUT_SIZE);

  t2 = fopen(argv[1], "wb");
  fwrite(p3, p4 - p3, 1, t2);

  fclose(t2);
  fclose(sf);

  free(p0);
  free(p3);

  return 0;
}


int prepr4(char const *argv[]) {
  FILE *sf, *t2;
  int j, k;
  char *p1 = (char *)malloc(P4_INPUT_SIZE), *p0 = p1,
                *p3 = (char *)malloc(P4_INPUT_SIZE), *p4 = p3;

  sf = fopen(argv[0], "rb");
  fread(p1, P4_INPUT_SIZE, 1, sf);

  do {
    j = *p1++;
    *p4++ = j;

    if (j == '&' && *(p1 - 2) == '&') { // if curent is & and the one before prev is & ...&&..
      k = *(int *)p1;
      if (k == 'touq' && *(p1 + 4) == ';') {
        *p4++ = '"', p1 += 5;
      } else if (k == 'psbn' && *(p1 + 4) == ';') {
        *p4++ = '}', p1 += 5;
      } else if (k == 'sadn' && *(p1 + 4) == 'h' && *(p1 + 5) == ';') {
        *p4++ = '@', p1 += 6;
      } else if (k == 'sadm' && *(p1 + 4) == 'h' && *(p1 + 5) == ';') {
        *p4++ = '`', p1 += 6;
      } else {
        k = k * 256 + ' ';
        if (k == ';tl ') {
          *p4++ = '<', p1 += 3;
        } else if (k == ';tg ') {
          *p4++ = '>', p1 += 3;
        }
      }
    }
  } while (p1 != p0 + P4_INPUT_SIZE);

  t2 = fopen(argv[1], "wb");
  fwrite(p3, p4 - p3, 1, t2);

  fclose(t2);
  fclose(sf);

  free(p0);
  free(p3);

  return 0;
}

int prepr5(char const* argv[])
{
  char s[16384];
  FILE *sf, *t1;
  int mi=0, tf=0;
  sf = fopen(argv[0],"rb");
  t1 = fopen(argv[1],"wb");

  while(1) {
    fgets(s, 16384, sf);
    char *p = strstr(s, "<text "), *w;
    if (p)  { tf=1, p = strchr(p, '>'); assert(p);  if(p[-1]=='/') tf=0;  ++p; }
    else  p = s;
    if (strstr(s, "</text>"))  tf=0;
    if (tf) {
        for(w=p; *p!=0; ++p)  {
          if (*p=='"' || *p=='<' || *p=='>')  { assert(p[-1]=='&');  --w; }
          *w++=*p;
        }
        *w = 0;
    }
    fputs(s, t1);
    ++mi;
    if (feof(sf))  break;
  }

  fclose(t1);
  fclose(sf);

  return 0;
}


int prepr6(char const* argv[])
{
  char s[16384], z[16384];
  FILE *sf, *t1;
  int mi=0;
  sf = fopen(argv[0],"rb");
  t1 = fopen(argv[1],"wb");

#define PROCESS(sym, src, dst, CONDITION) \
  {\
    char *t,  *p = src,  *q = dst,  *end = p + strlen(src);\
    while ( (t=strchr(p, sym)) != 0) {   \
        memcpy(q, p, t-p);  q+=t-p;      \
        int count = 0;                   \
        while(*t++ == sym)  ++count;     \
        if ((CONDITION) && (count==1 || count==2))  count = 3-count;\
        memset(q, sym, count);  q+=count;\
        p = t-1;\
    }\
    memcpy(q, p, end-p);  q[end-p]=0; \
  }

  while(1) {
    fgets(s, 16384, sf);
    if (feof(sf))  break;
    PROCESS('{', s, z, 1)
    PROCESS('}', z, s, 1)
    PROCESS('[', s, z, 1)
    PROCESS(']', z, s, 1)
    PROCESS('&', s, z, 1)
    fputs(z, t1);
    ++mi;
  }

  fclose(t1);
  fclose(sf);

  return 0;
}


bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
      return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

bool sed(std::string old_pattern, std::string new_pattern, char const * filename_from, char const * filename_to) {
  std::ifstream ifile(filename_from, std::ifstream::in); 
  std::ofstream ofile(filename_to); 

  std::string s;
  while (std::getline(ifile, s)) {
    replace(s, old_pattern, new_pattern);
    ofile << s << std::endl;
  }
  ofile.close();
  return true;
}

bool cat(char const * filename_from1, char const * filename_from2, char const * filename_to) {
  std::ifstream ifile1(filename_from1); 
  std::ifstream ifile2(filename_from2); 
  std::ofstream ofile(filename_to); 

  char ch;
  while (ifile1 >> std::noskipws >> ch) {
    ofile << ch; 
  }

  while (ifile2 >> std::noskipws >> ch) {
    ofile << ch; 
  }
  ofile.close();
  return true;
}

int phda9_prepr() {
  sed("&lt;/title&gt;", "&lt;/tiqqqtle&gt;", ".main_reordered", "out8");
  sed("&amp;<", "&amp;qqq<", "out8", "out7");

  {
    char const* argv[] = {"out7", "out3"};
    prepr3(argv);
  }
  {
    char const* argv[] = {"out3", "out4"};
    prepr4(argv);
  }
  {
    char const* argv[] = {"out4", "out5"};
    prepr5(argv);
  }
  {
    char const* argv[] = {"out5", "tmp2a", "tmp2b"};
    prepr2(argv);
  }
  cat("tmp2a", "tmp2b", "out2");

  sed("<textarea", "<tqqextarea", "out2", "out10");
  sed("</textarea", "</tqqextarea", "out10", "out11");
  sed("<textinput", "<tqqextinput", "out11", "out12");
  sed("</textinput", "</tqqextinput", "out12", "out9");

  {
    char const* argv[] = {"out9", "tmp1a", "tmp1b"};
    prepr1(argv);
  }
  cat("tmp1a", "tmp1b", "out1");
  {
    char const* argv[] = {"out1", ".main_phda9prepr"};
    prepr6(argv);
  }

  return 0;
}

int resto1(char const* argv[]) {
  unsigned char tmp[8192];

  setbuf(stdout, NULL);
  FILE *sf, *t2;
  int c = 0, i, j, k, lnu = 0, f = 0;
  unsigned char *p1 = (unsigned char *)malloc(TMP1A_SIZE + TMP1B_SIZE), *p2 = p1 + TMP1A_SIZE,
                *p0 = p2, *p3 = (unsigned char *)malloc(TMP1A_SIZE + TMP1B_SIZE), *p4 = p3, *p1_base = p1, *p3_base = p3;
  //p1 - pointer to start of the tmp1a 
  //p2 - pointer to start of the tmp1b
  //p3 - pointer to start of the output
  //p4 - ?
  //p0 = p2 - pointer to start of the tmp1b
  sf = fopen(argv[0], "rb");
  fread(p1, TMP1A_SIZE + TMP1B_SIZE, 1, sf);

  do {
    ++lnu; //counter of the number of articles

    j = (int)(strchr((char *)p1, 10) + 1 - (char *)p1);
    if (p1 + j > p0) {
      j = p0 - p1;
    }

    for (i = 0; i < j; i++) {
      if (*(int *)&p1[i] == 'xet<') {
        memcpy(tmp, p1, 40);
        tmp[40] = '\0';
        c = 1, f = lnu; // detected <text> in tmp1a
      }
    }
    for (i = 0; i < j; i++) {
      if (*(int *)&p1[i] == 'et/<') { //detected <\text> in tmp1a
        c = 0;
      }
    }

    if (memcmp(&p1[j - 12], "</revision>", 11) == 0 && c == 1 &&
        (lnu - f >= 4)) {
      do {
        k = (int)(strchr((char *)p2, 10) + 1 - (char *)p2);
        memcpy(p4, p2, k); //memcpy(dest, source, num)
        memcpy(tmp, p2, k);
        tmp[k+1] = '\0';
        p4 += k;
        p2 += k;
      } while (memcmp(&p2[-8], "</text>", 7));
    }

    memcpy(p4, p1, j); //memcpy(dest, source, num)
    p4 += j;
    p1 += j;
  } while (p4 < p3 + TMP1A_SIZE + TMP1B_SIZE);

  t2 = fopen(argv[1], "wb");
  fwrite(p3, TMP1A_SIZE + TMP1B_SIZE, 1, t2);

  fclose(sf);
  fclose(t2);

  free(p3_base);
  free(p1_base);

  return 0;
}


int resto2(char const* argv[])
{
  FILE *sf, *t2;
  int j, k, lastID = 0;
  char *p1=(char *)malloc(TMP2A_SIZE + TMP2B_SIZE), *p2=p1+TMP2A_SIZE, *p0=p2, *p3=(char *)malloc(OUT5_SIZE), *p4=p3, *p1_base = p1, *p3_base = p3;

  sf = fopen(argv[0],"rb");
  fread(p1, TMP2A_SIZE + TMP2B_SIZE, 1, sf);

  do {
    j = (int)(strchr((char *)p1,10) + 1 - (char*)p1);
    if (p1+j > p0)  j = p0-p1;
    memcpy(p4, p1, j);  p4+=j;  p1+=j;

    if (memcmp(&p1[-9],"</title>",8)==0 && *(int*)p1=='    ') {
      lastID += atoi(p2+1);
      sprintf(p4, "    <id>%d</id>%c", lastID, 10);
      p4 = strchr(p4, 10) + 1;
      p2 = strchr(p2, 10) + 1;
      int seenContributor = 0;
      do {
        unsigned int p2bytes =  *(unsigned int*)p2;
        if (p2bytes=='emit') {
          char *p = strchr(p2+=11, ':');
          int d = atoi(p2), hms = atoi(p+1), h = hms/3600;
          sprintf(p4, "      <timestamp>%d-%02d-%02dT%02d:%02d:%02dZ</timestamp>%c", p2[-1] - '0' + 2002, d/31+1, d%31+1, h, hms/60 - h*60, hms%60, 10);
          p4 = strchr(p4, 10) + 1;
          p2 = strchr(p2, 10) + 1;
          continue;
        }
        k = (int)(strchr(p2,10) + 1 - (char*)p2);
        *(int*)p4 = 0x20202020;
        p4 += 4;
        if (p2bytes!='iver' && p2bytes!='tser') {
          *(int*)p4 = 0x20202020;
          p4 += 2;
          if (seenContributor && p2bytes!='noc/')  p4 += 2;
          if (p2bytes=='tnoc')  seenContributor = 1;
        }
        *p4++ = '<';
        memcpy(p4, p2, k);  p4+=k;  p2+=k;
        if (p2[-2]=='>')  {
          continue;
        }
        char *p8 = p2-k;
        assert(p8);
        char *p9 = strchr(p8, '>');
        assert(p9);
        p4[-1]='<';
        *p4   ='/';
        k = p9-p8;
        memcpy(p4+1, p8, k+1);
        p4 += k+3;
        p4[-1] = 10;
      }
      while (memcmp(&p2[-14],"/contributor>",13));
    }
  }
  while (p4 < p3 + OUT5_SIZE);

  t2 = fopen(argv[1],"wb");
  fwrite(p3, OUT5_SIZE, 1, t2);

  fclose(sf);
  fclose(t2);

  free(p3_base);
  free(p1_base);

  return 0;
}



int resto3(char const* argv[]) {
  FILE *sf, *t2;
  int j, k;
  unsigned char *p1 = (unsigned char *)malloc(R3_INPUT_SIZE),
                *p0 = p1,
                *p3 = (unsigned char *)malloc(R3_OUTPUT_ARR_SIZE),
                *p4 = p3, *p1_base = p1, *p3_base = p3;

  sf = fopen(argv[0], "rb");
  fread(p1, R3_INPUT_SIZE, 1, sf);

  do {
    j = *p1++;
    *p4++ = j;

    if (j == '&') {
      k = *p1++;
      if (k == '&') {
        *(int *)p4 = ';pma', p4 += 4;
      } else if (k == '"') {
        *(int *)p4 = 'touq', p4 += 4, *p4++ = ';';
      } else if (k == '<') {
        *(int *)p4 = ' ;tl', p4 += 3;
      } else if (k == '>') {
        *(int *)p4 = ' ;tg', p4 += 3;
      } else {
        *p4++ = k;
      }
    }
  } while (p1 < p0 + R3_INPUT_SIZE);

  t2 = fopen(argv[1], "wb");
  fwrite(p3, p4 - p3, 1, t2);

  fclose(sf);
  fclose(t2);

  free(p3_base);
  free(p1_base);

  return 0;
}



int resto4(char const* argv[]) {
  FILE *sf, *t2;
  int j, k;
  unsigned char *p1 = (unsigned char *)malloc(R4_INPUT_SIZE),
                *p0 = p1,
                *p3 = (unsigned char *)malloc(R4_OUTPUT_ARR_SIZE), *p4 = p3, *p1_base = p1, *p3_base = p3;

  sf = fopen(argv[0], "rb");
  fread(p1, R4_INPUT_SIZE, 1, sf);

  do {
    j = *p1++;
    *p4++ = j;
    if (j == ';' && *(int *)(p1 - 5) == 'pma&') {
      k = *p1++;
      if (k == '"') {
        *(int *)p4 = 'touq', p4 += 4, *p4++ = ';';
      } else if (k == '<') {
        *(int *)p4 = ' ;tl', p4 += 3;
      } else if (k == '>') {
        *(int *)p4 = ' ;tg', p4 += 3;
      } else if (k == '}') {
        *(int *)p4 = 'psbn', p4 += 4, *p4++ = ';';
      } else if (k == '@') {
        *(int *)p4 = 'sadn', p4 += 4, *p4++ = 'h', *p4++ = ';';
      } else if (k == '`') {
        *(int *)p4 = 'sadm', p4 += 4, *p4++ = 'h', *p4++ = ';';
      } else {
        *p4++ = k;
      }
    }
  } while (p1 < p0 + R4_INPUT_SIZE);

  t2 = fopen(argv[1], "wb");
  fwrite(p3, p4 - p3, 1, t2);

  fclose(sf);
  fclose(t2);

  free(p3_base);
  free(p1_base);

  return 0;
}


int resto5(char const* argv[])
{
  char s[16384], z[16384];
  FILE *sf, *t1;
  int mi=0, tf=0;
  sf = fopen(argv[0],"rb");
  t1 = fopen(argv[1],"wb");

  while(1) {
    fgets(s, 16384, sf);
    char *p = strstr(s, "<text "), *w;
    if (p)  { tf=1, p = strchr(p, '>'); assert(p);  if(p[-1]=='/') tf=0;  ++p; }
    else  p = s;

    if (strstr(s, "</text>"))  tf=0;

    if (tf) {
        memcpy(z, s, p-&s[0]);
        for(w=&z[p-&s[0]]; *p!=0; ++p)  {
            if (*p=='"' || *p=='<' || *p=='>')  *w++ = '&';
            *w++ = *p;
        }
        *w = 0;
        fputs(z, t1);
    }
    else fputs(s, t1);

    ++mi;
    if (feof(sf))  break;
  }

  fclose(sf);
  fclose(t1);

  return 0;
}


int phda9_resto() {
  {
    char const* argv[] = {".main_decomp", "out6d"};
    prepr6(argv);
  }
  {
    char const* argv[] = {"out6d", "out1d"};
    resto1(argv);
  }

  sed("</tqqextinput", "</textinput", "out1d", "out13d");
  sed("<tqqextinput","<textinput",  "out13d", "out11d");
  sed("</tqqextarea", "</textarea", "out11d", "out12d");
  sed("<tqqextarea", "<textarea", "out12d", "out10d");

  {
    char const* argv[] = {"out10d", "out2d"};
    resto2(argv);
  }
  {
    char const* argv[] = {"out2d", "out5d"};
    resto5(argv);
  }
  {
    char const* argv[] = {"out5d", "out3d"};
    resto3(argv);
  }
  {
    char const* argv[] = {"out3d", "out4d"};
    resto4(argv);
  }

  sed("&amp;qqq<", "&amp;<", "out4d", "out15d");
  sed("&lt;/tiqqqtle&gt;", "&lt;/title&gt;", "out15d", ".main_decomp_restored");

  return 0;
}
