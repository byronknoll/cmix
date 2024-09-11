#ifndef ARTICLE_REORDER_H 
#define ARTICLE_REORDER_H 


#include <string.h>
#include <string>
#include <vector>
#include <cstdio>
#include <unordered_map>
#include <fstream>

#define NUM_OF_ARTICLES 243425

struct Accumulator {
  int id;
  int start;
  int end;
   #ifdef DUMPARTICLE
  std::string title;
  std::string infobox;
  std::string redirect;
  #endif
};

enum ParserState {
  expect_page = 0,
  expect_id,
  expect_pageend
};

int line_count = 0;
static char s[8192*8];

void bubblesort(std::vector<Accumulator>& mylist) {
    for (int i = 1; i < mylist.size(); i++)	{
        for (int j = 0; j < mylist.size() - i; j++) {
            if (mylist[j].id > mylist[j + 1].id) {
            std::swap(mylist[j], mylist[j + 1]);
            }
       }
   }
}
const char *patterns1[] = { "<page>", "<id>", "</page>" };
const int transitions1[3] = { expect_id, expect_pageend, expect_page };
std::vector<Accumulator> vec;
std::vector<std::string> lines;
// in phda9
int wfgets(char *str, int count, FILE  *fp);
void wfputs(const char *str,FILE *fp);

void loadFile(const char *fname) {
   
    line_count = 0;
    int  state = expect_page;
    Accumulator acc;
    FILE* file = fopen(fname, "rb");
    while (wfgets(s, 8192*8, file) )  {
         #ifdef DUMPARTICLE
        char *pt = strstr(s, "<title>");
        if (pt) acc.title=s,acc.title.pop_back();
        if (strstr(s, "infobox ") || strstr(s, "infobox_")  || strstr(s, "Infobox_") || strstr(s, "Infobox ") ) acc.infobox=s,acc.infobox.pop_back();
        if (strstr(s, "#REDIRECT") || strstr(s, "#redirect") || strstr(s, "softredirect")|| strstr(s, "#Redirect")|| strstr(s, "#REdirect")) acc.redirect=s,acc.redirect.pop_back();
        #endif
    char *p = strstr(s, patterns1[state]);
    if (p) {
      if (state==expect_page){
           acc.start = line_count;
      }else if(state==expect_id){
           char *p = strstr(s, patterns1[1]); //id
           if (p){
               p=p+4;
               acc.id = atoi(p);
           }
      }else if(state==expect_pageend){
           acc.end = line_count;
      }
      state = transitions1[state];
      if (state == expect_page){
        vec.push_back(acc);
         #ifdef DUMPARTICLE
        acc.title="",acc.infobox="",acc.redirect="";
        #endif
       }
    } 
    line_count++;
    std::string so=std::string(s);
    lines.push_back(so);
  }
  fclose(file);
  
}
void reorder() {
    std::vector<int> positions;
    std::vector<int> used(NUM_OF_ARTICLES, 0);
    
    loadFile(".main");

    // Article numbers need to be remapped (see article_remap.cpp).
    // Here we read enwik9 in order to generate the article number mapping.
    std::unordered_map<int, int> remap;
    std::ifstream infile(".main");
    std::string line;
    int count1 = -1, count2 = -1;
    bool redirect = false;
    std::vector<std::string> prefix = {
      "      <text xml:space=\"preserve\">#REDIRECT",
      "      <text xml:space=\"preserve\">#redirect",
      "      <text xml:space=\"preserve\">#Redirect",
      "      <text xml:space=\"preserve\">#REdirect",
      "      <text xml:space=\"preserve\">{{softredirect",};
    while (getline(infile, line)) {
      for (auto pre : prefix) {
        if (line.rfind(pre, 0) == 0) {
          redirect = true;
          break;
        }
      }
      if (line == "  <page>") {
        remap[count2] = count1;
        if (!redirect) {
          ++count2;
        }
        ++count1;
        redirect = false;
      }
    }

    std::ifstream infile2("new_article_order");
    while (getline(infile2, line)) {
      int res = remap[stoi(line)];  // Article numbers get remapped.
      positions.push_back(res);
      used[res] = 1;
    }
  
    if (positions.size() < NUM_OF_ARTICLES) {
        for (int i = 0; i < NUM_OF_ARTICLES; i++) {
            if (used[i] == 0) {
                positions.push_back(i);
            }
       }
    }
			  
  FILE* out = fopen(".main_reordered", "wb");
  std::string so;
  for(int i = 0; i < positions.size(); i++) {
    int pos = positions[i];
     #ifdef DUMPARTICLE
     printf("%d\t%d\t%d\t%d\t%s\t%s\t%s\n",pos,vec[pos].id,vec[pos].start,vec[pos].end,vec[pos].title.c_str(),vec[pos].infobox.c_str(),vec[pos].redirect.c_str());
     #else
    for(int j = vec[pos].start; j <= vec[pos].end; j++) {
      so=lines[j];
      wfputs(so.c_str(),out);
    }
    #endif
  } 
  #ifdef DUMPARTICLE
  exit(0);
  #endif
  fclose(out);
  vec.clear();
  std::vector<Accumulator>(vec).swap(vec);  // Why C++? why
  lines.clear();
  std::vector<std::string>(lines).swap(lines);
}

void sort() {
  loadFile(".main_decomp_restored");
 
  bubblesort(vec);
	  
  FILE* out = fopen(".main_decomp_restored_sorted", "wb");
  std::string so;
  for(int i =0; i < vec.size(); i++) {
    for(int j = vec[i].start; j <= vec[i].end; j++) {
      so=lines[j];
      wfputs(so.c_str(),out);
    }
  }
  
  fclose(out);
  vec.clear();
  std::vector<Accumulator>(vec).swap(vec);
  lines.clear(); 
  std::vector<std::string>(lines).swap(lines); 
}
#endif // ARTICLE_REORDER_H 
