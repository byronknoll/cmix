// This is used to remap the article numbers in the article order file. By
// skipping redirect article numbers, we can improve compression of the
// article order file.

// To use: ./remap article_path enwik9_path > dictionary/new_article_order

#include <cstdio>
#include <string>
#include <unordered_map>
#include <fstream>
#include <vector>

int main(int argc, char** argv) {
  if (argc != 3) {
    printf("Usage: ./remap article_path enwik9_path\n");
    return 0;
  }
  std::string order_path = argv[1];
  std::string enwik9_path = argv[2];

  std::ifstream infile(enwik9_path);
  std::string line;
  int count1 = -1, count2 = -1;
  bool redirect = false;
  std::vector<std::string> prefix = {
    "      <text xml:space=\"preserve\">#REDIRECT",
    "      <text xml:space=\"preserve\">#redirect",
    "      <text xml:space=\"preserve\">#Redirect",
    "      <text xml:space=\"preserve\">#REdirect",
    "      <text xml:space=\"preserve\">{{softredirect",};
  std::unordered_map<int, int> remap;
  while (getline(infile, line)) {
    for (auto pre : prefix) {
      if (line.rfind(pre, 0) == 0) {
        redirect = true;
        break;
      }
    }
    if (line == "  <page>") {
      if (!redirect) {
        remap[count1] = count2;
        ++count2;
      }
      ++count1;
      redirect = false;
    }
  }
  std::ifstream infile2(order_path);
  while (getline(infile2, line)) {
    int i = stoi(line);
    if (remap.find(i) != remap.end()) {
      printf("%d\n", remap[i]);
    }
  }
}
