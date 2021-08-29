
#define NUM_OF_ARTICLES 243425

struct Accumulator {
  int id;
  int start;
  int end;
};

enum ParserState {
  expect_page = 0,
  expect_id,
  expect_pageend
};

int line_count = 0;

int action_get_line_count(std::string s) {
  return line_count;
}

int action_get_id(std::string s) {
  std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
  s.erase(end_pos, s.end());

  std::string tok = "<id>";
  std::string::size_type i = s.find(tok);
  if (i != std::string::npos)
     s.erase(i, tok.length());

  tok = "</id>";
  i = s.find(tok);
  if (i != std::string::npos)
     s.erase(i, tok.length());

  return std::stoi(s);
}

void save_pagestart(int res, Accumulator* acc) {
  acc->start = res;
}

void save_id(int res, Accumulator* acc) {
  acc->id = res;
}

void save_pageend(int res, Accumulator* acc) {
  acc->end = res;
}

void bubblesort(std::vector<Accumulator>& mylist)
{
	for (int i = 1; i < mylist.size(); i++)
	{

	    for (int j = 0; j < mylist.size() - i; j++) 
	    {
		if (mylist[j].id > mylist[j + 1].id) 
		{
			std::swap(mylist[j], mylist[j + 1]);
		}
	    }
	}
}

int reorder() {
  line_count = 0;

  std::ifstream file(".main"); //file just has some sentences

  std::ifstream order_file(".new_article_order"); //file just has some sentences

  std::vector<std::string> lines;
  std::vector<int> positions;

  std::vector<std::string> patterns  = { "<page>", "<id>", "</page>" };
  std::vector<ParserState>  transitions    = { expect_id, expect_pageend, expect_page };
  int (*actions[3])(std::string)     = {action_get_line_count, action_get_id, action_get_line_count};
  void (*save[3])(int, Accumulator*) = {save_pagestart, save_id, save_pageend};

  ParserState state = expect_page;

  std::vector<Accumulator> vec;

  std::string s;
  std::string pattern;
  int res = 0;
  Accumulator acc;
  while (std::getline(file, s))
  {
    pattern = patterns[state];
    if (s.find(pattern) != std::string::npos) {
      res = actions[state](s);
      save[state](res, &acc);
      state = transitions[state];
      if (state == expect_page)
        vec.push_back(acc);
    } 
    line_count++;
    lines.push_back(s);
  }

  std::vector<int> used(NUM_OF_ARTICLES, 0);
  while (std::getline(order_file, s)) {
    positions.push_back(std::stoi(s));
    used[std::stoi(s)] = 1;
  }

  if (positions.size() < NUM_OF_ARTICLES) {
	for (int i = 0; i < NUM_OF_ARTICLES; i++) {
		if (used[i] == 0) {
			positions.push_back(i);
		}
	}
  }
			  

  std::ofstream out(".main_reordered");
  for(int i = 0; i < positions.size(); i++) {
    int pos = positions[i];
    for(int j = vec[pos].start; j <= vec[pos].end; j++) {
      out << lines[j] << std::endl;
    }
  } 
  out.close();

  return 0;
}


int sort() {
  line_count = 0;

  std::ifstream file(".main_decomp_restored"); 
  if (!file) {
    return -1;
  }

  std::vector<std::string> lines;

  std::vector<std::string> patterns  = { "<page>", "<id>", "</page>" };
  std::vector<ParserState>  transitions    = { expect_id, expect_pageend, expect_page };
  int (*actions[3])(std::string)     = {action_get_line_count, action_get_id, action_get_line_count};
  void (*save[3])(int, Accumulator*) = {save_pagestart, save_id, save_pageend};

  ParserState state = expect_page;

  std::vector<Accumulator> vec;

  std::string s;
  std::string pattern;
  int res = 0;
  Accumulator acc;
  while (std::getline(file, s))
  {
    pattern = patterns[state];
    if (s.find(pattern) != std::string::npos) {
      res = actions[state](s);
      save[state](res, &acc);
      state = transitions[state];
      if (state == expect_page)
        vec.push_back(acc);
    } 
    line_count++;
    lines.push_back(s);
  }

  bubblesort(vec);

  std::ofstream out(".main_decomp_restored_sorted");
  if (!out) {
    return -1;
  }
  for(int i =0; i < vec.size(); i++) {
    for(int j = vec[i].start; j <= vec[i].end; j++) {
      out << lines[j] << std::endl;
    }
  } 
  out.close();

  return 0;
}

