#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <fstream>
#include <filesystem>
#include <chrono>
// Smith-Waterman

#define inf std::numeric_limits<int>::min();

typedef std::vector<int> vi;
typedef std::vector<vi> vvi;
const std::string RESET_COLOR = "\033[0m";
const std::string RED_TEXT = "\033[31m";
const std::string GREEN_TEXT = "\033[32m";
const std::string YELLOW_TEXT = "\033[33m";
const std::string BLUE_TEXT = "\033[34m";
const std::string MAGENTA_TEXT = "\033[35m";
const std::string CYAN_TEXT = "\033[36m";
const std::string BOLD_TEXT = "\033[1m";
const int penaltyScore = 2;
typedef long long ll;

using namespace std;

inline int max4(int a, int b, int c, int d)
{
  return max(max(max(a, b), c), d);
}

int maxI, maxJ, countSeq;

void printAlignment(vector<int> &path, string &s, string &t)
{
  cout << "[" << ++countSeq << "]\n";
  int pathSize = path.size();
  int si = maxI - pathSize;
  int ti = maxJ - pathSize;

  for (int dir = pathSize - 1; dir >= 0; dir--)
  {
    if (path[dir] == 2)
      si++;
    if (path[dir] == 3)
      ti++;
  }

  for (int dir = pathSize - 1; dir >= 0; dir--)
  {
    if (path[dir] == 2)
      cout << "_";
    else
      cout << s[si++];
  }
  cout << '\n';

  for (int dir = pathSize - 1; dir >= 0; dir--)
  {
    if (path[dir] == 2 || path[dir] == 3)
      cout << " ";
    if (path[dir] == 0)
      cout << "*";
    if (path[dir] == 1)
      cout << "|";
  }
  cout << '\n';
  for (int dir = pathSize - 1; dir >= 0; dir--)
  {
    if (path[dir] == 3)
      cout << "_";
    else
      cout << t[ti++];
  }
  cout << '\n';
}

void traceback(int i, int j, vvi &dp, string &s, string &t, vector<int> &path)
{
  if (i == 0 && j == 0)
  {
    printAlignment(path, s, t);
    return;
  }
  else
  {
    int diag, up, left, max;
    bool match = s[i - 1] == t[j - 1];

    diag = up = left = inf;

    if (i > 0)
    {
      if (j > 0)
      {
        diag = dp[i - 1][j - 1] + (match ? +1 : -1);
        left = dp[i - 1][j] - penaltyScore;
        up = dp[i][j - 1] - penaltyScore;
      }
      else
      {
        left = dp[i - 1][j] - penaltyScore;
      }
    }
    else
    {
      up = dp[i][j - 1] - penaltyScore;
    }

    max = max4(diag, up, left, 0);

    if (max != 0)
    {
      if (diag == max)
      {
        path.push_back(!match);
        traceback(i - 1, j - 1, dp, s, t, path);
        path.pop_back();
      }
      if (up == max)
      {
        path.push_back(2);
        traceback(i, j - 1, dp, s, t, path);
        path.pop_back();
      }
      if (left == max)
      {
        path.push_back(3);
        traceback(i - 1, j, dp, s, t, path);
        path.pop_back();
      }
    }
    else
    {
      traceback(0, 0, dp, s, t, path);
    }
  }
}

void show(vvi &dp)
{
  for (auto v : dp)
  {
    for (auto i : v)
      std::cout << i << ' ';
    std::cout << '\n';
  }
}

void sw(std::string &s, std::string &t)
{
  int rows, cols;
  rows = s.size();
  cols = t.size();

  vector<pair<int, int>> maxIndexes;
  std::vector<std::vector<int>> dp;
  dp.resize(rows + 1, std::vector<int>(cols + 1, 0));

  auto start = chrono::high_resolution_clock::now();

  int maxScore = inf;

  for (int i = 1; i <= rows; i++)
    for (int j = 1; j <= cols; j++)
    {
      dp[i][j] = max4(
          dp[i - 1][j - 1] + (s[i - 1] == t[j - 1] ? +1 : -1),
          dp[i - 1][j] - penaltyScore,
          dp[i][j - 1] - penaltyScore,
          0);
      if (maxScore < dp[i][j])
      {
        maxScore = dp[i][j];
        maxIndexes.clear();
        maxIndexes.push_back(make_pair(i, j));
      }
      else if (maxScore == dp[i][j])
      {
        maxIndexes.push_back(make_pair(i, j));
      }
    }

  auto stop = chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> dpTime = stop - start;

  cout << BOLD_TEXT << GREEN_TEXT << "[ INFO | Original Sequences ]\n"
       << RESET_COLOR;
  cout << "s: " << s << '\n';
  cout << "t: " << t << "\n\n";

  cout << BOLD_TEXT << GREEN_TEXT << "[ INFO | Possible Alignments ]\n"
       << RESET_COLOR;
  start = chrono::high_resolution_clock::now();

  for (auto maxIndex : maxIndexes)
  {
    vector<int> tmpPath;
    maxI = maxIndex.first;
    maxJ = maxIndex.second;
    traceback(maxIndex.first, maxIndex.second, dp, s, t, tmpPath);
  }

  stop = chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> seqTime = stop - start;

  cout << BOLD_TEXT << GREEN_TEXT << "[ INFO | Summary]\n";
  cout << "max score: " << maxScore << '\n';
  cout << "n of alignments: " << countSeq << '\n';
  cout << "score matrix: " << dpTime.count() << " milliseconds\n";
  cout << "traceback: " << seqTime.count() << " milliseconds\n";
  cout << "total time: " << dpTime.count() + seqTime.count() << " milliseconds\n";
}

vector<string> getFilenames(string &dir)
{
  vector<string> filenames;
  for (const auto &entry : filesystem::directory_iterator(dir))
  {
    if (entry.is_regular_file())
      filenames.push_back(entry.path());
  }
  return filenames;
}

unordered_map<string, string> getSequences(vector<string> &filenames)
{
  unordered_map<string, string> sequences;

  for (auto filename : filenames)
  {
    fstream file(filename);
    string word, sequence, name;

    file >> name;

    while (file >> word)
    {
      if (word.size() == 10)
        sequence += word;
    }
    sequences[name] = sequence;
  }
  return sequences;
}

void printSequences(unordered_map<string, string> &sequences)
{
  for (auto seq : sequences)
  {
    cout << "Name: " << seq.first << '\n';
    cout << "Length: " << seq.second.size() << '\n';
    cout << "Sequence: " << seq.second << "\n\n";
  }
}

int main()
{
  countSeq = 0;
  string dir = "../sequences";
  vector<string> filenames = getFilenames(dir);
  unordered_map<string, string> sequences = getSequences(filenames);
  printSequences(sequences);

  std::string s, t;
  s = sequences["Bacteria"];
  t = sequences["Sars-Cov"];
  // s = "AATCGWWWGWGW";
  // t = "AACGGWWGWWGW";
  sw(s, t);
}