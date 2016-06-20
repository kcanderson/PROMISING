/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** coreroutines.cpp
** This class implements some helper functions.
**
** Author: kca
** -------------------------------------------------------------------------*/

#include "coreroutines.h"
#include <algorithm>
#include <sstream>
#include <iostream>

// trim from start
static inline std::string &ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
  return ltrim(rtrim(s));
}

int parseNamesLine(const std::string line, TIndexMap& map) {
  std::string entry;
  int i = 0;
  auto start = 0;
  auto end = line.find("\t");
  while (end != std::string::npos) {
    std::string s = line.substr(start, end - start);
    s = trim(s);
    map[s] = i++;
    //rmap[i] = s;
    start = end + 1;
    end = line.find("\t", start);
  }
  std::string s = line.substr(start, line.back()-start);
  s = trim(s);
  map[s] = i++;
  //  rmap[i] = s;
  
  return i;
}

bool readEntireNetwork(std::ifstream& matFile, const int numNodes, float* const mat) {
  float* curr(mat);
  std::string line;
  int row(0);
  while (!matFile.eof() && row < numNodes) {
    int col(0);
    float v;
    std::getline(matFile, line);
    std::stringstream stream(line);
    while ((stream >> v) && col < numNodes) {
      curr[col] = v;
      col++;
    }
    curr += numNodes;
    row++;
  }
  return true;
}

bool readEntriesFromNetwork(std::ifstream& matFile, TIndexMap& fullMap, const std::vector<std::string>& entries, TIndexMap& newNameMap, float* const mat, const int width) {

  // Get values
  std::vector<int> indices;
  for (auto const &s: entries) {
    if (fullMap.find(s) != fullMap.end()) {
      indices.push_back(fullMap[s]);
    }
  }

  float* curr_mat(mat);
  std::sort(indices.begin(), indices.end());
  std::vector<int>::const_iterator itrow = indices.begin();
  std::map<int, int> newIndexMapping;
  int row = 0, ii = 0, nextRow = *itrow;

  std::string line;
  while (!matFile.eof() && itrow != indices.end()) {
    std::getline(matFile, line);
    if (row == nextRow) {
      std::vector<int>::const_iterator curr_col = indices.begin();
      int nextCol(*curr_col);
      std::stringstream stream(line);
      float v;
      int col = 0, colNew = 0;
      while ((stream >> v) && curr_col != indices.end()) {
	if (col == nextCol) {
	  curr_mat[colNew++] = v;
	  curr_col++;
	  nextCol = *curr_col;
	}
	col++;
      }
      curr_mat += width;
      newIndexMapping[row] = ii++;
      itrow++;
      nextRow = *itrow;
    }
    row++;
  }

  // Add entries for new index mapping
  for (auto const &entry : entries) {
    if (fullMap.find(entry) != fullMap.end()) {
      int v = fullMap[entry];
      if (newIndexMapping.find(v) != newIndexMapping.end()) {
	int vv = newIndexMapping[v];
	newNameMap[entry] = vv;
      }
    }
  }
  
  return true;
}

bool readGroups(std::ifstream& groupStream, TGroups& groups) {
  std::string line;
  while (!groupStream.eof()) {
    std::getline(groupStream, line);
    std::vector<std::string> entries;

    std::stringstream stream(line);
    while (!stream.eof()) {
      std::string entry;
      stream >> entry;
      entry = trim(entry);
      if (!entry.empty()) {
	entries.push_back(entry);
      }
    }

    if (entries.size() > 0) {
      groups.push_back(entries);
    }
  }
  
  return true;
}

void flattenGroups(const TGroups& groups, std::vector<std::string>& entries)
{
  for (auto group : groups) {
    entries.insert(entries.end(), group.begin(), group.end());
  }
}

void mapGroupsToIndices(const TGroups& groups, TIndexMap& map, TIndicesGroups& indicesGroups)
{
  for (auto const &g : groups) {
    std::vector<int> inds;
    for (auto const &e : g) {
      if (map.find(e) != map.end()) {
	inds.push_back(map[e]);
      }
    }
    indicesGroups.push_back(inds);
  }
}

void printMatrix(const float* buffer, const int rows, const int cols) {
  const float* v = buffer;
  for (int y = 0; y < rows; ++y) {
    for (int x = 0; x < cols; ++x) {
      printf("%e\t", v[x]);
    }
    printf("\n");
    v += cols;
  }
}
