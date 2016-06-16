/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** coreroutines.h
** This file declares a few helper functions.
**
** Author: kca
** -------------------------------------------------------------------------*/

#ifndef COREROUTINES_H
#define COREROUTINES_H

#include "../include/IModuleScorer.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>

typedef std::map<std::string, int> TIndexMap;
typedef std::vector< std::vector<std::string> > TGroups;

int parseNamesLine(const std::string line, TIndexMap& map);
bool readEntriesFromNetwork(std::ifstream& matFile, TIndexMap& fullMap, const std::vector<std::string>& entries, TIndexMap& newNameMap, float* const mat, const int width);
bool readEntireNetwork(std::ifstream& matFile, const int numNodes, float* const mat);
bool readGroups(std::ifstream& groupFile, TGroups& groups);
void flattenGroups(const TGroups& groups, std::vector<std::string>& entries);
void mapGroupsToIndices(const TGroups& groups, TIndexMap& map, TIndicesGroups& indicesGroups);
void printMatrix(const float* buffer, const int rows, const int cols);

template <typename T> void printGroups(const std::vector< std::vector<T> >& groups) {
  for (auto const &g : groups) {
    for (auto item : g) {
      std::cout << item << " ";
    }
    std::cout << std::endl;
  }
}

template <typename C>
void sortMapByVal(const std::map<int, float>& map, std::vector<std::pair<int, float> >& sorted, C& comp) {
  for (auto it = map.begin(); it != map.end(); ++it) {
    sorted.push_back(*it);
  }
  std::sort(sorted.begin(), sorted.end(), comp);
}

#endif
