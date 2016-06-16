/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** PValueModuleScorer.cpp
** This class scores a gene based on an empirical p-value.
**
** Author: kca
** -------------------------------------------------------------------------*/

#include "PValueModuleScorer.h"
#include "coreroutines.h"
#include <algorithm>

void shuffleGroups(const std::vector<int> groupSizes, std::vector<int>& allIndices, TIndicesGroups& shuffledGroups) {
  std::random_shuffle(allIndices.begin(), allIndices.end());
  std::vector<int>::const_iterator it(allIndices.begin());
  
  for (auto size : groupSizes) {
    std::vector<int> newGroup;
    for (int i = 0; i < size && it != allIndices.end(); ++i, it++) {
      int item(*it);
      newGroup.push_back(item);
    }
    shuffledGroups.push_back(newGroup);
  }
}

bool PValueModuleScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, TScoreMap& scores) const
{
  std::map<int, std::vector<float> > allScores;
  
  for (auto const &g : groups) {
    std::vector<int> otherIndices;
    for (int i = 0; i < width; ++i) {
      otherIndices.push_back(i);
    }
    for (auto i : g) {
      auto it = std::find(otherIndices.begin(), otherIndices.end(), i);
      otherIndices.erase(it);
    }
    TIndicesGroups otherGroups(groups);
    auto it = std::find(otherGroups.begin(), otherGroups.end(), g);
    otherGroups.erase(it);
    std::vector<int> groupSizes;
    for (auto const &o : otherGroups) {
      groupSizes.push_back(o.size());
    }
    for (int i = 0; i < mNumIterations; ++i) {
      //if (i % 1000 == 0) printf("Iteration %d complete.\n", i);
      TIndicesGroups shuffledGroups;
      shuffleGroups(groupSizes, otherIndices, shuffledGroups);
      shuffledGroups.push_back(g);
      TScoreMap temp;
      mScorer->ScoreModule(similarities, width, shuffledGroups, temp);
      for (auto i : g) {
	allScores[i].push_back(temp[i]);
      }
    }
  }

  TScoreMap real;
  mScorer->ScoreModule(similarities, width, groups, real);
  for (auto const &item : allScores) {
    float myScore(real[item.first]);
    int numBetter(0);
    
    for (auto i : item.second) {
      if (i > myScore) {
	numBetter++;
      }
    }
    scores[item.first] = ((float)numBetter / mNumIterations);
  }

  return true;
}

bool pvalCompare(const std::pair<int, float>& firstElem, const std::pair<int, float>& secondElem) {
  return firstElem.second < secondElem.second;
}

void PValueModuleScorer::BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const
{
  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(scores, pairs, pvalCompare);
  
  out << std::endl << "Top 10 genes" << std::endl << "Gene\tP-value" << std::endl << "----------------" << std::endl;
  int i = 0;
  for (auto it = pairs.begin(); it != pairs.end() && i < 10; it++, i++) {
    out << rmap[it->first] << "\t" << it->second << std::endl;
  }
}

void PValueModuleScorer::LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const
{
  int l = 1;
  for (auto const &g : groups) {
    out << "----------------" << std::endl;
    out << "Locus " << l << std::endl;
    out << "Gene\tP-value" << std::endl;
    out << "----------------" << std::endl;
    TScoreMap locusScores;
    for (auto const &e : g) {
      locusScores[e] = scores[e];
    }
    std::vector<std::pair<int, float> > pairs;
    sortMapByVal(locusScores, pairs, pvalCompare);
    for (auto const &e : pairs) {
      out << rmap[e.first] << "\t" << e.second << std::endl;
    }
    l++;
  }

}
