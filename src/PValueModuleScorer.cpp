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
#include <math.h>
#include <ctime>

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

void pullOutSubMatrix(const float* const similarities, const int width, std::vector<int>& indices, float* const newMatrix, std::map<int, int>& newIndexMap) {
  // Sort indices from least to greatest
  std::sort(indices.begin(), indices.end());
  const int n = indices.size();
  float* wcurr = newMatrix;
  
  // Pull out submatrix
  for (auto i : indices) {
    const float* curr = similarities + width * i;
    int col = 0;
    for (auto j : indices) {
      wcurr[col++] = curr[j];
    }
    wcurr += n;
  }

  // Map new indices
  int ii = 0;
  for (auto i : indices) {
    newIndexMap[i] = ii++;
  }
}

float calculateZScore(const std::vector<float>& nullScores, const float actualScore) {
  float zscore = (actualScore - mean(nullScores)) / sqrtf(variance(nullScores));
  return zscore;
}

void populateNullDistribution(const int numIterations, std::vector<float>& values) {

}

bool PValueModuleScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, TScoreMap& scores) const
{
  // Seed the random number generator.
  std::srand(std::time(0));
    
  std::map<int, std::vector<float> > allScores;
  int gind = 1;

  int n = 0;
  for (auto const &g : groups) {
    n += g.size();
  }

  // Allocate buffer
  float* const subMatrix = (float*)malloc(sizeof(float) * n * n);
  if (subMatrix == 0) {
    std::cerr << "Could not allocate buffer of " << sizeof(float) * n * n << " bytes" << std::endl;
    exit(-1);
  }
  
  for (auto const &g : groups) {
    std::vector<int> otherIndices;
    for (int i = 0; i < width; ++i) {
      otherIndices.push_back(i);
    }
    for (auto const i : g) {
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
    
    printf("Calculating p-values for locus %d / %d\n", gind++, (int)groups.size());

    // ** Pseudo-randomly choose guys in loci. **
    // [1] Randomly choose starting point
    //int startingPoint = rand() % otherIndices.size();
    for (int i = 0; i < mNumIterations; ++i) {
      //if (i % 1000 == 0) printf("Iteration %d complete.\n", i);
      TIndicesGroups shuffledGroups;
      shuffleGroups(groupSizes, otherIndices, shuffledGroups);
      shuffledGroups.push_back(g);
      // int curr = startingPoint;
      // for (auto const s : groupSizes) {
      // 	std::vector<int> newGroup;
      // 	for (int v = 0; v < s; ++v) {
      // 	  const int me = curr % otherIndices.size();
      // 	  newGroup.push_back(otherIndices[me]);
      // 	  curr++;
      // 	}
      // 	shuffledGroups.push_back(newGroup);
      // }
      // shuffledGroups.push_back(g);
      // startingPoint += (otherIndices.size() - 19);
      // startingPoint %= otherIndices.size();

      // Make temporary matrix-- this is good for cache locality. Much faster.
      std::vector<int> flattenedIndices;
      flattenGroups(shuffledGroups, flattenedIndices);
      std::map<int, int> indexMap;
      pullOutSubMatrix(similarities, width, flattenedIndices, subMatrix, indexMap);

      // Map new indices
      TIndicesGroups newGroups;
      mapGroupsToIndices(shuffledGroups, indexMap, newGroups);
      
      TScoreMap temp;
      //mScorer->ScoreModule(similarities, width, shuffledGroups, temp);
      mScorer->ScoreModule(subMatrix, n, newGroups, temp);
      for (auto const i : g) {
	allScores[i].push_back(temp[indexMap[i]]);
      }
    }
  }

  TScoreMap real;
  mScorer->ScoreModule(similarities, width, groups, real);
  for (auto const &item : allScores) {
    float myScore(real[item.first]);
    // int numBetter(0);
    // for (auto const i : item.second) {
    //   if (i > myScore) {
    // 	numBetter++;
    //   }
    // }
    // scores[item.first] = ((float)numBetter / mNumIterations);
    scores[item.first] = calculateZScore(item.second, myScore);
  }

  free(subMatrix);
  return true;
}

bool pvalCompare(const std::pair<int, float>& firstElem, const std::pair<int, float>& secondElem) {
  return firstElem.second < secondElem.second;
}

bool zscoreCompare(const std::pair<int, float>& firstElem, const std::pair<int, float>& secondElem) {
  return firstElem.second > secondElem.second;
}

void calculatePValues(const TScoreMap& scores, TScoreMap& pvalScores) {
  for (auto const &e : scores) {
    pvalScores[e.first] = 1.0f - phi(e.second);
  }
}


void benjaminiHochbergCorrection(const TScoreMap& pvals, TScoreMap& pvalAdjusted) {
  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(pvals, pairs, pvalCompare);
  
  int n(pairs.size());
  int k = 1;
  for (auto const &e : pairs) {
    pvalAdjusted[e.first] = std::min(1.0f, e.second * (float)n / (float)k);
    k++;
  }
}

void bonferroniCorrection(const TScoreMap& pvals, TScoreMap& pvalAdjusted) {
  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(pvals, pairs, pvalCompare);
  
  int n(pairs.size());
  for (auto const &e : pairs) {
    pvalAdjusted[e.first] = std::min(1.0f, e.second * (float)n);
  }
}

void PValueModuleScorer::BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const
{
  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(scores, pairs, zscoreCompare);
  //sortMapByVal(scores, pairs, pvalCompare);
  TScoreMap pvals, pvalsAdjusted;
  //TScoreMap pvalsAdjusted;
  calculatePValues(scores, pvals);
  //benjaminiHochbergCorrection(pvals, pvalsAdjusted);
  bonferroniCorrection(pvals, pvalsAdjusted);
  //bonferroniCorrection(scores, pvalsAdjusted);
  
  out << std::endl << "Top 10 genes" << std::endl << "Gene\tZ score\tP value\tAdj. p value" << std::endl << "----------------" << std::endl;
  int i = 0;
  for (auto it = pairs.begin(); it != pairs.end() && i < 10; it++, i++) {
    out << rmap[it->first] << "\t" << it->second << "\t" << pvals[it->first] << "\t" << pvalsAdjusted[it->first] << std::endl;
    //out << rmap[it->first] << "\t" << it->second << "\t" << scores[it->first] << "\t" << pvalsAdjusted[it->first] << std::endl;
  }
}

void PValueModuleScorer::LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const
{
  //TScoreMap pvalsAdjusted;
  TScoreMap pvals, pvalsAdjusted;
  calculatePValues(scores, pvals);
  //benjaminiHochbergCorrection(pvals, pvalsAdjusted);
  bonferroniCorrection(pvals, pvalsAdjusted);
  //bonferroniCorrection(scores, pvalsAdjusted);
 
  int l = 1;
  for (auto const &g : groups) {
    out << "----------------" << std::endl;
    out << "Locus " << l << std::endl;
    out << "Gene\tZ score\tP value\tAdj. p value" << std::endl;
    out << "----------------" << std::endl;
    TScoreMap locusScores;
    for (auto const &e : g) {
      locusScores[e] = scores[e];
    }
    std::vector<std::pair<int, float> > pairs;
    //sortMapByVal(locusScores, pairs, zscoreCompare);
    sortMapByVal(locusScores, pairs, pvalCompare);
    for (auto const &e : pairs) {
      out << rmap[e.first] << "\t" << e.second << "\t" << pvals[e.first] << "\t" << pvalsAdjusted[e.first] << std::endl;
      //out << rmap[e.first] << "\t" << e.second << "\t" << scores[e.first] << "\t" << pvalsAdjusted[e.first] << std::endl;
    }
    l++;
  }

}
