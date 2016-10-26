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
#include <set>

// void PValueModuleScorer::ShuffleGroups(const std::vector<int> groupSizes, std::vector<int>& allIndices, TIndicesGroups& shuffledGroups) const {
//   std::random_shuffle(allIndices.begin(), allIndices.end());
//   std::vector<int>::const_iterator it(allIndices.begin());
  
//   for (auto size : groupSizes) {
//     std::vector<int> newGroup;
//     for (int i = 0; i < size && it != allIndices.end(); ++i, it++) {
//       int item(*it);
//       newGroup.push_back(item);
//     }
//     shuffledGroups.push_back(newGroup);
//   }
// }

void PValueModuleScorer::ShuffleGroups(const TIndicesGroups& groups, const std::vector<int>& excisedIndices, TIndicesGroups& shuffledGroups) const {
  std::set<int> indices(excisedIndices.begin(), excisedIndices.end());
  for (auto const& g : groups) {
    std::vector<int> newGroup;
    for (auto i : g) {
      auto const& it = mNodeDegreeMap.find(i);
      const int degreeGroup = mNodeDegreeMap.find(i)->second;

      const std::vector<int>* group = &(mDegreeGroupNodes.find(degreeGroup)->second);
      int j = rand() % group->size();
      int v = (*group)[j];
      while (indices.find(v) != indices.end()) {
      	j = rand() % group->size();
      	v = (*group)[j];
      }

      indices.insert(v);
      newGroup.push_back(v);
    }
    shuffledGroups.push_back(newGroup);
  }

}

// void PValueModuleScorer::ShuffleGroups(const TIndicesGroups& groups, const std::vector<int>& excisedIndices, TIndicesGroups& shuffledGroups) const {
//   std::set<int> indices(excisedIndices.begin(), excisedIndices.end());
//   for (auto const& g : groups) {
//     std::vector<int> newGroup;
//     for (auto i : g) {
//       newGroup.push_back(rand() % 19100);
//     }
//     shuffledGroups.push_back(newGroup);
//   }

// }



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

bool PValueModuleScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, const TIndices& indicesToScore, TScoreMap& scores) const
{
  // Seed the random number generator.
  std::srand(std::time(0));
  
  std::map<int, std::vector<float> > allScores;
  int gind = 1;

  int n = 0;
  for (auto const &g : groups) {
    n += g.size();
  }

  //std::map<int, float> degrees;
  //degreeMatrix(similarities, width, degrees);

  // Allocate buffer
#define TEMP_BUFFER 0
#if TEMP_BUFFER
  float* const subMatrix = (float*)malloc(sizeof(float) * n * n);
  if (subMatrix == 0) {
    std::cerr << "Could not allocate buffer of " << sizeof(float) * n * n << " bytes" << std::endl;
    exit(-1);
  }
#endif
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
    //std::vector<int> groupSizes;
    //for (auto const &o : otherGroups) {
    //groupSizes.push_back(o.size());
    //}
    printf("Calculating p-values for locus %d / %d\n", gind++, (int)groups.size());
    // ** Pseudo-randomly choose guys in loci. **
    // [1] Randomly choose starting point
    //int startingPoint = rand() % otherIndices.size();
    
    for (int i = 0; i < mNumIterations; ++i) {
      //if (i % 1000 == 0) printf("Iteration %d complete.\n", i);
      TIndicesGroups shuffledGroups;
      ShuffleGroups(otherGroups, g, shuffledGroups);
      shuffledGroups.push_back(g);
      
      TScoreMap temp;
#if TEMP_BUFFER
      // Make temporary matrix-- this is good for cache locality. Much faster.
      std::vector<int> flattenedIndices;
      flattenGroups(shuffledGroups, flattenedIndices);
      std::map<int, int> indexMap;
      pullOutSubMatrix(similarities, width, flattenedIndices, subMatrix, indexMap);

      // Map new indices
      TIndicesGroups newGroups;
      mapGroupsToIndices(shuffledGroups, indexMap, newGroups);
      TIndices newInds;
      flattenGroups(newGroups, newInds);
      mScorer->ScoreModule(subMatrix, n, newGroups, newInds, temp);
      
#else
      mScorer->ScoreModule(similarities, width, shuffledGroups, g, temp);
#endif
      for (auto const i : g) {
#if TEMP_BUFFER
      	allScores[i].push_back(temp[indexMap[i]]);
#else
	allScores[i].push_back(temp[i]);
#endif
      }
    }
  }

  TScoreMap real;
  mScorer->ScoreModule(similarities, width, groups, indicesToScore, real);
  for (auto const &item : allScores) {
    float myScore(real[item.first]);
#define ZSCORE 0
#if ZSCORE
    scores[item.first] = calculateZScore(item.second, myScore);
#else
    int numBetter(0);
    for (auto const i : item.second) {
      if (i > myScore) {
    	numBetter++;
      }
    }
    scores[item.first] = ((float)numBetter / mNumIterations);
#endif
  }

#if TEMP_BUFFER
  free(subMatrix);
#endif
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
#if ZSCORE
  TScoreMap pvals, pvalsAdjusted;
  sortMapByVal(scores, pairs, zscoreCompare);
  calculatePValues(scores, pvals);
  bonferroniCorrection(pvals, pvalsAdjusted);
#else
  sortMapByVal(scores, pairs, pvalCompare);
  TScoreMap pvalsAdjusted;
  bonferroniCorrection(scores, pvalsAdjusted);
#endif
  //benjaminiHochbergCorrection(pvals, pvalsAdjusted);

#if ZSCORE
  out << std::endl << "Top 10 genes" << std::endl << "Gene\tZ score\tP value\tAdj. p value" << std::endl << "----------------" << std::endl;
#else
  out << std::endl << "Top 10 genes" << std::endl << "Gene\tP value\tAdj. p value" << std::endl << "----------------" << std::endl;
#endif
  int i = 0;
  for (auto it = pairs.begin(); it != pairs.end() && i < 10; it++, i++) {
#if ZSCORE
    out << rmap[it->first] << "\t" << it->second << "\t" << pvals[it->first] << "\t" << pvalsAdjusted[it->first] << std::endl;
#else
    out << rmap[it->first] << "\t" << it->second << "\t" << pvalsAdjusted[it->first] << std::endl;
#endif
  }
}

void PValueModuleScorer::LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const
{
#if ZSCORE
  TScoreMap pvals, pvalsAdjusted;
  calculatePValues(scores, pvals);
  bonferroniCorrection(pvals, pvalsAdjusted);
#else
  TScoreMap pvals(scores);
  TScoreMap pvalsAdjusted;
  bonferroniCorrection(scores, pvalsAdjusted);
#endif
  
  //benjaminiHochbergCorrection(pvals, pvalsAdjusted);
  
#define OLD_FORMAT 0
#if OLD_FORMAT
  int l = 1;
  for (auto const &g : groups) {
    out << "----------------" << std::endl;
    out << "Locus " << l << std::endl;
#if ZSCORE
    out << "Gene\tZ score\tP value\tAdj. p value" << std::endl;
#else
    out << "Gene\tP value\tAdj. p value" << std::endl;
#endif
    out << "----------------" << std::endl;
    TScoreMap locusScores;
    for (auto const &e : g) {
      locusScores[e] = scores[e];
    }
    std::vector<std::pair<int, float> > pairs;
#if ZSCORE
    sortMapByVal(locusScores, pairs, zscoreCompare);
#else
    sortMapByVal(locusScores, pairs, pvalCompare);
#endif
    for (auto const &e : pairs) {
#if ZSCORE
      out << rmap[e.first] << "\t" << e.second << "\t" << pvals[e.first] << "\t" << pvalsAdjusted[e.first] << std::endl;
#else
      out << rmap[e.first] << "\t" << e.second << "\t" << pvalsAdjusted[e.first] << std::endl;
#endif
    }
    l++;
  }
#else
  // New formatting
  // Gene Locus P-val Adj. p-val
  out << "Gene\tLocus\tP-val\tAdj. p-val" << std::endl;
  std::map<int, std::string> groupMap;
  int i = 0;
  for (auto const &g : groups) {
    std::string locus = "Locus " + std::to_string(i++);
    for (auto gene : g) {
      groupMap[gene] = locus;
    }
  }

  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(pvals, pairs, pvalCompare);
  for (auto const &e : pairs) {
    out << rmap[e.first] << "\t" << groupMap[e.first] << "\t" << pvals[e.first] << "\t" << pvalsAdjusted[e.first] << std::endl;
  }
    
  
#endif
}
