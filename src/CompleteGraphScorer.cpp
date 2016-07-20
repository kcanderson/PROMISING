/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** CompleteGraphScorer.cpp
** This file implements the CompleteGraphScorer class.
**
** Author: kca
** -------------------------------------------------------------------------*/

#include "CompleteGraphScorer.h"
#include "coreroutines.h"
#include "string.h"
#include <algorithm>

typedef std::vector<int> TInts;

bool CompleteGraphScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, TScoreMap& scores) const
{
  // Count total number of nodes in all groups
  // int numInGroups = 0;
  // for (auto const &g : groups) {
  //   numInGroups += g.size();
  // }

  // Allocate memory for scores
  float* const rawScores = (float*)malloc(sizeof(float) * width);
  if (rawScores == 0) {
    printf("Could not allocate buffer of %d bytes\n", (int)sizeof(float) * width);
    return false;
  }

  float* const intermediateScores = (float*)malloc(sizeof(float) * width);
  if (intermediateScores == 0) {
    printf("Could not allocate buffer of %d bytes\n", (int)sizeof(float) * width);
    return false;
  }

  for (int i = 0; i < width; ++i) {
    rawScores[i] = 0.0f;
    intermediateScores[i] = 0.0f;
  }
  
  const int gsize = 3;
  const int numGroups = groups.size();
  
  for (int i = 0; i <= numGroups - gsize; ++i) {
    const TInts g1 = groups[i];
    for (int j = i + 1; j <= numGroups - gsize + 1; ++j) {
      const TInts g2 = groups[j];
      for (int k = j + 1; k <= numGroups - gsize + 2; ++k) {
	const TInts g3 = groups[k];

	for (const auto i1 : g1) {
	  const float* r1 = similarities + width * i1;
	  for (const auto i2 : g2) {
	    const float* r2 = similarities + width * i2;
	    const float v12 = *(r1 + i2);
	    for (const auto i3 : g3) {
	      // i1 - i2
	      float acc = v12;
	      // i1 - i3
	      acc += *(r1 + i3);
	      // i2 - i3
	      acc += *(r2 + i3);

	      if (acc > intermediateScores[i1]) intermediateScores[i1] = acc;
	      if (acc > intermediateScores[i2]) intermediateScores[i2] = acc;
	      if (acc > intermediateScores[i3]) intermediateScores[i3] = acc;
	    }
	  }
	}
	for (int i = 0; i < width; ++i) {
	  rawScores[i] += intermediateScores[i];
	  intermediateScores[i] = 0.0f;
	}
      }
    }
  }

  for (auto const &g : groups) {
    for (auto i : g) {
      scores[i] = rawScores[i];
    }
  }

  free(rawScores);
  free(intermediateScores);
  return true;
}
bool scoreCompare(const std::pair<int, float>& firstElem, const std::pair<int, float>& secondElem) {
  return firstElem.second > secondElem.second;
}

void CompleteGraphScorer::BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const {
  
  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(scores, pairs, scoreCompare);
  
  out << std::endl << "Top 10 genes" << std::endl << "Gene\tScore" << std::endl << "----------------" << std::endl;
  int i = 0;
  for (auto it = pairs.begin(); it != pairs.end() && i < 10; it++, i++) {
    out << rmap[it->first] << "\t" << it->second << std::endl;
  }
}

void CompleteGraphScorer::LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const {

  int l = 1;
  for (auto const &g : groups) {
    out << "----------------" << std::endl;
    out << "Locus " << l << std::endl << "----------------" << std::endl;
    TScoreMap locusScores;
    for (auto const &e : g) {
      locusScores[e] = scores[e];
    }
    std::vector<std::pair<int, float> > pairs;
    sortMapByVal(locusScores, pairs, scoreCompare);
    for (auto const &e : pairs) {
      out << rmap[e.first] << "\t" << e.second << std::endl;
    }
    l++;
  }
  
}

bool CompleteGraphScorer4::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, TScoreMap& scores) const
{
  // Allocate memory for scores
  float* const rawScores = (float*)malloc(sizeof(float) * width);
  if (rawScores == 0) {
    printf("Could not allocate buffer of %d bytes\n", (int)sizeof(float) * width);
    return false;
  }

  float* const intermediateScores = (float*)malloc(sizeof(float) * width);
  if (intermediateScores == 0) {
    printf("Could not allocate buffer of %d bytes\n", (int)sizeof(float) * width);
    return false;
  }

  for (int i = 0; i < width; ++i) {
    rawScores[i] = 0.0f;
    intermediateScores[i] = 0.0f;
  }

  printf("Starting complete graph 4 scoring.\n");
  
  const int gsize = 4;
  const int numGroups = groups.size();
  
  for (int i = 0; i <= numGroups - gsize; ++i) {
    const TInts g1 = groups[i];
    for (int j = i + 1; j <= numGroups - gsize + 1; ++j) {
      const TInts g2 = groups[j];
      for (int k = j + 1; k <= numGroups - gsize + 2; ++k) {
	const TInts g3 = groups[k];
	for (int l = k + 1; l <= numGroups - gsize + 3; ++l) {
	  const TInts g4 = groups[l];
	  for (const auto i1 : g1) {
	    const float* r1 = similarities + width * i1;
	    for (const auto i2 : g2) {
	      const float* r2 = similarities + width * i2;
	      const float v12 = *(r1 + i2);
	      for (const auto i3 : g3) {
		const float* r3 = similarities + width * i3;
		const float v123 = v12 + *(r1 + i3) + *(r2 + i3);

		for (const auto i4 : g4) {
		  // i1, i2, i3
		  float acc = v123;
		  // i1 - i4
		  acc += *(r1 + i4);
		  // i2 - i4
		  acc += *(r2 + i4);
		  // i3 - i4
		  acc += *(r3 + i4);
		  
		  if (acc > intermediateScores[i1]) intermediateScores[i1] = acc;
		  if (acc > intermediateScores[i2]) intermediateScores[i2] = acc;
		  if (acc > intermediateScores[i3]) intermediateScores[i3] = acc;
		  if (acc > intermediateScores[i4]) intermediateScores[i4] = acc;
		}
	      }
	    }
	  }
	  for (int i = 0; i < width; ++i) {
	    rawScores[i] += intermediateScores[i];
	    intermediateScores[i] = 0.0f;
	  }
	}
      }
    }
  }

  for (auto const &g : groups) {
    for (auto i : g) {
      scores[i] = rawScores[i];
    }
  }

  printf("Finishing complete graph 4 scoring.\n");
  
  free(rawScores);
  free(intermediateScores);
  return true;
}
