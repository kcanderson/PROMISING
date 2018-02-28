/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** FastScorer.cpp
** This file implements the FastScorer class.
**
** Author: kca
** -------------------------------------------------------------------------*/

#include "FastScorer.h"
#include "coreroutines.h"
#include "string.h"
#include <algorithm>
#include <cfloat>

typedef std::vector<int> TInts;

bool FastScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, const TIndices& indicesToScore, TScoreMap& scores) const
{
  // For each group, g
  //   For each node in g, n
  //     best_nodes = [g]
  //     For each other group, g_other
  //       max_node = nil; max_score = 0
  //       For each node in g_other, n_other
  //         if similarities(n, n_other) > max_score: update max_node, max_score
  //     append max_node to best_nodes
  //     score_g = 0
  //     for i in 0:length(best_nodes)-1
  //       for j in i+1:length(best_nodes)-1
  //         score_g += similarities(best_nodes[i], best_nodes[j])

  
  TIndicesGroups goodGroups;
  for (auto const& group : groups) {
    // intersection between indicesToScore and group not nil?
    TIndices g;
    for (auto i : group.second) {
      if (std::find(indicesToScore.begin(), indicesToScore.end(), i) != indicesToScore.end()) {
	g.push_back(i);
      }
    }
    if (g.size() > 0) {
      //goodGroups.push_back(g);
      goodGroups[group.first] = g;
    }
  }
  
  const int numGroups = groups.size();
  
  int* const indexBuffer = (int*)malloc(sizeof(int) * numGroups);
  if (indexBuffer == 0) {
    printf("Could not allocate buffer of %d bytes\n", (int)sizeof(int) * numGroups);
    return false;
  }

  for (auto const& group : goodGroups) {
    for (auto i : group.second) {
      scores[i] = 0.0f;
    }
  }
  
  for (auto const& group : goodGroups) {
    //for (int i = 0; i < group.size(); ++i) {
    //const auto me = group[i];
    for (auto const& me : group.second) {
      const float* const sim = similarities + width * me;
      TIndicesGroups otherGroups(groups);
      auto it = std::find(otherGroups.begin(), otherGroups.end(), group);
      otherGroups.erase(it);
      //printf("%i, %i\n", groups.size(), otherGroups.size());
      int numothers = 0;
      //indexBuffer[numothers++] = me;
      float score = 0.0f;
      for (auto const& othergroup : otherGroups) {
  	int maxnode = -1;
  	float maxscore = -FLT_MAX;
  	//for (const auto j : othergroup) {
	//for (int j = 0; j < othergroup.size(); ++j) {
	//const int ot = othergroup[j];
	for (const auto& ot : othergroup.second) {
  	  const float sc = sim[ot];
	  //maxnode = (sc > maxscore) ? ot : maxnode;
	  //maxscore = (sc > maxscore) ? sc : maxscore;
  	  if (sc > maxscore) {
  	    maxnode = ot;
  	    maxscore = sc;
  	  }
	}
  	if (maxnode != -1) {
  	  indexBuffer[numothers++] = maxnode;
	  score += maxscore;
  	}
      }
      
      const int endloop = numothers;
      const int* cindexBuffer = indexBuffer;
      for (int ii = 0; ii < endloop; ++ii) {
      	const int iii = cindexBuffer[ii];
      	const float* base_iii = similarities + width * iii;
	const float me_iii(sim[iii]);
      	for (int jj = ii + 1; jj < endloop; ++jj) {
      	  const int jjj = cindexBuffer[jj];
	  const float me_jjj(sim[jjj]);
      	  const float iii_jjj(base_iii[jjj]);
	  float val = (iii_jjj < me_iii) ? iii_jjj : me_iii;
	  val = (val < me_jjj) ? val : me_jjj;
	  //float val(MIN(iii_jjj, MIN(me_iii, me_jjj)));
      	  //float b = sim[jjj];
	  //printf("%i: %i %i= %f\n", me, iii, jjj, vvv);
      	  //b = (a < b) ? a : b;
      	  //b = (b < vvv) ? b : vvv;
      	  //score += b;
	  score += val;
      	}
      }
      scores[me] = score;
    }
  }

  free(indexBuffer);

  return true;
}
bool scoreCompare2(const std::pair<int, float>& firstElem, const std::pair<int, float>& secondElem) {
  return firstElem.second > secondElem.second;
}

void FastScorer::BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const {
  
  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(scores, pairs, scoreCompare2);
  
  out << std::endl << "Top 10 genes" << std::endl << "Gene\tScore" << std::endl << "----------------" << std::endl;
  int i = 0;
  for (auto it = pairs.begin(); it != pairs.end() && i < 10; it++, i++) {
    out << rmap[it->first] << "\t" << it->second << std::endl;
  }
}

void FastScorer::LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const {

#define OLD_FORMAT 0
#if OLD_FORMAT
  int l = 1;
  for (auto const &g : groups) {
    out << "----------------" << std::endl;
    out << "Locus " << l << std::endl << "----------------" << std::endl;
    TScoreMap locusScores;
    for (auto const &e : g.second) {
      locusScores[e] = scores[e];
    }
    std::vector<std::pair<int, float> > pairs;
    sortMapByVal(locusScores, pairs, scoreCompare);
    for (auto const &e : pairs) {
      out << rmap[e.first] << "\t" << e.second << std::endl;
    }
    l++;
  }
#else
  std::map<int, std::string> groupMap;
  int i = 0;
  for (auto const &g : groups) {
    //std::string locus = "Locus " + std::to_string(i++);
    for (auto gene : g.second) {
      groupMap[gene] = g.first;
    }
  }

  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(scores, pairs, scoreCompare2);
  out << "locus\tgene\tscore" << std::endl;
  for (auto const &e : pairs) {
    out << groupMap[e.first] << "\t" << rmap[e.first] << "\t" << scores[e.first] << std::endl;
  }
#endif

  
  // int l = 1;
  // for (auto const &g : groups) {
  //   out << "----------------" << std::endl;
  //   out << "Locus " << l << std::endl << "----------------" << std::endl;
  //   TScoreMap locusScores;
  //   for (auto const &e : g.second) {
  //     locusScores[e] = scores[e];
  //   }
  //   std::vector<std::pair<int, float> > pairs;
  //   sortMapByVal(locusScores, pairs, scoreCompare2);
  //   for (auto const &e : pairs) {
  //     out << rmap[e.first] << "\t" << e.second << std::endl;
  //   }
  //   l++;
  // }
  
}
