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
#include <cfloat>
#include <cmath>

typedef std::vector<int> TInts;

bool CompleteGraphScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, const TIndices& indicesToScore, TScoreMap& scores) const
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
  std::vector<TInts> grps;
  for (auto const &group : groups) {
    grps.push_back(group.second);
  }

  for (int i = 0; i <= numGroups - gsize; ++i) {
    const TInts g1 = grps[i];
    for (int j = i + 1; j <= numGroups - gsize + 1; ++j) {
      const TInts g2 = grps[j];
      for (int k = j + 1; k <= numGroups - gsize + 2; ++k) {
	const TInts g3 = grps[k];

	for (const auto i1 : g1) {
	  const float* r1 = similarities + width * i1;
	  for (const auto i2 : g2) {
	    const float* r2 = similarities + width * i2;
	    const float v12 = *(r1 + i2);
	    for (const auto i3 : g3) {
	      const float v13 = *(r1 + i3);
	      const float v23 = *(r2 + i3);
#define MIN_SCORING 1
#if MIN_SCORING
	      // Score i1
	      float acc = v12 + v13;
	      float mn = (v12 > v13) ? v12 : v13;
	      acc += (v23 > mn) ? mn : v23;
	      if (acc > intermediateScores[i1]) intermediateScores[i1] = acc;

	      // Score i2
	      acc = v12 + v23;
	      //mx = (v12 > v23) ? v12 : v23;
	      mn = (v12 > v23) ? v12 : v23;
	      acc += (v13 > mn) ? mn : v13;
	      if (acc > intermediateScores[i2]) intermediateScores[i2] = acc;

	      // Score i3
	      acc = v13 + v23;
	      //mx = (v13 > v23) ? v13 : v23;
	      mn = (v13 > v23) ? v13 : v23;
	      acc += (v12 > mn) ? mn : v12;
	      if (acc > intermediateScores[i3]) intermediateScores[i3] = acc;
#else
	      //i1 - i2
	      float acc = v12;
	      //i1 - i3
	      acc += *(r1 + i3);
	      acc += r1[i3];
	      //i2 - i3
	      acc += *(r2 + i3);
	      acc += r2[i3];

	      if (acc > intermediateScores[i1]) intermediateScores[i1] = acc;
	      if (acc > intermediateScores[i2]) intermediateScores[i2] = acc;
	      if (acc > intermediateScores[i3]) intermediateScores[i3] = acc;
#endif
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
    for (auto i : g.second) {
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

void CompleteGraphScorer::LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const
{

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
  sortMapByVal(scores, pairs, scoreCompare);
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
  //   sortMapByVal(locusScores, pairs, scoreCompare);
  //   for (auto const &e : pairs) {
  //     out << rmap[e.first] << "\t" << e.second << std::endl;
  //   }
  //   l++;
  // }
  
}

bool CompleteGraphScorer4::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, const TIndices& indicesToScore, TScoreMap& scores) const
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
  std::vector<TInts> grps;
  for (auto const &group : groups) {
    grps.push_back(group.second);
  }
  
  for (int i = 0; i <= numGroups - gsize; ++i) {
    const TInts g1 = grps[i];
    for (int j = i + 1; j <= numGroups - gsize + 1; ++j) {
      const TInts g2 = grps[j];
      for (int k = j + 1; k <= numGroups - gsize + 2; ++k) {
	const TInts g3 = grps[k];
	for (int l = k + 1; l <= numGroups - gsize + 3; ++l) {
	  const TInts g4 = grps[l];
	  for (const auto i1 : g1) {
	    const float* r1 = similarities + width * i1;
	    for (const auto i2 : g2) {
	      const float* r2 = similarities + width * i2;
	      const float v12 = *(r1 + i2);
	      for (const auto i3 : g3) {
		const float* r3 = similarities + width * i3;
		const float v13 = r1[i3];
		const float v23 = r2[i3];
		const float v123 = v12 + v13 + v23;

		for (const auto i4 : g4) {
#if MIN_SCORING
		  const float v14 = r1[i4];
		  const float v24 = r2[i4];
		  const float v34 = r3[i4];
		  
		  // i1
		  float mn = (v12 < v13) ? v12 : v13;
		  mn = (mn < v14) ? mn : v14;
		  float acc = v12 + v13 + v14;
		  acc += (v23 > mn) ? mn : v23;
		  acc += (v24 > mn) ? mn : v24;
		  acc += (v34 > mn) ? mn : v34;
		  if (acc > intermediateScores[i1]) intermediateScores[i1] = acc;

		  // i2
		  mn = (v12 < v23) ? v12 : v23;
		  mn = (mn < v24) ? mn : v24;
		  acc = v12 + v23 + v24;
		  acc += (v13 > mn) ? mn : v13;
		  acc += (v14 > mn) ? mn : v14;
		  acc += (v34 > mn) ? mn : v34;
		  if (acc > intermediateScores[i2]) intermediateScores[i2] = acc;

		  // i3
		  mn = (v13 < v23) ? v13 : v23;
		  mn = (mn < v34) ? mn : v34;
		  acc = v13 + v23 + v34;
		  acc += (v12 > mn) ? mn : v12;
		  acc += (v14 > mn) ? mn : v14;
		  acc += (v23 > mn) ? mn : v23;
		  if (acc > intermediateScores[i3]) intermediateScores[i3] = acc;

		  // i4
		  mn = (v14 < v24) ? v14 : v24;
		  mn = (mn < v34) ? mn : v34;
		  acc = v14 + v24 + v34;
		  acc += (v12 > mn) ? mn : v12;
		  acc += (v13 > mn) ? mn : v13;
		  acc += (v23 > mn) ? mn : v23;
		  if (acc > intermediateScores[i4]) intermediateScores[i4] = acc;
#else
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
#endif
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
    for (auto i : g.second) {
      scores[i] = rawScores[i];
    }
  }

  printf("Finishing complete graph 4 scoring.\n");
  
  free(rawScores);
  free(intermediateScores);
  return true;
}

float score_complete4(const int score_node, const std::vector< std::vector<int> >& other_groups, const float* const similarities, const int width)
{
  const float* const me_base = similarities + (width * score_node);
  std::vector< std::vector<int> > combinations;
  comb(other_groups.size(), 3, combinations);
  float score(0.0f);
  
  for (auto const& combination : combinations) {
    float maxscore = -FLT_MAX;
    const auto ig1 = combination[0];
    const auto g1 = other_groups[ig1];
    const auto ig2 = combination[1];
    const auto g2 = other_groups[ig2];
    const auto ig3 = combination[2];
    const auto g3 = other_groups[ig3];

    for (const auto n1 : g1) {
      const float me_n1 = me_base[n1];
      const float* n1_base = similarities + (width * n1);
      const float curr_outer(me_n1);
      for (const auto n2 : g2) {
	const float me_n2 = me_base[n2];
	const float* n2_base = similarities + (width * n2);
	const float n1_n2 = n1_base[n2];
	const float curr_middle(curr_outer + me_n2);
	for (const auto n3 : g3) {
	  const float me_n3 = me_base[n3];
	  float curr(curr_middle + me_n3);
	  const float n2_n3 = n2_base[n3];
	  const float n1_n3 = n1_base[n3];
	  float v((me_n1 < me_n2) ? me_n1 : me_n2);
	  v = (v < n1_n2) ? v : n1_n2;
	  //curr += v;
	  curr += n1_n2;
	  v = (me_n1 < me_n3) ? me_n1 : me_n3;
	  v = (v < n1_n3) ? v : n1_n3;
	  //curr += v;
	  curr += n1_n3;
	  v = (me_n2 < me_n3) ? me_n2 : me_n3;
	  v = (v < n2_n3) ? v : n2_n3;
	  //curr += v;
	  curr += n2_n3;
	  maxscore = (curr > maxscore) ? curr : maxscore;
	}
      }
    }
    if (maxscore > -FLT_MAX) score += maxscore;
  }
  return score;
}

// assignment into a
#define MIN(a, b, c)				\
  a = a;
//  a = (a < b) ? a : b;			\
//  a = (a < c) ? a : c;

float score_complete5(const int score_node, const std::vector< std::vector<int> >& other_groups, const float* const similarities, const int width)
{
  const float* const me_base = similarities + (width * score_node);
  std::vector< std::vector<int> > combinations;
  comb(other_groups.size(), 4, combinations);
  float score(0.0f);
  
  for (auto const& combination : combinations) {
    float maxscore = -FLT_MAX;
    const auto ig1 = combination[0];
    const auto g1 = other_groups[ig1];
    const auto ig2 = combination[1];
    const auto g2 = other_groups[ig2];
    const auto ig3 = combination[2];
    const auto g3 = other_groups[ig3];
    const auto ig4 = combination[3];
    const auto g4 = other_groups[ig4];

    for (const auto n1 : g1) {
      const float me_n1 = me_base[n1];
      const float* n1_base = similarities + (width * n1);
      const float curr_outer(me_n1);
      for (const auto n2 : g2) {
	const float me_n2 = me_base[n2];
	const float* n2_base = similarities + (width * n2);
	// n1_n2
	const float n1_n2 = n1_base[n2];
	float v(n1_n2);
	MIN(v, me_n1, me_n2);
	float curr_middle(curr_outer + me_n2 + v);
	for (const auto n3 : g3) {
	  const float me_n3 = me_base[n3];
	  const float* n3_base = similarities + (width * n3);
	  // n1_n3
	  const float n1_n3 = n1_base[n3];
	  v = n1_n3;
	  MIN(v, me_n1, me_n3);
	  // n2_n3
	  const float n2_n3 = n2_base[n3];
	  float vv(n2_n3);
	  MIN(vv, me_n2, me_n3);
	  
	  const float curr_middle2(curr_middle + me_n3 + v + vv);
	  
	  for (const auto n4 : g4) {
	    const float me_n4 = me_base[n4];
	    float curr(curr_middle2 + me_n4);
	    // n1_n4
	    v = n1_base[n4];
	    MIN(v, me_n1, me_n4);
	    curr += v;
	    // n2_n4
	    v = n2_base[n4];
	    MIN(v, me_n2, me_n4);
	    curr += v;
	    // n3_n4
	    v = n2_base[n4];
	    MIN(v, me_n3, me_n4);
	    curr += v;
	    
	    maxscore = (curr > maxscore) ? curr : maxscore;
	  }
	}
      }
    }
    if (maxscore > -FLT_MAX) score += maxscore;
  }
  
  return score;
}


float score_complete3(const int score_node, const std::vector< std::vector<int> >& other_groups, const float* const similarities, const int width)
{
  const float* const me_base = similarities + (width * score_node);
  std::vector< std::vector<int> > combinations;
  comb(other_groups.size(), 2, combinations);
  float score(0.0f);
  
  for (auto const& combination : combinations) {
    float maxscore = -FLT_MAX;
    const auto ig1 = combination[0];
    const auto g1 = other_groups[ig1];
    const auto ig2 = combination[1];
    const auto g2 = other_groups[ig2];

    for (const auto n1 : g1) {
      const float me_n1 = me_base[n1];
      const float* n1_base = similarities + (width * n1);
      const float curr_outer(me_n1);
      for (const auto n2 : g2) {
	const float me_n2 = me_base[n2];
	const float* n2_base = similarities + (width * n2);
	const float n1_n2 = n1_base[n2];
	float curr(curr_outer + me_n2);
	float v = (me_n1 < me_n2) ? me_n1 : me_n2;
	v = (v < n1_n2) ? v : n1_n2;
	//float v = (n1_n2 < me_n1) ? n1_n2 : me_n1;
	//v = (v < me_n2) ? v : me_n2;
	//curr += v;
	curr += n1_n2;
	maxscore = (curr > maxscore) ? curr : maxscore;
      }
    }
    if (maxscore > -FLT_MAX) score += maxscore;
  }
  return score;
}


void maximum_similarities_in_each_group(const int me, const float* const similarities, const int width, const std::vector< std::vector<int> >& other_groups, std::map<int, float>& max_values)
{
  const float* base_sim(similarities + (width * me));
  for (int i = 0; i < other_groups.size(); ++i) {
    const auto other_group = other_groups[i];
    float mx(-FLT_MAX);
    for (const auto other_node : other_group) {
      const float v(*(base_sim + other_node));
      mx = (v > mx) ? v : mx;
    }
    max_values[i] = mx;
  }
}

void top_groups_for_candidate(const int candidate, const float* const similarities, const int width, const std::vector< std::vector<int> >& other_groups, const int num_groups_to_consider, std::vector< std::vector<int> >& top_groups)
{
  std::map<int, float> max_values;
  maximum_similarities_in_each_group(candidate, similarities, width, other_groups, max_values);
  std::vector<std::pair<int, float> > pairs;
  sortMapByVal(max_values, pairs, scoreCompare);
  const int num_groups = std::min({(int)(pairs.size()), num_groups_to_consider});
  for (int i = 0; i < num_groups; ++i) {
    top_groups.push_back(other_groups[pairs[i].first]);
  }
}

#define NUMGROUPSTOCONSIDER3 150
#define NUMGROUPSTOCONSIDER4 100
#define NUMGROUPSTOCONSIDER5 35

bool CompleteGraphFasterScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, const TIndices& indicesToScore, TScoreMap& scores) const
{
  // For each group, g_score, to score
  //  For each node, n_score, in g_score
  //   Figure out N best other groups, g_others, among all_groups - g_score
  //
  std::vector< std::vector<int> > all_groups;
  for (const auto& g : groups) {
    all_groups.push_back(g.second);
  }
  
  for (const auto& g : groups) {
    const auto my_group = g.second;
    std::vector< std::vector<int> > others(all_groups);
    auto it = std::find(others.begin(), others.end(), my_group);
    others.erase(it);
    printf("group %s\n", g.first.c_str());
    for (const auto me : my_group) {
      // Find the strongest hits in each locus
      std::vector< std::vector<int> > top_groups;
      if (this->mScoreSize ==  3 || groups.size() < 4) {
	top_groups_for_candidate(me, similarities, width, others, NUMGROUPSTOCONSIDER3, top_groups);
	scores[me] = score_complete3(me, top_groups, similarities, width);
      }
      else if (this->mScoreSize == 4 or groups.size() < 5) {
	top_groups_for_candidate(me, similarities, width, others, NUMGROUPSTOCONSIDER4, top_groups);
	scores[me] = score_complete4(me, top_groups, similarities, width);
      }
      else {
	top_groups_for_candidate(me, similarities, width, others, NUMGROUPSTOCONSIDER5, top_groups);
	scores[me] = score_complete5(me, top_groups, similarities, width);
      }
    }
  }
  return true;
}

void randomize_matrix(float* const similarities, const int width, const int numShuffles)
{
  for (int i = 0; i < numShuffles; ++i) {
    const int i1 = rand() % width;
    const int i2 = rand() % width;
    const int j1 = rand() % width;
    const int j2 = rand() % width;
    const float v1 = similarities[width * i1 + i2];
    const float v2 = similarities[width * j1 + j2];
    similarities[width * i1 + i2] = v2;
    similarities[width * i2 + i1] = v2;
    similarities[width * j1 + j2] = v1;
    similarities[width * j2 + j1] = v1;
  }
}


bool PValCompleteScorer::ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, const TIndices& indicesToScore, TScoreMap& scores) const
{
  CompleteGraphFasterScorer scorer(3);

  TScoreMap myscores;
  scorer.ScoreModule(similarities, width, groups, indicesToScore, myscores);

  // Copy matrix
  const int numBytes(sizeof(float) * width * width);
  float* const scratch((float*)malloc(numBytes));
  if (scratch == 0) {
    printf("Could not allocate buffer of %i bytes.\n", numBytes);
    return false;
  }
  memcpy(scratch, similarities, numBytes);

  std::map<int, std::vector<float> > scoresMap;
  for (const auto &g : groups) {
    for (auto i : g.second) {
      //scoresMap[i].push_back(1.0);
    }
  }
  const int numIterations(30);
  for (int i = 0; i < numIterations; ++i) {
    printf("iteration %i\n", i);
    randomize_matrix(scratch, width, (width * width) / 2);
    TScoreMap permuted_scores;
    scorer.ScoreModule(scratch, width, groups, indicesToScore, permuted_scores);
    for (const auto &score : permuted_scores) {
      scoresMap[score.first].push_back(score.second);
    }
  }

  for (const auto &permutations : scoresMap) {
    const float real_score(myscores[permutations.first]);
    //int acc(0);
    float tot(0.0f);
    for (const auto s : permutations.second) {
      //if (real_score > s) acc++;
      tot += s;
    }
    const float m = tot / permutations.second.size();
    tot = 0.0f;
    for (const auto s : permutations.second) {
      //if (real_score > s) acc++;
      tot += fabs(s - m);
    }
    const float std = tot / permutations.second.size();
    
    scores[permutations.first] = (real_score - m) / std;
  }
  
  return true;
}

