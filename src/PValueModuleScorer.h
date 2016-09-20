/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** PValueModuleScorer.h
** This header declares the CompleteGraphScorer.
**
** Author: kca
** -------------------------------------------------------------------------*/

#include "../include/IModuleScorer.h"

class PValueModuleScorer : public IModuleScorer {
 public:
  // I take ownership of scorer
 PValueModuleScorer(const int numIterations, IModuleScorer* scorer, std::map<int, int>& nodeDegrees)
   : mNumIterations(numIterations), mScorer(scorer), mNodeDegreeMap(nodeDegrees)
  {
    for (auto const& p : nodeDegrees) {
      mDegreeGroupNodes[p.second].push_back(p.first);
    }
  }
  ~PValueModuleScorer(void) {
    delete mScorer;
  }
  
  bool ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups,
		   const TIndices& indicesToScore, TScoreMap& scores) const;
  void BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const;
  void LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const;
  void ShuffleGroups(const TIndicesGroups& groups, const std::vector<int>& excisedIndicies, TIndicesGroups& shuffledGroups) const;
  
 private:
  const int mNumIterations;
  IModuleScorer* mScorer;
  const std::map<int, int>& mNodeDegreeMap;
  std::map<int, std::vector<int> > mDegreeGroupNodes;
};
