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
 PValueModuleScorer(const int numIterations, IModuleScorer* scorer)
   : mNumIterations(numIterations), mScorer(scorer)
  {}
  ~PValueModuleScorer(void) {
    delete mScorer;
  }
  
  bool ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups,
		   TScoreMap& scores) const;
  void BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const;
  void LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const;
  
 private:
  const int mNumIterations;
  IModuleScorer* mScorer;
};
