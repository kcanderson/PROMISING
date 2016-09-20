/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** FastScorer.h
** This header declares the CompleteGraphScorer.
**
** Author: kca
** -------------------------------------------------------------------------*/

#include "../include/IModuleScorer.h"

class FastScorer : public IModuleScorer {
public:
  FastScorer(void) {}
  bool ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups,
		   const TIndices& indicesToScore, TScoreMap& scores) const;
  void BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const;
  void LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const;
};


