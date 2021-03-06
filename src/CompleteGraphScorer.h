/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** CompleteGraphScorer.h
** This header declares the CompleteGraphScorer.
**
** Author: kca
** -------------------------------------------------------------------------*/

#include "../include/IModuleScorer.h"

/* class CompleteGraphScorer : public IModuleScorer { */
/* public: */
/*  CompleteGraphScorer(bool clamp): mClamp(clamp) {} */
/*   bool ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, */
/* 		   const TIndices& indicesToScore, TScoreMap& scores) const; */
/*   void BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const; */
/*   void LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const; */
/*  private: */
/*   bool mClamp; */
/* }; */

/* class CompleteGraphScorer4 : public CompleteGraphScorer { */
/* public: */
/*  CompleteGraphScorer4(bool clamp): CompleteGraphScorer(clamp) {} */
/*   bool ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups, */
/* 		   const TIndices& indicestoScore, TScoreMap& scores) const; */
/*   //void BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const; */
/*   //void LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const; */
/* }; */

class BaseScorer: public IModuleScorer {
 public:
  void BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& out) const;
  void LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const;
};


class CompleteGraphFasterScorer : public BaseScorer {
 public:
 CompleteGraphFasterScorer(const int scoreSize, const bool clamp = false): mScoreSize(scoreSize), mClamp(clamp) {}
  bool ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups,
		   const TIndices& indicestoScore, TScoreMap& scores) const;
 private:
  int mScoreSize;
  bool mClamp;
};

class PValCompleteScorer : public BaseScorer {
 public:
  PValCompleteScorer() {}
  bool ScoreModule(const float* const similarities, const int width, const TIndicesGroups& groups,
		   const TIndices& indicestoScore, TScoreMap& scores) const;
};
