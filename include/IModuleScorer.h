/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** IModuleScorer.h
** Interface for generic prioritization scorer.
**
** Author: kca
** -------------------------------------------------------------------------*/

#ifndef IMODULESCORER_H
#define IMODULESCORER_H

#include <vector>
#include <map>
#include <string>

typedef std::vector<int> TIndices;
//typedef std::vector<TIndices> TIndicesGroups;
typedef std::map<std::string, TIndices> TIndicesGroups;
typedef std::map<int, float> TScoreMap;
typedef std::map<int, std::string> TReverseIndexMap;

class IModuleScorer {
 public:
  virtual ~IModuleScorer() {}
  virtual bool ScoreModule(const float* const similarities, const int width,
			   const TIndicesGroups& groups, 
			   const TIndices& indicesToScore,
			   TScoreMap& scores) const = 0;
  virtual void BriefSummary(TScoreMap& scores, TReverseIndexMap& rmap, std::ostream& outstream) const = 0;
  virtual void LongSummary(TScoreMap& scores, TReverseIndexMap& rmap, const TIndicesGroups& groups, std::ostream& out) const = 0;
};

#endif
