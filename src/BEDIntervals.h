#ifndef GINSENG_BEDINTERVALS_H__
#define GINSENG_BEDINTERVALS_H__

#include "SeqLib/GenomicRegionCollection.h"
#include <unordered_map>

typedef std::unordered_map<int, int> PosBin;
typedef std::unordered_map<int, PosBin> ChrPosBin;

class BEDIntervals {
  
 public:

  size_t size() const { return grc.size(); }
  
  /*
  void AddHits(const SeqLib::GenomicRegion& g);

  // query if a 2D point is in this region (1) or even more specifically, 
  // if they are both in the exact same interval (2). Otherwise 0
  int QueryPair(const SeqLib::GenomicRegion* r, const SeqLib::GenomicRegion* c) const;

  int QueryElem(const SeqLib::GenomicRegion* g) const;
  */

  // hold the intervals
  SeqLib::GRC grc;  

 private:

  // chr (pos, binid)
  //std::unordered_map<int, std::unordered_map<int, int>> m_hits; // chr, pos

    
};


#endif
