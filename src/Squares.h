#ifndef SWAP_SQUARES_H__
#define SWAP_SQUARES_H__

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "Matrix.h"

class Matrix; 

/** Container to hold a bunch of regions on the 2D matrix to query
 */
class Squares 
{
  
 public:
  Squares() {}

  void checkForOverlaps(const Matrix& m);

 private:

  SeqLib::GRC m_grv_a;
  SeqLib::GRC m_grv_b;

  // indexed by a*size(a) + b
  std::vector<size_t> m_overlaps;
  
};

#endif
