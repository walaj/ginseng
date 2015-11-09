#ifndef SWAP_SQUARES_H__
#define SWAP_SQUARES_H__

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"
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

  SnowTools::GRC m_grv_a;
  SnowTools::GRC m_grv_b;

  // indexed by a*size(a) + b
  std::vector<size_t> m_overlaps;
  
};

#endif
