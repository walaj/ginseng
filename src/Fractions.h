#ifndef SNOWMAN_FRACTIONS_H__
#define SNOWMAN_FRACTIONS_H__

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamHeader.h"
#include <string>

  /** @brief Extension of GenomicRegion with a fraction of reads to keep on an interval
   * 
   * Used in conjunction with Fractions and used for selectively sub-sampling BAM files.
   */
class FracRegion : public SeqLib::GenomicRegion {

 public:
  
  FracRegion() {}

  FracRegion(const std::string& c, const std::string& p1, const std::string& p2, const SeqLib::BamHeader& h, const std::string& f);

  friend std::ostream& operator<<(std::ostream& out, const FracRegion& f);

  double frac;
};

  /** @brief Genomic intervals and associated sampling fractions
   */
class Fractions : public SeqLib::GenomicRegionCollection<FracRegion> {

 public:
  
  // make an empty one
 Fractions() : SeqLib::GenomicRegionCollection<FracRegion>() {}

  // construct tiled across genome
  Fractions(int w, int o, const SeqLib::HeaderSequenceVector& v) : SeqLib::GenomicRegionCollection<FracRegion>(w, o, v) {}

  void readFromBed(const std::string& file, const SeqLib::BamHeader& h);

 private:

};

#endif
