#ifndef FISH_HOOK_INTERVAL_H__
#define FISH_HOOK_INTERVAL_H__

#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamHeader.h"
#include "Fractions.h"

class FishHookInterval : public SeqLib::GenomicRegion {

 public:
 
  /** Create an empty interval */
  FishHookInterval() {}

  /** Create an interval */
  FishHookInterval(int c, int p, int e) : SeqLib::GenomicRegion(c, p, e, '*') {}

  /** Construct a fish hook interval */
 FishHookInterval(const std::string& tchr, const std::string& tpos1, const std::string& tpos2, const SeqLib::BamHeader& hdr) : 
    SeqLib::GenomicRegion(tchr, tpos1, tpos2, hdr) {}

    /** Add a covariate (eg SINE density) by name and its value */
    void AddCovariate(const std::string& name, double val) {
      m_hash_table[name] = val;
    }
    
 private:
    
  // frac covered or score in this interval
  SeqHashMap<std::string, double> m_hash_table;

};

class FishHookTiles : public SeqLib::GenomicRegionCollection<FishHookInterval> {

 public:

  FishHookTiles(int width, int ovlp, const SeqLib::HeaderSequenceVector& h);

  void AddIntervalCovariate(const std::string& name, const Fractions& f);

  void AddScoreCovariate(const std::string& name, const Fractions& f);

};


#endif
