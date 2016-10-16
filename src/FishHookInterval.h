#ifndef FISH_HOOK_INTERVAL_H__
#define FISH_HOOK_INTERVAL_H__

#include <map>

#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamHeader.h"
#include "Fractions.h"
#include "Events.h"

class FishHookTiles;

class FishHookInterval : public SeqLib::GenomicRegion {

 public:
 
  friend class FishHookTiles;

  /** Create an interval from a simple GenomicRegion */
 FishHookInterval(const SeqLib::GenomicRegion& gr) : SeqLib::GenomicRegion(gr) {}

  /** Create an empty interval */
  FishHookInterval() {}

  /** Create an interval */
  FishHookInterval(int c, int p, int e) : SeqLib::GenomicRegion(c, p, e, '*') {}

  /** Construct a fish hook interval */
 FishHookInterval(const std::string& tchr, const std::string& tpos1, const std::string& tpos2, const SeqLib::BamHeader& hdr) : 
    SeqLib::GenomicRegion(tchr, tpos1, tpos2, hdr) {}

    /** Add a covariate (eg SINE density) by name and its value */
    void AddCovariate(const std::string& name, float val) {
      m_hash_table[name] += val;
    }

    void ScaleCovariate(const std::string& name, float scale) {
      m_hash_table[name] *= scale;
    }
    
    /** Print a single record */
    void PrintBED(std::ostream& out, const SeqLib::BamHeader& h) const;

    float covered = 0; ///< fraction of this interval that is covered by positive mask (eg mapq)
    
    int events = 0; ///< number of events in this region

    size_t NumCovariates() const { return m_hash_table.size(); }
    
    typename std::map<std::string, float>::iterator begin() { return m_hash_table.begin(); } 
    
    typename std::map<std::string, float>::iterator end() { return m_hash_table.end(); } 
    
    typename std::map<std::string, float>::const_iterator begin() const { return m_hash_table.begin(); } 
    
    typename std::map<std::string, float>::const_iterator end() const { return m_hash_table.end(); } 

 private:

  // frac covered or score in this interval
  std::map<std::string, float> m_hash_table;

};

class FishHookTiles : public SeqLib::GenomicRegionCollection<FishHookInterval> {

 public:

  FishHookTiles(int width, int ovlp, const SeqLib::HeaderSequenceVector& h);

  void AddIntervalCovariate(const std::string& name, const Fractions& f);

  void PrintBED(std::ostream& out, const SeqLib::BamHeader& h) const;

  size_t NumCovariates() const;

  /** Print a the BED header */
  void PrintBEDHeader(std::ostream& out) const;

  /** Tally number of events per tile */
  void CountEvents(const EventList& events);
    
};


#endif
