#include "FishHookInterval.h"

FishHookTiles::FishHookTiles(int width, int ovlp, const SeqLib::HeaderSequenceVector& h) : GenomicRegionCollection<FishHookInterval>(width, ovlp, h) {
 
  /*  
  size_t chr = 0;
  for (auto& i : h) {
    size_t start = 0; 
    for (; start + width < i.Length; start += width) {
      add(FishHookInterval(chr, start, start + width));
    }
    ++chr;
    }*/

}

void FishHookTiles::PrintBED(std::ostream& out, const SeqLib::BamHeader& h) const {
  for (const auto& i : *this)
    i.PrintBED(out, h);
}

void FishHookInterval::PrintBED(std::ostream& out, const SeqLib::BamHeader& h) const {
  out << h.IDtoName(chr) << "\t" << pos1 << "\t" << pos2 << "\t" << covered << "\t" << events;
  for (const auto& i : m_hash_table)
    out << "\t" << i.second;
  out << std::endl;
}

void FishHookTiles::PrintBEDHeader(std::ostream& out) const {
  
  for (const auto& i : *this) {
    out << "#CHROM\tPOS1\tPOS2\tcovered\tevents";
    for (const auto& j : i.m_hash_table)
      out << "\t" << j.first;
    out << std::endl;
    return;
  }
}

void FishHookTiles::AddIntervalCovariate(const std::string& name, const Fractions& f) { //, const SeqLib::GRC& cov_mask) {
  
  // do the overlaps
  std::vector<int32_t> q, s;
  SeqLib::GRC overlaps;
  if (f.size() > size())
    overlaps = FindOverlaps<FracRegion>(f, q, s, true);
  else 
    overlaps = f.FindOverlaps(*this, s, q, true);    

  assert(s.size() == q.size() && s.size() == overlaps.size());

  // zero all of the tiles first
  for (auto& i : *this)
    i.AddCovariate(name, 0);

  // query is fish tiles
  for (size_t i = 0; i < s.size(); ++i) {
    
    const FishHookInterval * fishr = &at(q[i]);
    const FracRegion * fracr = &f.at(s[i]);

    const double fishw = at(q[i]).Width();
    const double ovlpw = overlaps[i].Width();
    double val = (fracr->frac == DUMMY_FLOAT) ? 1 : fracr->frac; // if this is dummy val, then no score

    // Num bases in track&fish, divided by fish width, both sacled by mask coverage
    val = val * ovlpw / fishw / fishr->covered; 

    // set the value for this track
    (*this)[q[i]].AddCovariate(name, val);
  }
  
}

size_t FishHookTiles::NumCovariates() const {
  
  if (!size())
    return 0;
  
  return at(0).NumCovariates();

}

void FishHookTiles::CountEvents(const EventList& events) {

  // do the overlaps
  std::vector<int32_t> q, s;
  SeqLib::GRC overlaps;
  if (events.size() > size())
    overlaps = FindOverlaps(events, q, s, true);
  else 
    overlaps = events.FindOverlaps(*this, s, q, true);    
  
  // only one patient per bin
  std::unordered_map<int32_t, std::unordered_set<std::string>> event_donor_dedupe;

  for (size_t i = 0; i < overlaps.size(); ++i) {
    const Event * e = &events.at(s[i]);
    // if id seen before for this fish tile, skip
    if (!event_donor_dedupe[q[i]].count(e->id)) { 
      (*this)[q[i]].events++;
      event_donor_dedupe[q[i]].insert(e->id);
    }
  }
  
}
