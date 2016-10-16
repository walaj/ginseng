#ifndef GINSENG_EVENTS_H
#define GINSENG_EVENTS_H

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamHeader.h"
#include <string>

class Event : public SeqLib::GenomicRegion {

 public:
  
  Event() {}
  
  Event(const SeqLib::GenomicRegion& gr, const std::string& i) : SeqLib::GenomicRegion(gr), id(i) {}

  Event(const std::string& c, const std::string& p1, const std::string& p2, const SeqLib::BamHeader& h, const std::string& f);

  std::string id;
};

class EventList : public SeqLib::GenomicRegionCollection<Event> {

 public:
  
  // make an empty one
 EventList() : SeqLib::GenomicRegionCollection<Event>() {}
  
 bool readFromBed(const std::string& file, const SeqLib::BamHeader& h);

 private:

};

#endif
