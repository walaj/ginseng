#include "Events.h"
#include <zlib.h>

Event::Event(const std::string& c, const std::string& p1, const std::string& p2, const SeqLib::BamHeader& h, const std::string& f) : SeqLib::GenomicRegion(c, p1, p2, h) {
  // if no id, give random id
  id = f.empty() ? std::to_string(rand()) : f;
}

bool EventList::readFromBed(const std::string& file, const SeqLib::BamHeader& h) {
  
  gzFile fp = NULL;
  fp = strcmp(file.c_str(), "-")? gzopen(file.c_str(), "r") : gzdopen(fileno(stdin), "r");
  
  if (file.empty() || !fp) {
    std::cerr << "BED file not readable: " << file << std::endl;
    return false;
  }
  
  // http://www.lemoda.net/c/gzfile-read/
  while (1) {

    int err;                    
    char buffer[GZBUFFER];
    gzgets(fp, buffer, GZBUFFER);
    int bytes_read = strlen(buffer);

    // get one line
    if (bytes_read < GZBUFFER - 1) {
      if (gzeof (fp)) break;
      else {
	const char * error_string;
	error_string = gzerror (fp, &err);
	if (err) {
	  fprintf (stderr, "Error: %s.\n", error_string);
	  exit (EXIT_FAILURE);
	}
      }
    }

    // prepare to loop through each field of BED line
    size_t counter = 0;
    std::string chr, pos1, pos2, id;
    std::string line(buffer);
    std::istringstream iss_line(line);
    std::string val;
    if (line.empty() || line.at(0) == '#')
      continue;
    
    // loop VCF lines
    while(std::getline(iss_line, val, '\t')) {
      switch (counter) { 
      case 0 : chr = val; break; 
      case 1 : pos1 = val; break;
      case 2 : pos2 = val; break;
      case 3 : id = val; break;
      }
      if (counter >= 3)
	break;
      ++counter;
    }
    
    // construct the Event
    Event gr(chr, pos1, pos2, h, id);
    if (gr.chr >= 0) 
      m_grv->push_back(gr);
  }

  return true;

}

