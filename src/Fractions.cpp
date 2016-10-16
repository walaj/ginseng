#include "Fractions.h"

using namespace SeqLib;

FracRegion::FracRegion(const std::string& c, const std::string& p1, const std::string& p2, const SeqLib::BamHeader& h, const std::string& f) : SeqLib::GenomicRegion(c, p1, p2, h)
  {
    
    if (f.empty()) {
      frac = DUMMY_FLOAT;
      return;
    }

    // convert frac to double
    try { 
      frac = std::stod(f);
    } catch (...) {
      std::cerr << "FracRegion::FracRegion - Error converting fraction " << f << " to double "  << std::endl;
      exit(EXIT_FAILURE);
    } 
  }

std::ostream& operator<<(std::ostream& out, const FracRegion& f) {
  out << f.chr << ":" << SeqLib::AddCommas<int32_t>(f.pos1) << "-" << 
    SeqLib::AddCommas<int32_t>(f.pos2) << " Frac: " << f.frac;
  return out;
}

void Fractions::readFromBed(const std::string& file, const SeqLib::BamHeader& h, bool scored) {
  
  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "BED file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string line;
  std::string curr_chr = "-1";
  while (std::getline(iss, line, '\n')) {
    
    size_t counter = 0;
    std::string chr, pos1, pos2, f;
    std::istringstream iss_line(line);
    std::string val;
    
    if (line.find("#") == std::string::npos) {
      while(std::getline(iss_line, val, '\t')) {
	switch (counter) { 
	case 0 : chr = val; break; 
	case 1 : pos1 = val; break;
	case 2 : pos2 = val; break;
	case 3 : f = val; break;
	}
	if (counter >= 3 || (counter == 2 && !scored))
	  break;
	++counter;
	
	if (chr != curr_chr) {
	  //std::cerr << "...reading from BED - chr" << chr << std::endl;
	  curr_chr = chr;
	}
	
      }

      //remove chr
      // https://www.safaribooksonline.com/library/view/c-cookbook/0596007612/ch04s12.html
      const std::string s = "chr";
      std::string::size_type i = chr.find(s);
      if (i != std::string::npos)
	chr.erase(i, s.length());

      // construct the GenomicRegion
      FracRegion ff(chr, pos1, pos2, h, f);
      add(ff);
	
	//}
    } // end "keep" conditional
  } // end main while


}
