#include "Matrix.h"
#include "swap.h"

#include <regex>
#include <sstream>
#include <random>
#include <cassert>
#include <fstream>
#include <algorithm>

#include "gzstream.h"
#include "SeqLib/SeqLibUtils.h"

#include <prng_engine.hpp>
#include <set>

// define a mask so we can store two chr in one uint16_t
#define CHR1_MASK = 0xFF;
#define CHR2_MASK = 0xFF00;

#define MAX_RAR_SIZE 200e6

// min distance to be valid. Gets rid of Sanger high FP rate at small events
#define SANGER_DIST_LIM 1e6
// dont include VCFs files with more than this many non-comment lines
#define SANGER_PER_FILE_LIMIT 5000

#define SITMO_RNG 1

#define MIN_BIN_WIDTH 1000

#define ANIMATION 1

#define MIN_RAR_SIZE 5000

static const size_t INTER = 25;

int __countLines(const std::string& file) {

  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(file);

  while (std::getline(myfile, line)) {
    if (line.find("#") == std::string::npos && line.find("vcf") != std::string::npos)
      ++number_of_lines;
  }
  
  return number_of_lines;

}

int drawFromPower(double x0, double x1, double power) {

  assert(power != -1);

  const int PRECISION = 1e6;

  double r = (double)(rand() % PRECISION)/(double)PRECISION;
  double t1 = std::pow(x1, power+1);
  double t2 = std::pow(x0, power+1);
  double tsum = (t1-t2) * r + t2;

  return std::floor(std::pow(tsum, 1 / (power + 1)));
}


std::vector<int> drawFromPower(double x0, double x1, double power, int n_draws) {

  assert(power != -1);

  const int PRECISION = 1e6;

  std::vector<int> rpower(n_draws, 0);
  
  for (int i = 0; i < n_draws; ++i) {
    double r = (double)(rand() % PRECISION)/(double)PRECISION;
    double t1 = std::pow(x1, power+1);
    double t2 = std::pow(x0, power+1);
    double tsum = (t1-t2) * r + t2;
    rpower[i] = std::floor(std::pow(tsum, 1 / (power + 1)));
  }

  return rpower;
}

void Matrix::add() { //const Matrix& m) {

  std::set<std::string> hash;

  // hash the existing matrix
  for (auto& i : m_vec) 
    for (auto& j : i)
      hash.insert(std::to_string(j.c.chr) + ":" + std::to_string(j.c.pos1) + 
		  std::to_string(j.r.chr) + ":" + std::to_string(j.r.pos1));

  // add the element if not there, otherwise update counter
  for (auto& i : m_vec) {
    for (auto& j : i) {
      std::string thash = std::to_string(j.c.chr) + ":" + std::to_string(j.c.pos1) + 
	std::to_string(j.r.chr) + ":" + std::to_string(j.r.pos1);
      if (hash.count(thash))
	++j.count;
      else 
	addMatrixValue(j);
    }
  }
}
  
void Matrix::fillQuantileHistograms(size_t num_bins) {

  std::vector<S>* pspanv = new std::vector<S>();

  // fill a vector of the spans
  getSpans(pspanv);

  // fill the quantile histogram
  hist.initialSpans(num_bins, pspanv, MIN_BIN_WIDTH);
  m_hist_smallbins.initialSpans(num_bins*10, pspanv, 0);
  //hist.toCSV("hist.csv");

  delete pspanv;

  // copy initial histogram to swap histogram
  hist_swap = hist;  

}

void Matrix::getSpans(std::vector<S>* pspanv) {

  pspanv->clear();
  for (auto& i : m_vec) 
    for (auto& j : i)
      pspanv->push_back(j.distance());

}

void Matrix::allSwaps() { //pthread_mutex_t * lock, std::vector<Matrix*> * allm) {
  
  if (!(m_intra+m_inter))
    return;

#ifdef ANIMATION
  // open the animation file if need to
  if (m_anim_step > 0 && m_id == 0) {
    std::string anim_file = analysis_id + ".animation.csv";
    std::string anim_hist_file = analysis_id + ".animation.histogram.csv";
    std::string anim_hist_small_file = analysis_id + ".animation.histogram.small.csv";
    of_anim.open(anim_file.c_str());
    of_anim_hist.open(anim_hist_file.c_str());
    of_anim_hist_small.open(anim_hist_small_file.c_str());
  }
#endif
  
  // hash the original matrix, for later comparison
  for (auto& i : m_vec) 
    for (auto& j : i)
      m_orig_map[j.c.PointString() + j.r.PointString()] = true;
  
  // make the random values
  generateRandomVals();

  // do the swaps
  for (S i = 0; i < m_num_steps; i++) 
    doSwap();

  // print it out if need be
  std::cerr << "Swapped " << m_id << " matrices -- " << " shared " << ((double)shared()/(double)(m_intra+m_inter)) << " ";
  std::cerr << printMCMC() << std::endl;
  
#ifdef ANIMATION
  if (m_anim_step > 0 && m_id == 0) {
    of_anim.close();
    of_anim_hist.close();
    of_anim_hist_small.close();
  }
#endif

  //pthread_mutex_unlock(lock);
  
}

void Matrix::doSwap() {

#ifdef ANIMATION
  // output the animation
  if (m_anim_step > 0 && m_id == 0)
    if ( (m_mcmc.swap_tried % m_anim_step) == 0 || m_mcmc.swap_tried == 0) {
      std::cerr << "...writing animation for step " << m_mcmc.swap_tried << " matrix id " << id << std::endl;
      this->toCSV(of_anim, of_anim_hist, of_anim_hist_small, m_mcmc.swap_tried);
    }
#endif


  // new and faster
  size_t chr = INTER; 

  if (m_intra && !inter_only)
    chr = rand_chr[m_mcmc.swap_tried];
  size_t i1 = rand_rows[m_mcmc.swap_tried] % m_vec[chr].size();
  size_t i2 = rand_cols[m_mcmc.swap_tried] % m_vec[chr].size();
 
  MatrixValue mo1 = m_vec[chr][i1];
  MatrixValue mo2 = m_vec[chr][i2];

  // make the swapped vals
  MatrixValue ms1, ms2;
  MatrixValue::swapIntra(mo1, mo2, ms1, ms2);

  S d3 = ms1.distance();
  S d4 = ms2.distance();

  // inter switching to intra. NOT VALUD
  if (chr == INTER && (d3 != INTERCHR || d4 != INTERCHR)) {
    ++m_mcmc.swap_tried;
    return;
  }

  bool valid = chr == INTER || (d3 >= min_size && d4 >= min_size && d3 <= max_size && d4 <= max_size);  
  if (!valid) {
    ++m_mcmc.swap_tried;
    return;
  }

  S d1 = mo1.distance();
  S d2 = mo2.distance();

  // positive energy shift is move AWAY from optimal
  int es = 0;
  if (chr != INTER) 
    es = energyShift(d1, d2, d3, d4); 

  valid = (es <= 0 || probs[0][m_mcmc.swap_tried] == TRAND) ? true: false;
  if (!valid && probs[0][m_mcmc.swap_tried] > 10) { // must be more than 10/65535 to even try
    assert(es > 0 && es <= 4);
    size_t randt = rand_cols[m_mcmc.swap_tried] & TRAND;
    // proceed if probability greater than random tval number
    if (probs[es-1][m_mcmc.swap_tried] > randt)
      valid = true;
  }

  ++m_mcmc.swap_tried; 
  
  if (valid) {
    
    //assert(temps[m_mcmc.swap_tried-1] >= 0.001 || es <= 0);
    ++m_mcmc.accepted;
    
    // update the swap histogram
    if (chr != INTER) {

      ++hist_swap.m_bins[bin_table[d3]];
      ++hist_swap.m_bins[bin_table[d4]];
      --hist_swap.m_bins[bin_table[d1]];
      --hist_swap.m_bins[bin_table[d2]];
      
#ifdef ANIMATION
      if (id == 0 && m_anim_step > 0) {
	m_hist_smallbins.addElem(d3);
	m_hist_smallbins.addElem(d4);
	m_hist_smallbins.removeElem(d1);
	m_hist_smallbins.removeElem(d2);
      }
#endif
    }

    // add the new ones
    m_vec[chr][i1] = ms1;
    m_vec[chr][i2] = ms2;
  }
  
}


void Matrix::toSimpleCSV(const std::string& file) {

  std::ofstream fs;
  fs.open(file);
  
  for (auto& i : m_vec) 
    for (auto& j : i)
      fs << j.r.chr << "," << j.r.pos1 << "," << j.c.chr << "," << j.c.pos1 << "," << j.count << std::endl;
}

void Matrix::toCSV(std::ofstream &fs, std::ofstream &fh, std::ofstream &fh_small, size_t step) {
  
  // get the number of shared object
  //int shared_count = -1;
  //if (m_orig)
  //  shared_count = shared();

  for (auto& i : m_vec) 
    for (auto& j : i)
      fs << j.r.chr << "," << j.r.pos1 << "," << j.c.chr << "," << j.c.pos1 << "," << j.count << "," << step << "," << j.id << std::endl;

  m_mcmc.old_accepted = m_mcmc.accepted;
  
  // print out the histogram
  hist_swap.toCSV(fh);

  // print out the small histogram
  m_hist_smallbins.toCSV(fh_small);

}

void Matrix::writeGzip(ogzstream * out) const {

  char sep = '\t';
  for (auto& i : m_vec) {
    for (auto& j : i) 
      (*out) << j.r.chr << sep << j.r.pos1 << sep 
			 << j.c.chr << sep << j.c.pos1 << sep
			 << m_id << std::endl;
  }
}

void Matrix::writeBinary() const {
  /*
  FILE* binout;
  std::string name = "matrix" + std::to_string(m_id) + ".bin";
  binout = fopen(name.c_str(),"wb");
  for (auto& i : m_vec) {
    for (auto& j : i) {
      uint32_t buffer[] = {j.r.chr, j.r.pos1, j.c.chr, j.c.pos1};
      fwrite(buffer, sizeof(uint32_t), sizeof(buffer), binout);
    }
  }
  fclose(binout);
  */
}

int Matrix::energyShift(S odist1, S odist2, S pdist1, S pdist2) {

  // negative energy means more favorable histogram is being made

  // calculate energy shift
  //int shift = hist.binCount(bin_table[odist1]) - hist_swap.binCount(bin_table[odist1]) +
  //  hist.binCount(bin_table[odist2]) - hist_swap.binCount(bin_table[odist2]) + 
  //  hist_swap.binCount(bin_table[pdist1]) - hist.binCount(bin_table[pdist1]) + 
  //  hist_swap.binCount(bin_table[pdist2]) - hist.binCount(bin_table[pdist2]);

  int shift = 0;
  int weight1 = abs(hist_swap.binCount(bin_table[odist1]) - hist.binCount(bin_table[odist1]));
  shift += (hist_swap.binCount(bin_table[odist1]) >  hist.binCount(bin_table[odist1]) ? -1 : 1)*weight1;
  int weight2 = abs(hist_swap.binCount(bin_table[odist2]) - hist.binCount(bin_table[odist2]));
  shift += (hist_swap.binCount(bin_table[odist2]) >  hist.binCount(bin_table[odist2]) ? -1 : 1)*weight2;
  int weight3 = abs(hist_swap.binCount(bin_table[pdist1]) - hist.binCount(bin_table[pdist1]));
  shift += (hist_swap.binCount(bin_table[pdist1]) >= hist.binCount(bin_table[pdist1]) ? 1 : -1)*weight3;
  int weight4 = abs(hist_swap.binCount(bin_table[pdist2]) - hist.binCount(bin_table[pdist2]));
  shift += (hist_swap.binCount(bin_table[pdist2]) >= hist.binCount(bin_table[pdist2]) ? 1 : -1)*weight4;

  shift = shift > 4 ? 4 : shift;
  //shift = shift < -4 ? -4 : shift;

  //if (shift < -4)
  //  shift = -4;
  //else if (shift > 4)
  //  shift = 4;
      
  return shift;

}

Matrix::Matrix(const std::string &file_list, size_t nb, size_t ns, SeqLib::GRC &mk, bool inter_only, const std::vector<std::string>& identifiers, const std::string& tid, int tmin_size, int tmax_size, SeqLib::GRC& black, bool intra_only) : m_num_bins(nb), m_num_steps(ns), analysis_id(tid), m_intra_only(intra_only) {

  min_size = tmin_size >= 0 ? tmin_size : 0;
  max_size = tmax_size >= 0 ? tmax_size : 0;

  this->inter_only = inter_only;

  size_t good_files = 0;
  
  //open the file
  igzstream inFile(file_list.c_str());
  if (!inFile) {
    std::cerr << "Can't read file " << file_list << " for parsing VCF" << std::endl;
    return;
  }

  // get count of files
  size_t num_files = __countLines(file_list);

  // counter for number of masked events
  size_t masked = 0;

  // initialize the m_vec
  __initialize_mvec();

  // loop through the files and add the events
  std::string fileline, file;
  size_t file_counter = 0;
  while (std::getline(inFile, fileline)) {
    
    // get the first field, its the file name
    std::istringstream issf(fileline);
    while(std::getline(issf, file, '\t')) 
      break;
    
    ++file_counter;
    size_t new_events = 0;
    if (!SeqLib::read_access_test(file)) {
      std::cerr << "VCF File: "  << file << " not readable or does not exist" << std::endl;
      continue;
    }

    // check the identifiers
    bool good = identifiers.size() == 0;
    for (auto& id : identifiers) {
      if (fileline.find(id) != std::string::npos) {
	good = true;
	break;
      }
    }
    if (!good)
      continue;

    ++good_files;

    std::string val;
    igzstream this_file(file.c_str());
    std::string event_line;
    
    size_t num_events_in_file = __countLines(file);
    if (num_events_in_file < SANGER_PER_FILE_LIMIT) // block if too many lines
      while (std::getline(this_file, event_line)) {
	if (event_line.length() > 0 && event_line.at(0) != '#') {
	  size_t count = 0;
	  std::string chr1 = "-1", pos1 = "-1", chr2 = "-1", pos2 = "-1";
	  std::istringstream iss(event_line);
	  
	  // remove non human
	  if (event_line.find("GL00") != std::string::npos || event_line.find("gi|") != std::string::npos)
	    continue;

	  // regex to get mate information
	  while (std::getline(iss, val, '\t')) {
	    switch(++count) {
	    case 1: chr1 = val; break;
	    case 2: pos1 = val; break;
	    case 5: 
	      std::regex reg(".*?(\\]|\\[)([0-9A-Z]+):([0-9]+).*?$");	  
	      std::smatch match;
	      if (std::regex_search(val, match, reg) ) {
		chr2 = match[2];
		pos2 = match[3];
	      } else {
		std::cerr << "No match on line "  << val << std::endl;
	      }
	      break;
	    } // end switch
	  } // end intra-line loop
	  
	  // check that we have the data
	  if (chr1 == "-1" || pos1 == "-1" || chr2 == "-1" || pos2 == "-1") {
	    std::cerr << "Failed to parse VCF on line " << event_line << std::endl;
	    continue;
	  }
	  
	  // convert to numbers
	  try { 
	    
	    MatrixValue mv(chr1, pos1, chr2, pos2);
	    mv.id = good_files;// give unique id
	    
	    // check the mask
	    bool keep = true;
	    if (mk.size()) {
	      size_t ovl = 0;
	      ovl += mk.CountOverlaps(GenomicRegion(mv.r));	      
	      if (!ovl) 
		ovl += mk.CountOverlaps(mv.c);
	      if (ovl) {
		keep = false;
		++masked;
	      }
	    }
	    
	    // sanger conditional on distance
	    bool blacklist_pass = !black.CountOverlaps(mv.r) && !black.CountOverlaps(mv.c);
	    if (blacklist_pass && mv.r.chr >= 0 && mv.c.chr >= 0 && (mv.r.chr != mv.c.chr || (mv.distance() >= (int)min_size && mv.distance() <= (int)max_size)) && keep && (mv.r < mv.c) && (!inter_only || mv.r.chr != mv.c.chr) && (!m_intra_only || (m_intra_only && mv.r.chr == mv.c.chr))) {
	      
	      // add the event
	      addMatrixValue(mv);
	      
	      // increment the event counter
	      ++new_events;
	    } 

	  } catch(...) {
	    std::cerr << "********************************" << std::endl;
	    std::cerr << "********************************" << std::endl;
	    std::cerr << "Error converting VCF line from std::string to number. Line " << event_line << std::endl;
	    std::cerr << "chr1 " << chr1 << " chr2 " << chr2 << std::endl;
	    std::cerr << "********************************" << std::endl;
	    std::cerr << "********************************" << std::endl;
	  }
	} // end ## conditional
      } // end intra-file loop
    else 
      std::cerr << "!!! Removed file " << file << " with too many events: " << num_events_in_file<< std::endl;
    
      // remove events if too many
      //if (new_events > MAX_EVENTS && false) {
      //	std::cerr << "!!! Too many events: " << new_events << " -- ignoring all in this file !!!" << std::endl;
      //	m.erase(m.end() - new_events, m.end());
      //}
      
    if (file_counter % 100 == 0)
      std::cerr << "...imported VCF file " << file_counter << " of " << num_files << " with " << new_events << " events " << std::endl;
    
  } // end VCF file list
  
  if (m_verbose) 
    std::cerr << "Imported " << (m_inter+m_intra) << " breakpoints from " << good_files << " file. Masked out " << masked << " events " << std::endl;

  if ((m_intra+m_inter) == 0) {
    std::cerr << "********************************" << std::endl;
    std::cerr << "********************************" << std::endl;
    std::cerr << "ERROR: Didn't import any events" << std::endl;
    std::cerr << "********************************" << std::endl;
    std::cerr << "********************************" << std::endl;
    exit(EXIT_FAILURE);
  }

  //dedupe();

  // scramble to make even around the diagonal
  for (auto& i : m_vec) 
    for (auto& j : i) 
      if (rand() % 2) {
	int idd = j.id;
	j = MatrixValue(j.c, j.r);
	j.id = idd;
      }

  // fill event histograms
  fillQuantileHistograms(m_num_bins);

  m_orig = this;

}

void MatrixValue::swapIntra(const MatrixValue &m1, const MatrixValue &m2, MatrixValue &n1, MatrixValue &n2) {

  n1.r = m1.r;
  n1.c = m2.c;
  n2.r = m2.r;
  n2.c = m1.c;

}

void MatrixValue::swapInter(const MatrixValue &m1, const MatrixValue &m2, MatrixValue &n1, MatrixValue &n2) {

  n1.r = m1.r;
  n1.c = m2.c;
  n2.r = m2.r;
  n2.c = m1.c;

}

MatrixValue::MatrixValue(int chr1, int pos1, int chr2, int pos2) {
  r = GenomicRegion(chr1, pos1, pos1);
  c = GenomicRegion(chr2, pos2, pos2);
}

MatrixValue::MatrixValue(const std::string &chr1, const std::string &pos1, const std::string &chr2, const std::string &pos2) {
  
  r = GenomicRegion(chr1, pos1, pos1, SeqLib::BamHeader());
  c = GenomicRegion(chr2, pos2, pos2, SeqLib::BamHeader());
  
}

void Matrix::addMatrixValue(const MatrixValue &mv) {
  
  MatrixValue mr = mv;

  // put in order with row smaller
  if (mv.c < mv.r) {
    mr.c = mv.r;
    mr.r = mv.c;
  }
    
  // intra chr
  if (mv.r.chr == mv.c.chr) {
    ++m_intra;
    m_vec[mr.r.chr].push_back(mr);

  // inter chr
  } else {
    ++m_inter;
    m_vec[INTER].push_back(mr);
    assert(m_vec[INTER].size() == m_inter);
  }
}

void Matrix::generateRandomVals() {

  //clock_t t = clock();

  // destructor is responsible for freeing
  rand_rows = (uint32_t*) malloc(m_num_steps * sizeof(uint32_t));
  rand_cols = (uint32_t*) malloc(m_num_steps * sizeof(uint32_t));
  //rand_tval = (uint16_t*) malloc(m_num_steps * sizeof(uint32_t));

#ifdef SITMO_RNG
  // sitmo
  sitmo::prng_engine eng1;
  sitmo::prng_engine eng2;
  eng1.seed((m_id+1)*1000);
  eng2.seed((m_id+1)*10000);
  size_t sumr = m_intra+m_inter;
  for (size_t i = 0; i < m_num_steps; ++i){
    //size_t r1 = eng1();
    rand_rows[i] = eng1() % sumr;
    rand_cols[i] = eng2() % sumr;
    //rand_tval[i] = r1 % TRAND;
  }
#else

  // fill the values
  unsigned short seed[3] = {155,0,155};
  for (size_t i = 0; i < m_num_steps; ++i){
    rand_rows[i] = nrand48(&seed[0]) % (m_intra+m_inter);
    rand_cols[i] = nrand48(&seed[0]) % (m_intra+m_inter);
    //rand_tval[i] = rand_rows[i] % TRAND;
  }
#endif

}

size_t Matrix::shared() {

  size_t count = 0;
  
  // run through swapped matrix and query hash from original
  for (auto& i : m_vec)
    for (auto& j : i)
      count += m_orig_map.count(j.c.PointString() + j.r.PointString());

  return count;

}

std::ostream& operator<<(std::ostream &out, const MCMCData &m) 
{
    out << "   Swap tried: " << m.swap_tried << " accepted " << m.accepted << " (" << SeqLib::percentCalc<size_t>(m.accepted, m.swap_tried) << "%)";
    return out;
}

std::string Matrix::printMCMC() const {
  std::stringstream ss;
  ss << m_mcmc;
  return ss.str();
}

void Matrix::dedupe() {

  size_t intra = 0;
  size_t inter = 0;
  std::vector<MVec> new_vec;

  for (auto& i : m_vec) {

    std::sort(i.begin(), i.end());
    new_vec.push_back({});

    for (size_t j = 0; j < (i.size() - 1); j++)
      if (!(i[j].r == i[j+1].r) && !(i[j].c == i[j+1].c)) {
	new_vec.back().push_back(i[j]);
	//if (i == -1)
	//	  ++inter;
	//else
	//  ++intra;
      }

  }
  
  std::cerr << "Events before dedupe: " << (m_intra+m_inter) << " after " << (inter+intra) << std::endl;
  m_vec = new_vec;
  m_intra = intra;
  m_inter = inter;
}

OverlapResult Matrix::checkIntraUnitOverlaps(SeqLib::GRC * grvA) {

  if (!grvA->size() || inter_only)
    return OverlapResult(0,m_intra+m_inter);

  size_t overlap = 0;
  size_t no_overlap = 0;

  // faster way to do it
  for (auto& i : m_vec) {
    for (auto& j : i) {
      if (grvA->OverlapSameInterval(j.r, j.c))
	++overlap;
      else 
	++no_overlap;
    }
  }
  
  return OverlapResult(overlap, no_overlap);
  
  // make grc
  //bool ff = true;
  
  int ii = 0;
  if (!grc1.size() || !grc2.size()) {  // only make this once
    grc1.clear(); grc2.clear();
    //std::cerr << " creating grc for " << m_id << std::endl;
    for (auto& i : m_vec) {
      for (auto& j : i) {
	//if (ff){ std::cerr << j.r << " - " << j.c << std::endl; ff = false; } 
	grc1.add(GR(j.r, ii));
	grc2.add(GR(j.c, ii));
	assert(!inter_only || (j.r.chr != j.c.chr));
	++ii;
      }
    }
    grc1.CreateTreeMap(); // sorts it, so query order is messed up
    grc2.CreateTreeMap();
    
  }

  std::vector<int32_t> sub1, que1, sub2, que2;
  grvA->FindOverlaps(grc1, que1, sub1, true);
  grvA->FindOverlaps(grc2, que2, sub2, true);

  // make subject to query hash
  std::unordered_map<size_t,size_t> sub1_to_que1;
  std::unordered_map<size_t,size_t> sub2_to_que2;
  for (size_t i = 0; i < que2.size(); ++i) {
    sub2_to_que2[grc2.at(sub2[i]).id] = que2[i]; // each subject (event) should have unique query (bed region). Not vice versa
  }

  // want to count number of unique subject ids that share the same query id
  // loop through bed track regions
  for (size_t i =0; i < sub1.size(); ++i) {
    std::unordered_map<size_t, size_t>::const_iterator tt = sub2_to_que2.find(grc1.at(sub1[i]).id);
    if (tt != sub2_to_que2.end() && (int)tt->second == que1[i]) { // q1 (bed region for Row) and q2 (bed region for Column) must be equal for a given subject (event)
      ++overlap;
      assert(!inter_only); // impossible for intra-unit overlaps if inter-only
    } else {
      ++no_overlap;
    }
  }

  return OverlapResult(overlap, no_overlap);

}

OverlapResult Matrix::checkOverlaps(SeqLib::GRC * grvA, SeqLib::GRC * grvB) {

  if (!grvA->size() || !grvB->size())
    return OverlapResult(0,m_intra+m_inter);

  size_t overlap = 0;
  size_t no_overlap = 0;

  for (auto& i : m_vec)
    for (auto& j : i) {

      if (grvA->CountOverlaps(j.r) && grvB->CountOverlaps(j.c)) {
	++overlap;
	if (grvA->size() == grvB->size() && grvA->OverlapSameInterval(j.r, j.c)) { // take it back if same bin for same event
	  --overlap;
	  ++no_overlap;
	}
	continue;
      }
      
      if (grvA->CountOverlaps(j.c) && grvB->CountOverlaps(j.r)) {
	++overlap;
	if (grvA->size() == grvB->size() && grvA->OverlapSameInterval(j.r, j.c)) { // take it back if same bin for same event
	  ++no_overlap;
	  --overlap;
	}
	continue;
      }
      ++no_overlap;
    }

  return OverlapResult(overlap, no_overlap);
 }

void Matrix::__initialize_mvec() {
  
  assert(m_vec.size() == 0);
  for (size_t i = 0; i < 26; i++)
    m_vec.push_back({});

}

