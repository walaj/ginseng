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

// should we check length 
#define LENGTH_CHECK 1

// define a mask so we can store two chr in one uint16_t
#define CHR1_MASK = 0xFF;
#define CHR2_MASK = 0xFF00;

#define MAX_RAR_SIZE 200e6

// dont include VCFs files with more than this many non-comment lines
#define PER_FILE_LIMIT 5000

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

void Matrix::FillQuantileHistograms(size_t num_bins) {

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
  if (m_anim_step > 0 && id == 0) {
    std::string anim_file = "animation.csv";
    std::string anim_hist_file = "animation.histogram.csv";
    std::string anim_hist_small_file = "animation.histogram.small.csv";
    of_anim.open(anim_file.c_str());
    of_anim_hist.open(anim_hist_file.c_str());
    of_anim_hist_small.open(anim_hist_small_file.c_str());
  }
#endif
  
  // hash the original matrix, for later comparison
  for (auto& i : m_vec) 
    for (auto& j : i)
      m_orig_map[j.c->PointString() + j.r->PointString()] = true;
  
  // make the random values
  generateRandomVals();

  // do the swaps
  for (S i = 0; i < m_num_steps; i++) {
    size_t chr = rand_chr[m_mcmc.swap_tried];
    if (chr == INTER)
      doTransSwap();
    else
      doIntraSwap(chr);
    ++m_mcmc.swap_tried;
  }

  // print it out if need be
  std::cerr << "Swapped " << id << " matrices -- " << " shared " << ((double)shared()/(double)(m_intra+m_inter)) << " ";
  std::cerr << printMCMC() << std::endl;
  
#ifdef ANIMATION
  if (m_anim_step > 0 && id == 0) {
    of_anim.close();
    of_anim_hist.close();
    of_anim_hist_small.close();
  }
#endif

  //pthread_mutex_unlock(lock);
  
}

void Matrix::doTransSwap() {

  size_t i1 = rand_rows[m_mcmc.swap_tried] % m_vec[INTER].size();
  size_t i2 = rand_cols[m_mcmc.swap_tried] % m_vec[INTER].size();
 
  const MatrixValue * mo1 = &m_vec[INTER][i1];
  const MatrixValue * mo2 = &m_vec[INTER][i2];

  // inter switching to intra not valid
  if (mo1->c->chr != mo2->r->chr || mo1->r->chr != mo2->c->chr)
    return;

  // make the swapped vals
  MatrixValue ms1, ms2;
  MatrixValue::swapIntra(mo1, mo2, ms1, ms2);

  // add the new ones
  m_vec[INTER][i1] = ms1;
  m_vec[INTER][i2] = ms2;

}

void Matrix::doIntraSwap(size_t chr) {

#ifdef ANIMATION
  // output the animation
  if (m_anim_step > 0 && id == 0)
    if ( (m_mcmc.swap_tried % m_anim_step) == 0 || m_mcmc.swap_tried == 0) {
      std::cerr << "...writing animation for step " << m_mcmc.swap_tried << " matrix id " << id << std::endl;
      this->toCSV(of_anim, of_anim_hist, of_anim_hist_small, m_mcmc.swap_tried);
    }
#endif

  size_t i1 = rand_rows[m_mcmc.swap_tried] % m_vec[chr].size();
  size_t i2 = rand_cols[m_mcmc.swap_tried] % m_vec[chr].size();
 
  const MatrixValue * mo1 = &m_vec[chr][i1];
  const MatrixValue * mo2 = &m_vec[chr][i2];

  // make the swapped vals
  MatrixValue ms1, ms2;
  MatrixValue::swapIntra(mo1, mo2, ms1, ms2);

  // faster distance calc
  S d3 = std::abs(ms1.r->pos1 - ms1.c->pos1); //ms1.distance();
  S d4 = std::abs(ms2.r->pos1 - ms2.c->pos1); //ms1.distance();

#ifdef LENGTH_CHECK
  // length valid for intra-chrom?
  if (d3 < m_min_size || d4 < m_min_size ||  d3 > m_max_size || d4 > m_max_size)
    return;
#endif

  S d1 = std::abs(mo1->r->pos1 - mo1->c->pos1); //mo1->distance();
  S d2 = std::abs(mo2->r->pos1 - mo2->c->pos1); //mo2->distance();

  // positive energy shift is move AWAY from optimal
  int es = energyShift(d1, d2, d3, d4); 

  bool valid = (es <= 0 || probs[0][m_mcmc.swap_tried] == TRAND) ? true: false;
  if (!valid && probs[0][m_mcmc.swap_tried] > 10) { // must be more than 10/65535 to even try
    assert(es > 0 && es <= 4);
    size_t randt = rand_cols[m_mcmc.swap_tried] & TRAND;
    // proceed if probability greater than random tval number
    if (probs[es-1][m_mcmc.swap_tried] > randt)
      valid = true;
  }

  if (valid) {
    
    ++m_mcmc.accepted;
    
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
    
    // add the new ones
    m_vec[chr][i1] = ms1;
    m_vec[chr][i2] = ms2;
  }
  
}


void Matrix::toSimpleCSV(const std::string& file) {

  std::ofstream fs;
  fs.open(file);
  
  for (const auto& i : m_vec) 
    for (auto& j : i)
      fs << j.r->chr << "," << j.r->pos1 << "," << j.c->chr << "," << j.c->pos1 << std::endl;
}

void Matrix::toCSV(std::ofstream &fs, std::ofstream &fh, std::ofstream &fh_small, size_t step) {
  
  // get the number of shared object
  //int shared_count = -1;
  //if (m_orig)
  //  shared_count = shared();

  for (const auto& i : m_vec) 
    for (const auto& j : i)
      fs << j.r->chr << "," << j.r->pos1 << "," << j.c->chr << "," << j.c->pos1 << "," << step << "," << std::endl;

  m_mcmc.old_accepted = m_mcmc.accepted;
  
  // print out the histogram
  hist_swap.toCSV(fh);

  // print out the small histogram
  m_hist_smallbins.toCSV(fh_small);

}

void Matrix::writeGzip(ogzstream * out) const {

  char sep = '\t';
  for (auto& i : m_vec) {
    for (const auto& j : i) 
      (*out) << j.r->chr << sep << j.r->pos1 << sep 
			 << j.c->chr << sep << j.c->pos1 << sep
			 << id << std::endl;
  }
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

bool Matrix::LoadBEDPE(const std::string& file) {
  
  std::string val, event_line;
  igzstream this_file(file.c_str());
  
  if (!this_file)
    return false;

  // loop through the lines of thie BEDPE
  while (std::getline(this_file, event_line)) {

    // skip headers and weird chr
    if (!ValidateLine(event_line)) 
      continue;

    std::string chr1, pos1, chr2, pos2;
    std::istringstream iss(event_line);// for holdoing each line of the BEDPE
    
    size_t count = 0;
    while (std::getline(iss, val, '\t')) {
      switch(++count) {
      case 1: chr1 = val; break;
      case 2: pos1 = val; break;
      case 4: chr2 = val; break;
      case 5: pos2 = val; continue; //break;
      } // end switch
    } // end intra-line loop
    
    if (!AddNewMatrixValue(chr1, chr2, pos1, pos2)) {
      std::cerr << "Error converting line: " << event_line 
		<< " in file " << file << std::endl;
      return false;
    }
    
  } // end while read loop
  
  return true;

}

bool Matrix::LoadVCF(const std::string& file) {
  
  std::string val, event_line;
  igzstream this_file(file.c_str());
  
  if (!this_file)
    return false;

  // loop through the lines of thie BEDPE
  while (std::getline(this_file, event_line)) {

    // skip headers and weird chr
    if (!ValidateLine(event_line)) 
      continue;

    std::string chr1, pos1, chr2, pos2;
    std::istringstream iss(event_line);// for holdoing each line of the BEDPE
    
    size_t count = 0;
    while (std::getline(iss, val, '\t')) {
      switch(++count) {
      case 1: chr1 = val; break;
      case 2: pos1 = val; break;
      case 4: chr2 = val; break;
      case 5: 
	std::regex reg(".*?(\\]|\\[)([0-9A-Z]+):([0-9]+).*?$");	  
	std::smatch match;
	if (std::regex_search(val, match, reg) ) {
	  chr2 = match[2];
	  pos2 = match[3];
	} else {
	  std::cerr << "No match on line "  << val << std::endl;
	}
      } // end switch
    } // end intra-line loop
    
    if (!AddNewMatrixValue(chr1, chr2, pos1, pos2)) {
      std::cerr << "Error converting line: " << event_line 
		<< " in file " << file << std::endl;
      return false;
    }
    
  } // end while read loop
  
  return true;
}

void MatrixValue::swapIntra(const MatrixValue * m1, const MatrixValue * m2, MatrixValue &n1, MatrixValue &n2) {

  n1.r = m1->r; // switch the pointers
  n1.c = m2->c;
  n2.r = m2->r;
  n2.c = m1->c;

}

MatrixValue::MatrixValue(int chr1, int pos1, int chr2, int pos2) {
  r = std::shared_ptr<MatrixPoint>(new MatrixPoint(chr1, pos1, pos1));
  c = std::shared_ptr<MatrixPoint>(new MatrixPoint(chr2, pos2, pos2));
}

MatrixValue::MatrixValue(const std::string &chr1, const std::string &pos1, const std::string &chr2, const std::string &pos2) {
  r = std::shared_ptr<MatrixPoint>(new MatrixPoint(chr1, pos1, pos1, SeqLib::BamHeader()));
  c = std::shared_ptr<MatrixPoint>(new MatrixPoint(chr2, pos2, pos2, SeqLib::BamHeader()));
}

Matrix::Matrix() {
  __initialize_mvec();
}

void Matrix::__initialize_mvec() {
  
  assert(m_vec.size() == 0);
  for (size_t i = 0; i < 26; i++)
    m_vec.push_back({});

}

void Matrix::addMatrixValue(const MatrixValue &mv) {
  
  MatrixValue mr = mv;

  // put in order with row smaller
  if (*mv.c < *mv.r) 
    mr.Flip();
    
  // intra chr
  if (mr.r->chr == mr.c->chr) {
    ++m_intra;
    m_vec[mr.r->chr].push_back(mr);

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
  //rand_rows = std::shared_ptr<(uint32_t>(malloc(m_num_steps * sizeof(uint32_t)));
  //rand_cols = std::shared_ptr<uint32_t>(malloc(m_num_steps * sizeof(uint32_t)));
  rand_rows = (uint32_t*) malloc(m_num_steps * sizeof(uint32_t));
  rand_cols = (uint32_t*) malloc(m_num_steps * sizeof(uint32_t));

  //rand_tval = (uint16_t*) malloc(m_num_steps * sizeof(uint32_t));

#ifdef SITMO_RNG
  // sitmo
  sitmo::prng_engine eng1;
  sitmo::prng_engine eng2;
  eng1.seed((id+1)*1000);
  eng2.seed((id+1)*10000);
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

size_t Matrix::shared() const {

  size_t count = 0;
  
  // run through swapped matrix and query hash from original
  for (auto& i : m_vec)
    for (auto& j : i)
      count += m_orig_map.count(j.c->PointString() + j.r->PointString());

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

/*void Matrix::dedupe() {

  size_t intra = 0;
  size_t inter = 0;
  std::vector<MVec> new_vec;

  for (auto& i : m_vec) {

    std::sort(i.begin(), i.end());
    new_vec.push_back({});

    for (size_t j = 0; j < (i.size() - 1); j++)
      if (!(i[j]->r == i[j+1].r) && !(i[j].c == i[j+1].c)) {
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
*/

OverlapResult Matrix::checkIntraUnitOverlaps(SeqLib::GRC * grvA) {

  if (!grvA->size() || m_inter_only)
    return OverlapResult(0,m_intra+m_inter);

  size_t overlap = 0;
  size_t no_overlap = 0;

  // faster way to do it
  for (auto& i : m_vec) {
    for (auto& j : i) {
      if (grvA->OverlapSameInterval(*j.r, *j.c))
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
    for (auto& i : m_vec) {
      for (auto& j : i) {
	grc1.add(GR(*j.r, ii));
	grc2.add(GR(*j.c, ii));
	assert(!m_inter_only || (j.r->chr != j.c->chr));
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
      assert(!m_inter_only); // impossible for intra-unit overlaps if inter-only
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

      if (grvA->CountOverlaps(*j.r) && grvB->CountOverlaps(*j.c)) {
	++overlap;
	if (grvA->size() == grvB->size() && grvA->OverlapSameInterval(*j.r, *j.c)) { // take it back if same bin for same event
	  --overlap;
	  ++no_overlap;
	}
	continue;
      }
      
      if (grvA->CountOverlaps(*j.c) && grvB->CountOverlaps(*j.r)) {
	++overlap;
	if (grvA->size() == grvB->size() && grvA->OverlapSameInterval(*j.r, *j.c)) { // take it back if same bin for same event
	  ++no_overlap;
	  --overlap;
	}
	continue;
      }
      ++no_overlap;
    }

  return OverlapResult(overlap, no_overlap);
 }

bool Matrix::ValidateLine(const std::string& l) const {

  // skip header
  if (l.empty() || l.at(0) == '#')
    return false;
  
  // remove non human and weird
  if (l.find("GL00") != std::string::npos || l.find("gi|") != std::string::npos || 
      l.find("chrom1") != std::string::npos)
    return false;
  
  return true;

}

bool Matrix::AddNewMatrixValue(const std::string& c1, const std::string& c2,
			       const std::string& p1, const std::string& p2) {

  // check non empty
  if (c1.empty() || p1.empty() || c2.empty() || p2.empty()) 
    return false;
  
  // convert to numbers
  try { 
    MatrixValue mv(c1, p1, c2, p2);
    
    // only add if within distance bounds
    if (mv.r->chr >= 0 && mv.c->chr >= 0 && 
	(mv.r->chr != mv.c->chr || (mv.distance() >= (int)m_min_size && mv.distance() <= (int)m_max_size)) && 
	(mv.r < mv.c) && (!m_inter_only || mv.r->chr != mv.c->chr) && (!m_intra_only || (m_intra_only && mv.r->chr == mv.c->chr))) {
      
      // add the event
      addMatrixValue(mv);
    } 
    
  } catch(...) {
    return false;
  }
  
  return true;
}

void Matrix::ScrambleAroundDiagonal() {

  // scramble to make even around the diagonal
  for (auto& i : m_vec) 
    for (auto& j : i) 
      if (rand() % 2) {
	j.Flip();
	//int idd = j.id;
	//j = MatrixValue(j.c, j.r);
	//j.id = idd;
      }
}

size_t Matrix::size() const {

  size_t c = 0;
  for (auto& i : m_vec)
    c += i.size();
  return c;
}
