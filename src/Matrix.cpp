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

#include "gif.h"

// should we check length 
#define LENGTH_CHECK 1

// animation params for gif
#define WIDTH 200
#define COLSTEP 15 // make smaller to increase dynamic range
#define DELAY 20 // in 100th of seconds

// which rng to use
#define SITMO_RNG 1

#define MIN_BIN_WIDTH 1000

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

void Matrix::allSwaps() { 
  
  if (!(m_intra+m_inter)) {
    std::cerr << "NO EVENTS TO SWAP" << std::endl;
    return;
  }

#ifdef WIDTH
  // open the animation file if need to
  if (m_anim_step > 0 && id == 1) {
    std::string anim_file = "animation.csv";
    std::string anim_hist_file = "animation.histogram.csv";
    std::string anim_hist_small_file = "animation.histogram.small.csv";
    of_anim.open(anim_file.c_str());
    of_anim_hist.open(anim_hist_file.c_str());
    of_anim_hist_small.open(anim_hist_small_file.c_str());

    // open for writing
    gw = new GifWriter();
    GifBegin(gw, "matrix.gif", WIDTH, WIDTH, 1);
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
  
#ifdef WIDTH
  if (m_anim_step > 0 && id == 1) {
    of_anim.close();
    of_anim_hist.close();
    of_anim_hist_small.close();
    GifEnd(gw); // close the gif
  }
#endif

}

void Matrix::doTransSwap() {
  
#ifdef WIDTH
  if (m_anim_step > 0 && id == 1)
    Animate();
#endif

  size_t i1 = rand_rows[m_mcmc.swap_tried] % m_vec[INTER].size();
  size_t i2 = rand_cols[m_mcmc.swap_tried] % m_vec[INTER].size();
 
  const MatrixValue * mo1 = &m_vec[INTER][i1];
  const MatrixValue * mo2 = &m_vec[INTER][i2];

  // inter switching to intra not valid
  if (mo1->c->chr == mo2->r->chr || mo1->r->chr == mo2->c->chr)
    return;

  // make the swapped vals
  MatrixValue ms1, ms2;
  MatrixValue::swap(mo1, mo2, ms1, ms2);

  // add the new ones
  m_vec[INTER][i1] = ms1;
  m_vec[INTER][i2] = ms2;
}

void Matrix::doIntraSwap(size_t chr) {

#ifdef WIDTH
  if (m_anim_step > 0 && id == 1)
    Animate();
#endif

  size_t i1 = rand_rows[m_mcmc.swap_tried] % m_vec[chr].size();
  size_t i2 = rand_cols[m_mcmc.swap_tried] % m_vec[chr].size();
 
  const MatrixValue * mo1 = &m_vec[chr][i1];
  const MatrixValue * mo2 = &m_vec[chr][i2];

  // make the swapped vals
  MatrixValue ms1, ms2;
  MatrixValue::swap(mo1, mo2, ms1, ms2);

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
    
#ifdef WIDTH
    if (id == 1 && m_anim_step > 0) {
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

void MatrixValue::swap(const MatrixValue * m1, const MatrixValue * m2, MatrixValue &n1, MatrixValue &n2) {

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

Matrix::Matrix() : rand_chr(nullptr), rand_rows(nullptr), rand_cols(nullptr) {
  
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
  out << "   Swap tried: " << SeqLib::AddCommas(m.swap_tried)
      << " accepted " << SeqLib::AddCommas(m.accepted) 
      << " (" << SeqLib::percentCalc<size_t>(m.accepted, m.swap_tried) << "%)";
    return out;
}

std::string Matrix::printMCMC() const {
  std::stringstream ss;
  ss << m_mcmc;
  return ss.str();
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
    if (mv.r->chr >= 0 && mv.c->chr >= 0 && // valid chromosome
	(mv.r->chr != mv.c->chr || (mv.distance() >= (int)m_min_size && mv.distance() <= (int)m_max_size)) && // within distance bounds or translocation
	//(mv.r < mv.c) && // only add one half of the events to not duplicate. For BEDBE, it's already deduped
	(!m_inter_only || mv.r->chr != mv.c->chr) && // if inter_only, exclude non-translocations 
	(!m_intra_only || mv.r->chr == mv.c->chr)) { // if intra-only, exclude non-intras
      
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

uint8_t * Matrix::AsRGB8Image(uint32_t width, int step) const {

  // chr sizes for hg19
  const uint32_t csum[24] = {0, 249250621, 492449994, 690472424, 881626700, 1062541960,
                             1233657027, 1392795690, 1539159712, 1680373143,
                             1815907890, 1950914406, 
			     2084766301, 2199936179, 2307285719, 2409817111, 2500171864, 2581367074,
			     2659444322, 2718573305, 2781598825, 2829728720, 2881033286, 3036303846};
  
  uint8_t * out = (uint8_t*) calloc(width*width*4, sizeof(uint8_t));

  // convert to absolute coordinate
  const double maxlen = 3036303846;
  for (const auto& c : m_vec) {
    for (const auto& m : c) {
      size_t a = std::floor((double)(csum[m.c->chr] + m.c->pos1) / maxlen * width);
      size_t b = std::floor((double)(csum[m.r->chr] + m.r->pos1) / maxlen * width);
      a = a >= width ? (width - 1) : a;
      b = b >= width ? (width - 1) : b;
      size_t hit = (a * width + b) * 4;
      out[hit  ] = out[hit  ] < (256-step) ? out[hit  ] + step : 255; //r
      out[hit+1] = 0;//255;//out[hit+1] < 255 ? out[hit+1] + 1 : 255; //g
      out[hit+2] = 0;//255;//out[hit+2] < 255 ? out[hit+2] + 1 : 255; //b
      out[hit+3] = 0; //out[hit+3] =  < 255 ? out[hit+3] + 1 : 255; //a
    }
  }

  return out;
}

uint8_t * Matrix::AsBMPImage(uint32_t width, int step) const {

  // chr sizes for hg19
  const uint32_t csum[24] = {0, 249250621, 492449994, 690472424, 881626700, 1062541960,
                             1233657027, 1392795690, 1539159712, 1680373143,
                             1815907890, 1950914406, 
			     2084766301, 2199936179, 2307285719, 2409817111, 2500171864, 2581367074,
			     2659444322, 2718573305, 2781598825, 2829728720, 2881033286, 3036303846};
  
  uint8_t * out = (uint8_t*) calloc(width*width*3, sizeof(uint8_t));

  // convert to absolute coordinate
  const double maxlen = 3036303846;
  for (const auto& c : m_vec) {
    for (const auto& m : c) {
      size_t a = std::floor((double)(csum[m.c->chr] + m.c->pos1) / maxlen * width);
      size_t b = std::floor((double)(csum[m.r->chr] + m.r->pos1) / maxlen * width);
      a = a >= width ? (width - 1) : a;
      b = b >= width ? (width - 1) : b;
      size_t hit = (a * width + b) * 3;
      out[hit  ] = 0; //out[hit  ] < (256-step) ? out[hit  ] + step : 255; //b
      out[hit+1] = 0; //out[hit+1] < 255 ? out[hit+1] + 1 : 255;           //g
      out[hit+2] = out[hit+2] < (256-step) ? out[hit+2] + step : 255;      //r
    }
  }

  return out;
}

Matrix::~Matrix() { 
  if (rand_rows)
    free(rand_rows); 
  if (rand_cols)
    free(rand_cols); 
  if (gw)
      delete gw;
}

void Matrix::Animate() {

  // output the animation
  if ( (m_mcmc.swap_tried % m_anim_step) == 0 || m_mcmc.swap_tried == 0) {
    uint8_t * image = AsRGB8Image(WIDTH, COLSTEP);
    GifWriteFrame(gw, image, WIDTH, WIDTH, DELAY);
    std::cout << "animating_step\t" << m_mcmc.swap_tried 
	      << "\thist_distance\t" << hist_swap.EuclideanDistance(m_orig->hist_swap) << std::endl;
    hist_swap.toCSV(of_anim_hist); // write the histograms
    m_hist_smallbins.toCSV(of_anim_hist_small);
  }
}

void Matrix::AddBedElements(const std::string& b, BEDIntervals& bi) {

  // track which bed tracks we are tracking
  bed_list.push_back(b);

  // trackwhat is in and out of this BED track
  size_t in = 0, out = 0;

  for (auto& chr : m_vec) {
    for (auto& d : chr) {
      d.r->AddBED(b, bi); // eg SINE and grc of SINE track
      d.c->AddBED(b, bi);

      // sanity check. maybe cut later
      assert((d.r->olap_element[b] >= 0) == (bi.grc.CountOverlaps(*d.r) > 0));
      assert((d.c->olap_element[b] >= 0) == (bi.grc.CountOverlaps(*d.c) > 0));
      
      // update
      if (d.r->olap_element[b] >= 0)
	++in;
      else
	++out;

      if (d.c->olap_element[b] >= 0)
	++in;
      else
	++out;
      
    }
  }
      
  std::cerr << "\tTrack " << b << " In: " << SeqLib::AddCommas(in)
	    << " Out: " << out << " Total track width: " 
	    << SeqLib::AddCommas(bi.grc.TotalWidth()) << std::endl;
}


std::string Matrix::OutputOverlapsIntraExclusive() const {

  if (!size()) // empty matrix
    return std::string();

  std::stringstream ss;

  // run the A==A tracks (eg SINE-SINE)
  for (const auto& a : bed_list) {
    size_t in = 0, out = 0;
    for (const auto& chr : m_vec) {
      for (const auto& d : chr) {
	if (d.r->olap_element[a] == d.c->olap_element[a] && d.r->olap_element[a] >= 0)
	  ++in; // overlap is good
	else
	  ++out;
      }
    }
    ss << a << "\t" << a << "\t" << in << "\t" << out << "\t" << id << std::endl;
  }

  return ss.str();

}

std::string Matrix::OutputOverlapsInterExclusive() const {

  if (!size()) // empty matrix
    return std::string();

  std::stringstream ss;

  // run the A==A tracks (eg SINE-SINE)
  for (const auto& a : bed_list) {
    for (const auto& b : bed_list) {
      if (a == b) { 
	size_t in = 0, out = 0;
	for (const auto& chr : m_vec) {
	  for (const auto& d : chr) {
	    if (d.r->olap_element[a] >= 0 && d.c->olap_element[a] >= 0)
	      ++in; // overlap is good
	    else
	      ++out;
	  }
	}
	ss << a << "\t" << a << "\t" << in << "\t" << out << "\t" << id << std::endl;
      }
    }
  }
  
  // not A != B ones
  for (const auto& a : bed_list) {
    for (const auto& b : bed_list) {
      if (a < b) { // only take half of them (so don't do both SINE/LINE and LINE/SINE)
	size_t in = 0, out = 0;
	for (const auto& chr : m_vec) {
	  for (const auto& d : chr) {
	    // check if (point R in bedA && C in bedB) || (C in bedA && R in bedB)
	    // aelso since this is exclusive, check that they are not in the SAME elem
	    bool AB = d.r->olap_element[a] >= 0 && d.c->olap_element[b] >= 0;
	    bool BA = d.r->olap_element[b] >= 0 && d.c->olap_element[a] >= 0;
	    if (AB || BA)
	      ++in; // overlap is good
	    else
	      ++out;
	  }
	}
	ss << a << "\t" << b << "\t" << in << "\t" << out << "\t" << id << std::endl;
      }
    }
  }

  return ss.str();
  
}
