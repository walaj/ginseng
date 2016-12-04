#include "swap.h"

#include <random>
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "gzstream.h"

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include <prng_engine.hpp>
#include "bmp.h"

#define EXIT_ERROR(msg) { std::cerr << msg << std::endl; exit(EXIT_FAILURE); }

using namespace std;

vector<Matrix*> all_mats;
static pthread_mutex_t swap_lock;

static SeqLib::GRC blacklist;

namespace opt {
  
  static std::vector<std::string> identifiers;
  static std::string analysis_id = "no_id";

  static int min_rar_size = 1000;
  static int max_rar_size = 250e6;

  static std::string header_bam;

  static S num_steps = 10;
  static size_t verbose = 1;
  static size_t numThreads = 1;
  static S num_matrices = 1;
  static S num_bins = 100;
  static size_t half_life = 1000;
  static size_t anim_step = 0;
  static std::string input;

  static std::string blacklist_file;
  static bool intra_only = false;

  static std::string bed_list;

  static int seed = 42;

  static std::string mask;

  static bool inter_only = false;
};

enum {
  OPT_MINSIZE,
  OPT_MAXSIZE,
  OPT_BLACKLIST
};

static const char* shortopts = "w:S:hr:v:P:F:n:k:e:c:n:p:b:t:a:i:x:y:m:s:B:IA:R";
static const struct option longopts[] = {
  { "help",               no_argument, NULL, 'h' },
  { "input-vcf-list",     required_argument, NULL, 'i' },
  { "store-matrices",     required_argument, NULL, 's' },
  { "sub-id",             required_argument, NULL, 'w' },
  { "verbose",            required_argument, NULL, 'v' },
  { "anim-step",          required_argument, NULL, 'A' },
  { "analysis-id",        required_argument, NULL, 'a' },
  { "num-steps",          required_argument, NULL, 'k' },
  { "frac-inter",         required_argument, NULL, 'F' },
  { "num-matrices",       required_argument, NULL, 'n' },
  { "min-span",           required_argument, NULL, OPT_MINSIZE },
  { "max-span",           required_argument, NULL, OPT_MAXSIZE },
  { "inter-chr-only",     no_argument, NULL, 'I' },
  { "intra-chr-only",     no_argument, NULL, 'R' },
  { "num-threads",        required_argument, NULL, 'p' },
  { "power-law",          required_argument, NULL, 'P' },
  { "seed",               required_argument, NULL, 'S' },
  { "half-life",          required_argument, NULL, 't' },
  { "num-bins",           required_argument, NULL, 'b' },
  { "bed-list",           required_argument, NULL, 'B' },
  { "num-events",         required_argument, NULL, 'e' },
  { "bed-file-mask",      required_argument, NULL, 'm' },
  { "blacklist",             required_argument, NULL, OPT_BLACKLIST },
  { NULL, 0, NULL, 0 }
};

static const char *MATRIX_USAGE_MESSAGE =
"Usage: swap -r <num_rows> -c <num_cols>\n\n"
"  Description: Permute rearrangement matrices with preserved row/column marginals and approximately preserved \n"
"\n"
"  General options\n"
"  -h, --help                           Display this help and exit\n"
"  -a, --analysis-id                    Unique string ID to prepend to output files\n"
"  -i, --input-vcf-list                 Txt file with list of VCF files\n"
"  -w, --sub-id                         A string identifer for the VCF files list. If -w specified, line from file from -i must have one of the -w strings to continue.\n"
"  -k, --num-steps                      Number of steps\n"
"  -n, --num-matrices                   Number of matrices to produce\n"
"  -b, --num-bins                       Number of histogram bins for binning length distribution\n"
"  -t, --half-life                      Half-life (in number of steps) for the temperature function.\n"
"  -a, --anim-step                      Number of steps before outputting an animation snapshot (animX.csv). Default 0 (OFF)\n"
"  -m, --event-mask                     BED file to mask events from. If a breakpoint falls in this mask, it is not included.\n"
"  -s, --stored-matrices                Load matrices stored on disk and test hypotheses\n"
"  -I, --inter-chr-only                 Run only inter-chromosomal events. Much much faster, but lower power.\n" 
"      --min-span                       Minimum rearrangement span to accept [1000]\n" 
"      --max-span                       Maximum rearrangement span to accept [250e6]\n" 
"  -I, --inter-chr-only                 Run only inter-chromosomal events. Much much faster, but lower power.\n" 
"  -B, --blacklist                      BED-file with blacklisted regions to not extract variants reads from.\n"
"\n";


bool __header_has_chr_prefix(bam_hdr_t * h) {
  for (int i = 0; i < h->n_targets; ++i) 
    if (h->target_name[i] && std::string(h->target_name[i]).find("chr") != std::string::npos) 
      return true;
  return false;
}

int runSwap(int argc, char** argv) {

  // open a mutex
  if (pthread_mutex_init(&swap_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return false;
  }
  
  parseMatrixOptions(argc, argv);
  
  printinfo();

  // set the random seed
  std::srand(opt::seed);

  // read in the BED files
  BEDMap all_bed;
  import_bed_files(opt::bed_list, all_bed);

  Matrix *m = nullptr;
  
  // read in events from a list of BEDPE
  m = new Matrix(); //(opt::input, opt::num_bins, opt::num_steps, grv_m, opt::inter_only, opt::identifiers, opt::analysis_id, opt::min_rar_size, opt::max_rar_size, blacklist, opt::intra_only);
  m->SetInterOnly(opt::inter_only);
  m->SetInterOnly(opt::inter_only);
  m->SetMinSize(opt::min_rar_size);
  m->SetMaxSize(opt::max_rar_size);
  m->SetNumSwapSteps(opt::num_steps);
  m->id = 0;

  int ret = m->LoadBEDPE(opt::input);
  if (!ret)
    EXIT_ERROR("Unabled to load file: " + opt::input);
  
  // loop the list file and add the BEDPE/VCFs
  /*igzstream inFile(opt::input.c_str());
  if (!inFile) 
    EXIT_ERROR("Can't load the file list: " + opt::input);
  size_t fc = 0;
  std::string fileline;
  while (std::getline(inFile, fileline)) {
    bool ret = true;
    if (fileline.empty() || fileline.at(0) == '#')
      continue;
    if (fileline.find(".bedpe") != std::string::npos)
      ret = m->LoadBEDPE(fileline);
    else if (fileline.find(".vcf") != std::string::npos)
      ret = m->LoadVCF(fileline);
    else 
      std::cerr << "skipping non BEDPE/VCF file: " << fileline << std::endl;
    
    if (!ret)
      EXIT_ERROR("Unabled to load file: " + fileline);
    
    ++fc;
    if (fc % 500 == 0 && opt::verbose)
      std::cerr << "...loaded " << SeqLib::AddCommas(fc) << "th file: " << fileline << std::endl;
  }
  */

  // scramble it
  m->ScrambleAroundDiagonal();
  // fill the initial histogram
  m->FillQuantileHistograms(opt::num_bins);
  std::cerr << "...read in " << SeqLib::AddCommas(m->size()) << " rearrangements" << std::endl;

  // add bed tracks (get overlaps per point)
  std::cerr << "...overlapping BED tracks with 1D points" << std::endl;
  for (auto& b : all_bed)
    m->AddBedElements(b.first, b.second);

  std::cerr << "...creating image" << std::endl;
  const int wid = 2000;
  const int step = 40;
  uint8_t* image = m->AsBMPImage(wid, step);
  const std::string bmpname = opt::analysis_id + ".orig.bmp";
  write_bmp(image, wid, wid, bmpname.c_str());
  
  /*
  GifWriter gw;
  GifBegin(&gw, "test.gif", wid, wid, 1);
  GifWriteFrame(&gw, image, wid, wid, 1);
  GifEnd(&gw);
  */

  // Pre-randomize the chromosomes
  if (!opt::inter_only && m->GetFractionInterChromosomal() < 1) {
    std::cerr << "...pre-computing which chromsomes to swap on for each step" << std::endl;
    PrerandomizeChromosomes(m);
  }  

  // calculate the temps
  if (!opt::inter_only && m->GetFractionInterChromosomal() < 1) {
    std::cerr << "...precomputing probabilities" << std::endl;
    PrecalculateTemps(m);
  }
  
  // calculate the histogram bins
  std::cerr << "...precomputing histogram bins" << std::endl;
  PrecalculateHistogramBins(m);
  
  // write the matrix out and original histogram
  std::cerr << "...writing matrices" << std::endl;
  std::ofstream initial;
  initial.open(opt::analysis_id + ".original.csv");
  std::ofstream initial_h;
  initial_h.open(opt::analysis_id + ".original.histogram.csv");
  std::ofstream initial_h_small;
  initial_h_small.open(opt::analysis_id + ".original.histogram.small.csv");
  m->toCSV(initial, initial_h, initial_h_small);
  initial.close();
  initial_h.close();
  initial_h_small.close();
  
  // set up the result files
  std::ofstream inter_results(opt::analysis_id + ".results.interbin.csv"); // results for points between bed intervals
  std::ofstream intra_results(opt::analysis_id + ".results.intrabin.csv");  // results for points contained in bed interval

  // get the overlaps for the original matrix
  inter_results << m->OutputOverlapsInter(true); // true for exclusive (e.g. don't count SINE-SINE when SINE is exact same one)
  intra_results << m->OutputOverlapsIntraExclusive();

  // set the animation step
  if (opt::anim_step > 0)
    m->setAnimationStep(opt::anim_step);

  // setup the work queue (one matrix per job)
  WorkQueue<SwapWorkItem*>  queue;
  for (size_t i = 0; i < opt::num_matrices; i++) {
    SwapWorkItem * item = new SwapWorkItem(m, i+1, &all_bed);
    queue.add(item);  
  }

  // Create the queue and consumer (worker) threads
  std::vector<ConsumerThread<SwapWorkItem, SwapThreadItem>* > threadqueue;
  
  std::cerr << "...starting the swaps" << std::endl;
  // start the threads and add the jobs
  for (size_t i = 0; i < opt::numThreads; ++i) {
    SwapThreadItem * tu = new SwapThreadItem(i);
    ConsumerThread<SwapWorkItem, SwapThreadItem>* threadr = 
      new ConsumerThread<SwapWorkItem, SwapThreadItem>(queue, tu);
    threadr->start();
    threadqueue.push_back(threadr);
  }
  
  // wait for the threads to finish
  for (size_t i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();

  // write the outputs from each thread 
  for (const auto& i : threadqueue) {
    inter_results << i->GetThreadData()->inter_results.str();
    intra_results << i->GetThreadData()->intra_results.str();
  }
    
  // close the results text files
  inter_results.close();
  intra_results.close();
  
  // deallocate the precomputed probabilities
  if (m->probs) {
    for (size_t i = 0; i < 4; ++i)
      if (m->probs[i])
	free(m->probs[i]);
    free (m->probs);
  }
  
  return 0;
}

void parseMatrixOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'h': die = true; break;
    case 'I': opt::inter_only = true; break;
    case 'R': opt::intra_only = true; break;
      case 'k': arg >> opt::num_steps; break;
      case OPT_BLACKLIST: arg >> opt::blacklist_file; break;
      case OPT_MINSIZE: arg >> opt::min_rar_size; break;
      case OPT_MAXSIZE: arg >> opt::max_rar_size; break;
      case 'B': arg >> opt::bed_list; break;
      case 'v': arg >> opt::verbose; break;
      case 'i': arg >> opt::input; break;
      case 'A': arg >> opt::anim_step; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'w': 
	tmp = std::string();
	arg >> tmp;
	opt::identifiers.push_back(tmp);
	break;
      case 'S': arg >> opt::seed; break;
      case 'n': arg >> opt::num_matrices; break;
      case 'b': arg >> opt::num_bins; break;
      case 't': arg >> opt::half_life; break;
      case 'p': arg >> opt::numThreads; break;
      case 'm': arg >> opt::mask; break;
    }
  }
  
  opt::numThreads = opt::numThreads > opt::num_matrices ? opt::num_matrices : opt::numThreads;

  //if (opt::num_events >= opt::nr*opt::nc*0.5 && !die && opt::input.length() == 0) {
  //  cerr << "Too many events. Lower to less than 50% of bins in matrix" << endl;
  //  exit(EXIT_FAILURE);
  //}

  /*if (opt::anim_step > 0 && opt::num_matrices > 1) {
    cerr << "Animation only available when generating a single matrix" << endl;
    exit(EXIT_FAILURE);
    }*/


  if (die) {
      cerr << "\n" << MATRIX_USAGE_MESSAGE;
      exit(1);
    }
}



void readStoredMatrices(const std::string& file) 
{

  igzstream infile(file.c_str());
  if (!infile) {
    std::cerr << "Can't read stored matrix file " << file << std::endl;
    return;
  }

  std::string line;
  int curr_id = -1, id = 0;
  Matrix * mym = nullptr;
  while(std::getline(infile, line, '\n')) 
    {
      std::istringstream iss(line);
      std::string val;
      size_t count = 0;
      int c_chr = 0, r_chr = 0, r = 0, c = 0;
      
      while(std::getline(iss, val, '\t'))
	{
	  try {
	    switch(count++) {
	    case 0: r_chr = std::stoi(val); break;
	    case 1: r = std::stoi(val); break;
	    case 2: c_chr = std::stoi(val); break;
	    case 3: c = std::stoi(val); break;
	    case 4: id = std::stoi(val); break;
	    }
	  } catch (...) {
	    std::cerr << "Caught stoi error on line " << line << " val " << val << std::endl;
	  }
	}

      //if (id > 10)
      //	return;

      // should be make a new matrix?
      if (id != curr_id) {
	
	// add the old one to pile
	if (id != -1)
	  all_mats.push_back(mym);
	
	// create a new one
	std::cerr << "...creating new matrix with id " << id << std::endl;
	mym = new Matrix();
	for (size_t i = 0; i < 26; i++)
	  mym->m_vec.push_back({});
  
	curr_id = id;
	mym->id = id;
      } 
      
      // add the point
      MatrixValue mv(r_chr, r, c_chr, c);
      mym->addMatrixValue(mv);

    }
  
  

}

// import BED files
void import_bed_files(const std::string& bed_list, BEDMap& all_bed) {

  if (bed_list.empty())
    return;

  std::cerr << "...Importing track BED files" << std::endl;
  
  igzstream ibl(bed_list.c_str());
  if (!ibl) 
    EXIT_ERROR("Can't read bed list file: " + opt::bed_list); 
  
  std::string line;
  while(std::getline(ibl, line, '\n')) {
    if (line.find("#") != std::string::npos) 
      continue;
    std::istringstream iss(line);
    std::string val;
    std::string bedid;
    size_t count=0;
    
    while(getline(iss, val, ',')) {
      ++count;
      if (count == 1)
	bedid = val;
      else {
	all_bed[bedid] = BEDIntervals();
	if (!all_bed[bedid].grc.ReadBED(val, SeqLib::BamHeader())) {
	  std::cerr << "BED file: " << val << " could not be read" << std::endl;
	  exit(EXIT_FAILURE);
	}
	if (!all_bed[bedid].grc.size()) {
	  std::cerr << "BED file: " << val << " is empty before tree creation" << std::endl;
	  exit(EXIT_FAILURE);
	}
	all_bed[bedid].grc.MergeOverlappingIntervals(); // don't allow overlaps. create interval tree
	all_bed[bedid].grc.CreateTreeMap();
	std::cerr << "\tread " << bedid << " with " << SeqLib::AddCommas(all_bed[bedid].size()) << " regions " << std::endl;
      }
    } // end intra-line loop
  } // end line loop
    
}

void PrecalculateTemps(Matrix* m) {

  double frac_inter = m->GetFractionInterChromosomal();
  
  // pre-compute the probabilities given a shift (of 4 possible shifts away from optimal)
  // deallocation will take place at end of swap.cpp
  uint16_t** probs = (uint16_t**)malloc(4 * sizeof(uint16_t*));  
  probs[0] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[1] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[2] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[3] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  
  // pre-compute the temps
  size_t non_trans = 0;
  double tempr = TMAX; 
  if (frac_inter < 1 && opt::half_life && (double)opt::half_life/double(opt::num_steps) < 1) { // if half_life is zero, no temperature decay
    
    for (size_t i = 0; i < opt::num_steps; ++i) {
      
      if (m->rand_chr[i] == 25) // translocation swap, so skip cooling etc
	continue;
      else
	++non_trans;
      
      if (non_trans == 1) // first one
	std::cerr << "...swaps 0 to " << SeqLib::AddCommas(i) << " are translocation swaps" << std::endl;

      if (tempr >= 1) {
	double frac = (double)non_trans/(double)opt::half_life;
	tempr = TMAX*pow(2,-frac);
      } else {       // if too cold, then skip rest of compute
	std::cerr << "...filling cold (0 prob) from " << SeqLib::AddCommas(i)
		  << " to " << SeqLib::AddCommas(opt::num_steps) << std::endl;
	m->probs = probs; // the rest already 0 bc of calloc
	return;
      }
      probs[0][i] = (uint16_t)std::min(std::floor(exp(-1/tempr)*TRAND), (double)TRAND);
      probs[1][i] = (uint16_t)std::min(std::floor(exp(-2/tempr)*TRAND), (double)TRAND);
      probs[2][i] = (uint16_t)std::min(std::floor(exp(-3/tempr)*TRAND), (double)TRAND);
      probs[3][i] = (uint16_t)std::min(std::floor(exp(-4/tempr)*TRAND), (double)TRAND);
    }
    // half life is huge, so we don't want to decay at all (permanently hot)
  } else if ((double)opt::half_life/double(opt::num_steps) >= 1) { 
    for (size_t i = 0; i < opt::num_steps; ++i) {
      probs[0][i] = TRAND;
      probs[1][i] = TRAND;
      probs[2][i] = TRAND;
      probs[3][i] = TRAND;
    }    
  }
  
  m->probs = probs;

}

void PrecalculateHistogramBins(Matrix* m) {
  // precompute the histogram bins
  uint32_t max_dist = 250000000;
  uint8_t * bin_table = (uint8_t*) calloc(max_dist, sizeof(uint8_t));
  for (size_t i = 1; i < max_dist; ++i) {
    std::vector<int32_t>::const_iterator it = std::upper_bound(m->hist.begin(), m->hist.end(), i);
    bin_table[i] = (it - m->hist.begin() - 1);
  }
  m->bin_table = bin_table;  
}

void PrerandomizeChromosomes(Matrix* m) {
  
  sitmo::prng_engine eng1;
  eng1.seed(1337);
  
  double frac_inter = m->GetFractionInterChromosomal();
  const float inter_factor = 0.5; // the higher this number, the more translocation swaps
  size_t num_intra = frac_inter < 1 ? std::floor((double)opt::num_steps * ((double)1 - inter_factor * frac_inter)) : opt::num_steps;
  uint8_t * rand_chr = (uint8_t*) calloc(opt::num_steps, sizeof(uint8_t));
  std::cerr << "...fraction inter " << frac_inter << " -- number of swaps to make intra-chromosomal: " << SeqLib::AddCommas(num_intra) << std::endl;
  
  for (size_t i = 0; i < opt::num_steps; ++i) {
    size_t rv = eng1() % m->GetNumIntraChromosomal(); 
    size_t running_count = 0;
    
    int chr = 25; // start with inter
    if (i > (opt::num_steps - num_intra) || opt::intra_only) 
      {
	for (size_t j = 0; j < 24; ++j) 
	  {
	    size_t rc2 = running_count + m->m_vec[j].size();
	    if (rv >= running_count && rv <= rc2 && m->m_vec[j].size()) 
	      chr = j;
	    running_count = rc2;
	  }
      }
    rand_chr[i] = chr;
  }
  m->rand_chr = rand_chr;

}


void printinfo() {

  std::cerr << "Input VCF list: " << opt::input << endl;
  std::cerr << "ID:             " << opt::analysis_id << std::endl;
  std::cerr << "Matrices:       " << SeqLib::AddCommas(opt::num_matrices) << std::endl;
  std::cerr << "Steps:          " << SeqLib::AddCommas(opt::num_steps) << std::endl;
  std::cerr << "Histogram bins: " << opt::num_bins << endl;
  std::cerr << "Temp 1/2 life:  " << SeqLib::AddCommas(opt::half_life) << endl;
  std::cerr << "Animation:      " << (opt::anim_step > 0 ? SeqLib::AddCommas(opt::anim_step) : "OFF")  << endl;
  std::cerr << "BED list:       " << (opt::bed_list) << endl;
  
  if (opt::identifiers.size()) {
    std::cerr << "Identifiers to trim to" << std::endl;
    for (auto& i : opt::identifiers)
      std::cerr << "\t" << i;
    std::cerr << std::endl;
  }
  
  if (opt::half_life == 0) 
    std::cerr << "TEMPERATURE SET TO ZERO" << std::endl;
  else if ((double)opt::half_life/(double)opt::num_steps >= 1)
    std::cerr << "TEMPERATURE SET TO INFINITE" << std::endl;      
  
  if (opt::inter_only)
    std::cerr << "--------------------------------------------\n" << 
      "      RUNNING AS INTER CHROMOSOMAL ONLY     " << std::endl 
	      << "--------------------------------------------" << std::endl;
  if (opt::intra_only)
    std::cerr << "--------------------------------------------\n" << 
      "      RUNNING AS INTRA CHROMOSOMAL ONLY     " << std::endl 
	      << "--------------------------------------------" << std::endl;
  else 
    std::cerr << " Rearrangement size bounds: [" << SeqLib::AddCommas(opt::min_rar_size) << " - " << SeqLib::AddCommas(opt::max_rar_size) << "]" << std::endl;
  
  if (!opt::blacklist_file.empty()) {
    blacklist = SeqLib::GRC(opt::blacklist_file, SeqLib::BamHeader());
    blacklist.CreateTreeMap();
    if (opt::verbose) std::cerr << "...read in blacklist " << blacklist.size() << std::endl;
  }

}
