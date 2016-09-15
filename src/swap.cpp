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
  static std::string input = "";
  static string stored_matrices = "";

  static std::string blacklist_file;
  static bool intra_only = false;

  static std::string bed_list = "";

  static int mode = -1;

  static int seed = 42;

  static double frac_inter = 0.2;

  static std::string bedA;
  static std::string bedB;
  static std::string mask;

  static bool inter_only = false;

  static double power_law = 1;

  static int num_events = 1000;
};

enum {
  OPT_SIMULATE,
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
  { "simulate",           no_argument, NULL, OPT_SIMULATE },
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
"  Simulation options (--simulate)\n"
"  -e, --num-events                     Number of events to simulate.\n" 
"  -P, --power-law                      Power law parameter from which to draw lengths from (1^(-P)). Default 1.0\n" 
"  -F, --frac-inter                     Fraction of events to be inter-chromosomal. Default 0.2\n" 
"  -S, --seed                           Seed for the RNG. Default 42\n" 
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

  if (opt::verbose) {
    if (opt::input.length())
      cerr << "Input VCF list: " << opt::input << endl;
    else if (opt::stored_matrices.length()) {
      cerr << "Stored matrix file " << opt::stored_matrices << endl;
    } 
    cerr << "ID: " << opt::analysis_id << "Matrices: " << SeqLib::AddCommas(opt::num_matrices) << "\tSteps: " << SeqLib::AddCommas(opt::num_steps) << "\tHistBins: " << opt::num_bins << endl;
    cerr << "Temp 1/2 life: " << SeqLib::AddCommas(opt::half_life) << endl;
    cerr << "Animation:     " << (opt::anim_step > 0 ? SeqLib::AddCommas(opt::anim_step) : "OFF")  << endl;
    cerr << "BED list:      " << (opt::bed_list) << endl;

    if (opt::identifiers.size()) {
      std::cerr << "Identifiers to trim to" << std::endl;
      for (auto& i : opt::identifiers)
	std::cerr << "\t" << i;
      std::cerr << std::endl;
    }
    //cerr << "BED file mask: " << (opt::mask.length() ? opt::mask : "NONE") << endl;

    if (opt::mode == OPT_SIMULATE) {
      cerr << "--------------- Simulating events ----------------" << std::endl;
      cerr << "Num events: " << SeqLib::AddCommas((opt::num_events)) << "\tPowerLaw: " << opt::power_law << " FracInter: " << opt::frac_inter << std::endl;
      cerr << "--------------------------------------------------" << std::endl;
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
      
  }

  //blacklist.add(SeqLib::GenomicRegion(1,33139671,33143258)); // really nasty region
  if (!opt::blacklist_file.empty()) {
    blacklist = SeqLib::GRC(opt::blacklist_file, SeqLib::BamHeader());
    blacklist.CreateTreeMap();
    if (opt::verbose) std::cerr << "...read in blacklist " << blacklist.size() << std::endl;
  }

  // read in the BED files
  BEDMap all_bed;
  import_bed_files(opt::bed_list, all_bed);

  Matrix *m = nullptr;

  // set the random seed
  std::srand(opt::seed);

  if (opt::input.length() && opt::input.find("csv") == std::string::npos && opt::mode != OPT_SIMULATE) {

    // read in events from a list of BEDPE
    m = new Matrix(); //(opt::input, opt::num_bins, opt::num_steps, grv_m, opt::inter_only, opt::identifiers, opt::analysis_id, opt::min_rar_size, opt::max_rar_size, blacklist, opt::intra_only);
    m->SetInterOnly(opt::inter_only);
    m->SetInterOnly(opt::inter_only);
    m->SetMinSize(opt::min_rar_size);
    m->SetMaxSize(opt::max_rar_size);
    m->SetNumSwapSteps(opt::num_steps);
    m->id = 0;

    // loop the list file and add the BEDPE/VCFs
    igzstream inFile(opt::input.c_str());
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
    // scramble it
    m->ScrambleAroundDiagonal();
    // fill the initial histogram
    m->FillQuantileHistograms(opt::num_bins);

  } else if (opt::stored_matrices.length()) {
    std::cerr << "...reading stored matrices" << std::endl;
    readStoredMatrices(opt::stored_matrices);
  }


  if (opt::verbose)
    std::cerr << "...read in " << SeqLib::AddCommas(m->size()) << " rearrangements" << std::endl;

  double frac_inter = m->GetFractionInterChromosomal();

  // pre-compute the probabilities given a shift (of 4 possible shifts away from optimal)
  uint16_t* probs[4];  
  probs[0] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[1] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[2] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[3] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  
  if (all_mats.size() == 0) {
    
    // pre-compute the temps
    //double * temps = (double*) calloc(opt::num_steps, sizeof(double));
    if (frac_inter < 1 && opt::half_life && (double)opt::half_life/double(opt::num_steps) < 1) { // if half_life is zero, no temperature decay

     for (size_t i = 0; i < opt::num_steps; ++i) {
       
       double frac = (double)i/(double)opt::half_life;
       double tempr = TMAX*pow(2,-frac);
       // if too cold, then skip rest of compute
       if (tempr < 1) {
	 std::cerr << "...filling cold (0 prob) to end of " << SeqLib::AddCommas(opt::num_steps) << std::endl;
	 for (size_t j = i; j < opt::num_steps; ++j) {
	   probs[0][j] = 0; probs[1][j] = 0; probs[2][j] = 0; probs[3][j] = 0;
	 }
	 break;
	 }
       probs[0][i] = (uint16_t)std::min(std::floor(exp(-1/tempr)*TRAND), (double)TRAND);
       probs[1][i] = (uint16_t)std::min(std::floor(exp(-2/tempr)*TRAND), (double)TRAND);
       probs[2][i] = (uint16_t)std::min(std::floor(exp(-3/tempr)*TRAND), (double)TRAND);
       probs[3][i] = (uint16_t)std::min(std::floor(exp(-4/tempr)*TRAND), (double)TRAND);
       if (i % 5000000 == 0)
	 std::cerr << "P(least-non-optimal) " << probs[0][i] << " P(most-non-optimal) " << probs[3][i] << " Temp " << tempr << " Step " << SeqLib::AddCommas<size_t>(i) << std::endl;
     }
      // half life is huge, so we don't want to decay at all (permanently hot)
    } else if ((double)opt::half_life/double(opt::num_steps) >= 1) { 
      for (size_t i = 0; i < opt::num_steps; ++i) {
	//temps[i] = TMAX;
	probs[0][i] = TRAND;
	probs[1][i] = TRAND;
	probs[2][i] = TRAND;
	probs[3][i] = TRAND;
      }    
    }
    m->probs = probs;
    
    // precompute the histogram bins
    if (!opt::inter_only && frac_inter < 1) {
      std::cerr << "...precomputing bin indicies (which span goes to which histogram bin)" << std::endl;
      uint32_t max_dist = 250000000;
      uint8_t * bin_table = (uint8_t*) calloc(max_dist, sizeof(uint8_t));
      for (size_t i = 1; i < max_dist; ++i)
	{
	  std::vector<int32_t>::const_iterator it = std::upper_bound(m->hist.begin(), m->hist.end(), i);
	  bin_table[i] = (it - m->hist.begin() - 1);
	}
      m->bin_table = bin_table;  
    }

    // precompute which chrom to acces
    if (!opt::inter_only && frac_inter < 1) {
      std::cerr << "...pre-computing which chromsomes to swap on for each step" << std::endl;
      sitmo::prng_engine eng1;
      eng1.seed(1337);

      size_t num_intra = frac_inter < 1 ? std::floor((double)opt::num_steps * ((double)1 - 0.5 * frac_inter)) : opt::num_steps;
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
  } // done with check all_mats.size()
  // output the original matrix
  std::ofstream initial;
  initial.open(opt::analysis_id + ".original.csv");
  

  // output the size histograms
  std::ofstream initial_h;
  initial_h.open(opt::analysis_id + ".original.histogram.csv");
  std::ofstream initial_h_small;
  initial_h_small.open(opt::analysis_id + ".original.histogram.small.csv");
  m->toCSV(initial, initial_h, initial_h_small);
  initial.close();
  initial_h.close();
  initial_h_small.close();

  // set up the result files
  std::ofstream results(opt::analysis_id + ".results.csv");
  std::ofstream results2(opt::analysis_id + ".results.intra.csv");
  // write the original overlaps
  std::unordered_map<std::string, OverlapResult> all_overlaps;
  std::unordered_map<std::string, bool> ovl_ovl_seen;
  for (auto& i : all_bed) {
    for (auto& j : all_bed) {
      if (i.first < j.first || (!ovl_ovl_seen.count(i.first) && i.first==j.first) ) { // don't need to do each one twice
	OverlapResult ovl = m->checkOverlaps(&i.second, &j.second);

	std::string ovl_name = i.first + "," + j.first;
	all_overlaps[ovl_name] = ovl;
	//std::cerr << " ORIGINAL OVERLAP for " << ovl_name << " is "  << ovl.first << " ORIGINAL NO OVERLAP " << ovl.second << std::endl;      
	results << ovl_name << "," << ovl.first << "," << ovl.second << ",-1" << std::endl; 
	if (i.first==j.first)
	  ovl_ovl_seen[i.first] = true;
      }
    }
  }
  
  // check intra unit overlaps
  for (auto& i : all_bed) {
    if (do_intra_overlap(i.first)) {
      std::cerr << "...checking intra overlaps for " << i.first;
      OverlapResult ovl = m->checkIntraUnitOverlaps(&i.second);
      results2 << i.first << "," << ovl.first << "," << ovl.second << ",-1" << std::endl; 
      
    }
  }
  
  // Create the queue and consumer (worker) threads
  wqueue<SwapWorkItem*>  queue;
  vector<ConsumerThread<SwapWorkItem>*> threadqueue;
  
  if (all_mats.size() == 0) {
    
    std::cerr << "...running new matrix swaps" << std::endl;
    
    // if num threads too high, set to lower
    opt::numThreads = opt::numThreads > opt::num_matrices ? opt::num_matrices : opt::numThreads;
    
    for (unsigned i = 0; i < opt::numThreads; i++) {
      ConsumerThread<SwapWorkItem>* threadr = new ConsumerThread<SwapWorkItem>(queue, opt::verbose > 0);
      threadr->start();
      threadqueue.push_back(threadr);
    }
    
    // set the animation step
    if (opt::anim_step > 0)
      m->setAnimationStep(opt::anim_step);
    
  } else {
    
    std::cerr << "...checking results for stored matrices" << std::endl;

    for (auto& mmm : all_mats) {
      // loop through all of the stored matrices and do the overlaps
      std::unordered_map<std::string, OverlapResult> all_overlaps;
      std::unordered_map<std::string, bool> ovl_ovl_seen;
      for (auto& i : all_bed) {
	for (auto& j : all_bed) {
	  if (i.first < j.first || (!ovl_ovl_seen.count(i.first) && i.first==j.first) ) { // don't need to do each one twice
	    OverlapResult ovl = mmm->checkOverlaps(&i.second, &j.second);
	    std::string ovl_name = i.first + "," + j.first;
	    all_overlaps[ovl_name] = ovl;
	    //std::cerr << " Overlap for "  << ovl_name << " is "  << ovl.first << std::endl;      
	    results << ovl_name << "," << ovl.first << "," << ovl.second << "," << mmm->id  << std::endl; 
	    if (i.first==j.first)
	      ovl_ovl_seen[i.first] = true;
	  }
	}
      }

      // check intra unit overlaps
      for (auto& i : all_bed) {
	if (do_intra_overlap(i.first)) {
	  OverlapResult ovl = m->checkIntraUnitOverlaps(&i.second);
	  results2 << i.first << "," << ovl.first << "," << ovl.second << mmm->id << std::endl; 
	}
      }

      
    }
  }
    
  // set the output file
  ogzstream * oz_matrix = nullptr;
  if (all_mats.size() == 0) {
    std::string namr = opt::analysis_id + ".all.matrices.csv.gz";
    oz_matrix = new ogzstream(namr.c_str(), std::ios::out);
  }

  // make a new empty summed results matrix
  Matrix * s_results = new Matrix();
  
  vector<SwapWorkItem*> tmp_queue;
  m->grc1.clear();
  m->grc2.clear();
  for (S i = 0; i < opt::num_matrices; ++i) {
    SwapWorkItem * item = new SwapWorkItem(m, &all_mats, &swap_lock, i, &all_bed, &results, &results2, oz_matrix, s_results);
    tmp_queue.push_back(item);
  }

  if (opt::verbose)
    cerr << "Sending matrices to threads" << endl;
  for (size_t i = 0; i < tmp_queue.size(); i++)
    queue.add(tmp_queue[i]);  
  
  // wait for the threads to finish
  for (unsigned i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();

  results.close();
  if (oz_matrix) {
    oz_matrix->close();
    delete oz_matrix;
  }

  // write the summed results
  //s_results->toSimpleCSV("summed_out.csv");
  
  // free the temps and probabilities
  //free(temps);
  free(probs[0]);
  free(probs[1]);
  free(probs[2]);
  free(probs[3]);
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
      case OPT_SIMULATE: opt::mode = OPT_SIMULATE; break;
      case OPT_MINSIZE: arg >> opt::min_rar_size; break;
      case OPT_MAXSIZE: arg >> opt::max_rar_size; break;
      case 'B': arg >> opt::bed_list; break;
      case 'v': arg >> opt::verbose; break;
      case 's': arg >> opt::stored_matrices; break;
      case 'i': arg >> opt::input; break;
	//case 'r': arg >> opt::nr; break;
      case 'x': arg >> opt::bedA; break;
      case 'y': arg >> opt::bedB; break;
	//case 'c': arg >> opt::nc; break;
      case 'A': arg >> opt::anim_step; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'e': arg >> opt::num_events; break;
      case 'w': 
	tmp = std::string();
	arg >> tmp;
	opt::identifiers.push_back(tmp);
	break;
      case 'F': arg >> opt::frac_inter; break;
      case 'S': arg >> opt::seed; break;
      case 'n': arg >> opt::num_matrices; break;
      case 'b': arg >> opt::num_bins; break;
      case 't': arg >> opt::half_life; break;
      case 'p': arg >> opt::numThreads; break;
      case 'm': arg >> opt::mask; break;
      case 'P': arg >> opt::power_law; break;
    }
  }
  
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

  if (opt::verbose) 
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
	all_bed[bedid] = SeqLib::GRC(val, SeqLib::BamHeader());
	all_bed[bedid].CreateTreeMap();
	if (opt::verbose)
	  std::cerr << "\tread " << bedid << " with " << SeqLib::AddCommas(all_bed[bedid].size()) << " regions " << std::endl;
      }
    } // end intra-line loop
  } // end line loop
    
}

