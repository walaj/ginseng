#include "swap.h"

#include <random>
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "SnowTools/gzstream.h"
#include <regex>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"

#include <prng_engine.hpp>

using namespace std;

vector<Matrix*> all_mats;
static pthread_mutex_t swap_lock;

namespace opt {

  static S num_steps = 10;
  static size_t verbose = 1;
  static size_t numThreads = 1;
  static S num_matrices = 1;
  static S num_bins = 100;
  static size_t half_life = 1000;
  static size_t anim_step = 0;
  static string input = "";
  static string stored_matrices = "";

  static string bed_list = "";

  static string bedA;
  static string bedB;
  static string mask;

  static bool inter_only = false;

};

static const char* shortopts = "hr:v:n:k:e:c:n:p:b:t:a:i:x:y:m:s:B:I";
static const struct option longopts[] = {
  { "help",               no_argument, NULL, 'h' },
  { "input-vcf-list",     required_argument, NULL, 'i' },
  { "store-matrices",     required_argument, NULL, 's' },
  { "verbose",            required_argument, NULL, 'v' },
  { "anim-step",          required_argument, NULL, 'a' },
  { "num-steps",          required_argument, NULL, 'k' },
  { "num-matrices",       required_argument, NULL, 'n' },
  { "inter-chr-only",        no_argument, NULL, 'I' },
  { "num-threads",        required_argument, NULL, 'p' },
  { "half-life",          required_argument, NULL, 't' },
  { "num-bins",           required_argument, NULL, 'b' },
  { "bed-list",           required_argument, NULL, 'B' },
  { "bed-file-mask",           required_argument, NULL, 'm' },
  { NULL, 0, NULL, 0 }
};

static const char *MATRIX_USAGE_MESSAGE =
"Usage: swap -r <num_rows> -c <num_cols>\n\n"
"  Description: Permute rearrangement matrices with preserved row/column marginals and approximately preserved \n"
"\n"
"  General options\n"
"  -h, --help                           Display this help and exit\n"
"  -i, --input-vcf-list                 Txt file with list of VCF files\n"
"  -k, --num-steps                      Number of steps\n"
"  -n, --num-matrices                   Number of matrices to produce\n"
"  -b, --num-bins                       Number of histogram bins for binning length distribution\n"
"  -t, --half-life                      Half-life (in number of steps) for the temperature function.\n"
"  -a, --anim-step                      Number of steps before outputting an animation snapshot (animX.csv). Default 0 (OFF)\n"
"  -m, --event-mask                     BED file to mask events from. If a breakpoint falls in this mask, it is not included.\n"
"  -s, --stored-matrices                Load matrices stored on disk and test hypotheses\n"
"  -I, --inter-chr-only                 Run only inter-chromosomal events. Much much faster, but lower power.\n" 
"\n";

int main(int argc, char** argv) {

#ifdef __APPLE__
  cout << "CLOCK_GETTIME NOT AVAILABLE ON MAC" << endl;
#else
  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // open a mutex
  if (pthread_mutex_init(&swap_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return false;
  }

  parseMatrixOptions(argc, argv);

  if (opt::verbose) {
    if (opt::input.length())
      cout << "Input VCF list: " << opt::input << endl;
    else if (opt::stored_matrices.length()) {
      cout << "Stored matrix file " << opt::stored_matrices << endl;
    } 
    cout << "Num matrices:  " << opt::num_matrices << endl;
    cout << "Num steps:     " << opt::num_steps << endl;
    cout << "Num hist bins: " << opt::num_bins << endl;
    cout << "Num threads:   " << opt::numThreads << endl;
    cout << "Temp 1/2 life: " << opt::half_life << endl;
    cout << "Animation:     " << (opt::anim_step > 0 ? to_string(opt::anim_step) : "OFF")  << endl;
    cout << "BED list:      " << (opt::bed_list) << endl;
    cout << "BED file mask: " << (opt::mask.length() ? opt::mask : "NONE") << endl;
    
    if (opt::half_life == 0) 
      std::cout << "TEMPERATURE SET TO ZERO" << std::endl;
    else if ((double)opt::half_life/(double)opt::num_steps >= 1)
      std::cout << "TEMPERATURE SET TO INFINITE" << std::endl;      

    if (opt::inter_only)
      std::cerr << "--------------------------------------------\n" << 
	"        RUNNING AS INTER CHROMOSOMAL ONLY     " << std::endl 
		<< "--------------------------------------------" << std::endl;
  }

  std::unordered_map<string, SnowTools::GRC> all_bed;
  // read the bed files
  if (opt::bed_list.length()) {
    if (opt::verbose) 
      std::cout << "...Importing BED files" << std::endl;

    igzstream ibl(opt::bed_list.c_str());
    if (!ibl) { std::cerr << "Can't read bed list file: " << opt::bed_list << std::endl; return 1; }
    
    string line;
    while(std::getline(ibl, line, '\n')) 
      {
	if (line.find("#") == std::string::npos) {
	  std::istringstream iss(line);
	  string val;
	  string bedid;
	  size_t count=0;
	  while(getline(iss, val, ',')) {
	    ++count;
	    if (count == 1)
	      bedid = val;
	    else {
	      //std::cout << "...reading in " << val << std::endl;
	      all_bed[bedid] = SnowTools::GRC();
	      all_bed[bedid].regionFileToGRV(val);
	      all_bed[bedid].createTreeMap();
	      std::cout << "...read in " << bedid << " with " << all_bed[bedid].size() << " regions " << std::endl;
	    }
	  }
	}
      }
  }
  
  // read in the mask
  SnowTools::GRC grv_m;
  if (opt::mask.length()) {
    grv_m.regionFileToGRV(opt::mask);
    grv_m.createTreeMap();
    std::cout << "Read in " << grv_m.size() << " region from the mask file " << opt::mask << endl;
  }

  Matrix *m = nullptr;

  if (opt::input.length())
    // read in events from a list of VCFs
    m = new Matrix(opt::input, opt::num_bins, opt::num_steps, grv_m, opt::inter_only);
  //else if (!opt::stored_matrices.length())
    // make a random matrix
  //  m = new Matrix(opt::nr, opt::nc, opt::num_events, opt::num_bins, opt::num_steps);
  else if (opt::stored_matrices.length()) {
    std::cerr << "...reading stored matrices" << std::endl;
    readStoredMatrices(opt::stored_matrices);
  }

  //if (all_mats.size() == 0)
  //  opt::num_steps = 1;

  // pre-compute the probabilities given a shift (of 4 possible shifts away from optimal)
  uint16_t* probs[4];  
  probs[0] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[1] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[2] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  probs[3] = (uint16_t*) calloc(opt::num_steps, sizeof(uint16_t));
  
  std::cerr << "ALL_MATS.SIZE() " << all_mats.size() << std::endl;
  std::cerr << "INTER_ONLY: " << opt::inter_only << std::endl;

  if (all_mats.size() == 0) {

    // pre-compute the temps
    //double * temps = (double*) calloc(opt::num_steps, sizeof(double));
    if (opt::half_life && (double)opt::half_life/double(opt::num_steps) < 1) { // if half_life is zero, no temperature decay

      std::cerr << "...computing probabilities" << std::endl;

      for (size_t i = 0; i < opt::num_steps; ++i) {
	double frac = (double)i/(double)opt::half_life;
	//temps[i] = TMAX*pow(2,-frac);
	double tempr = TMAX*pow(2,-frac);
	probs[0][i] = (uint16_t)std::min(std::floor(exp(-1/tempr)*TRAND), (double)TRAND);
	probs[1][i] = (uint16_t)std::min(std::floor(exp(-2/tempr)*TRAND), (double)TRAND);
	probs[2][i] = (uint16_t)std::min(std::floor(exp(-3/tempr)*TRAND), (double)TRAND);
	probs[3][i] = (uint16_t)std::min(std::floor(exp(-4/tempr)*TRAND), (double)TRAND);
	if (i % 2000000 == 0)
	  std::cout << "P(least-non-optimal) " << probs[0][i] << " P(most-non-optimal) " << probs[3][i] << " Temp " << tempr << " Step " << SnowTools::AddCommas<size_t>(i) << std::endl;
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
    if (!opt::inter_only) {
      std::cout << "...precomputing bin indicies" << std::endl;
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
    if (!opt::inter_only) {
      std::cout << "...pre-computing which chromsomes to swap on for each step" << std::endl;
      sitmo::prng_engine eng1;
      eng1.seed(1337);

      size_t num_intra = (size_t)std::floor((double)opt::num_steps *0.97);
      uint8_t * rand_chr = (uint8_t*) calloc(opt::num_steps, sizeof(uint8_t));
      std::cerr << "...number of swaps to make intra-chromosomal: " << SnowTools::AddCommas(num_intra) << std::endl;

      for (size_t i = 0; i < opt::num_steps; ++i) {
	size_t rv = eng1() % m->m_intra;    
	size_t running_count = 0;
	int chr = 25; // start with intra
	if (i < num_intra) 
	  {
	    for (size_t j = 0; j < 24; ++j) 
	      {
		size_t rc2 = running_count + m->m_vec[j].size();
		if (rv >= running_count && rv <= rc2 && m->m_vec[j].size()) 
		  {
		    chr = j;
		    if (chr > 25) {//debug 
		      std::cerr << "error " << chr << std::endl; exit(1); }
		    break;
		  } 
		running_count = rc2;
	      }
	  }
	rand_chr[i] = chr;
      }
      m->rand_chr = rand_chr;
    }  
  } // done with check all_mats.size()
  
  // output the original matrix
  std::cout << "...outputting original data to csv" << std::endl;
  std::ofstream initial;
  initial.open("original.csv");
  std::ofstream initial_h;
  initial_h.open("original.histogram.csv");
  std::ofstream initial_h_small;
  initial_h_small.open("original.histogram.small.csv");

  m->toCSV(initial, initial_h, initial_h_small);
  initial.close();
  initial_h.close();
  initial_h_small.close();

  if (opt::verbose) {
    cout << "...done importing" << endl;
#ifndef __APPLE__
    SnowTools::displayRuntime(start);
    cout << endl;
#endif
  }

  // set up the result files
  std::ofstream results("results.csv");
  
  // write the original overlaps
  std::unordered_map<std::string, OverlapResult> all_overlaps;
  std::unordered_map<std::string, bool> ovl_ovl_seen;
  for (auto& i : all_bed) {
    for (auto& j : all_bed) {
      if (i.first < j.first || (!ovl_ovl_seen.count(i.first) && i.first==j.first) ) { // don't need to do each one twice
	OverlapResult ovl = m->checkOverlaps(&i.second, &j.second);
	std::string ovl_name = i.first + "," + j.first;
	all_overlaps[ovl_name] = ovl;
	std::cout << " ORIGINAL OVERLAP for " << ovl_name << " is "  << ovl.first << " ORIGINAL NO OVERLAP " << ovl.second << std::endl;      
	results << ovl_name << "," << ovl.first << "," << ovl.second << ",-1" << std::endl; 
	if (i.first==j.first)
	  ovl_ovl_seen[i.first] = true;
      }
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
	    //std::cout << " Overlap for "  << ovl_name << " is "  << ovl.first << std::endl;      
	    results << ovl_name << "," << ovl.first << "," << ovl.second << "," << mmm->m_id  << std::endl; 
	    if (i.first==j.first)
	      ovl_ovl_seen[i.first] = true;
	  }
	}
      }
      
    }
  }
    
  // set the output file
  ogzstream * oz_matrix = nullptr;
  std::cerr << "all_mats " << all_mats.size() << std::endl;
  if (all_mats.size() == 0) {
    std::string namr = "all.matrices.csv.gz";
    oz_matrix = new ogzstream(namr.c_str(), std::ios::out);
  }

  // make a new empty summed results matrix
  Matrix * s_results = new Matrix();
  
  vector<SwapWorkItem*> tmp_queue;
  for (S i = 0; i < opt::num_matrices; ++i) {
    SwapWorkItem * item = new SwapWorkItem(m, &all_mats, &swap_lock, i, &all_bed, &results, oz_matrix, s_results);
    tmp_queue.push_back(item);
  }

  if (opt::verbose)
    cout << "Sending matrices to threads" << endl;
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
  s_results->toSimpleCSV("summed_out.csv");
  
  // free the temps and probabilities
  //free(temps);
  free(probs[0]);
  free(probs[1]);
  free(probs[2]);
  free(probs[3]);

  if (opt::verbose) {
#ifndef __APPLE__
    SnowTools::displayRuntime(start);
    cout << endl;
#endif
  }
  //if (opt::verbose)
  //  cout << "Generating Random Matrices using Self-Loop method" << endl;

  //for (size_t i = 0; i < opt::num_matrices; i++) {
  //  if (opt::verbose)
  //    cout << "...generating matrix " << (i+1) << " of " << opt::num_matrices << endl;
  //  m->allSwaps();
  //  all_mats.push_back(m);
  //}
  
  //delete m;
}

void parseMatrixOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'h': die = true; break;
    case 'I': opt::inter_only = true; break;
      case 'k': arg >> opt::num_steps; break;
      case 'B': arg >> opt::bed_list; break;
      case 'v': arg >> opt::verbose; break;
      case 's': arg >> opt::stored_matrices; break;
      case 'i': arg >> opt::input; break;
	//case 'r': arg >> opt::nr; break;
      case 'x': arg >> opt::bedA; break;
      case 'y': arg >> opt::bedB; break;
	//case 'c': arg >> opt::nc; break;
      case 'a': arg >> opt::anim_step; break;
	//case 'e': arg >> opt::num_events; break;
      case 'n': arg >> opt::num_matrices; break;
      case 'b': arg >> opt::num_bins; break;
      case 't': arg >> opt::half_life; break;
      case 'p': arg >> opt::numThreads; break;
      case 'm': arg >> opt::mask; break;
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
      cout << "\n" << MATRIX_USAGE_MESSAGE;
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
  Matrix * mym;
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
	mym->m_id = id;
      } 
      
      // add the point
      MatrixValue mv(r_chr, r, c_chr, c);
      mym->addMatrixValue(mv);

    }
  
  

    }
