#include "sim.h"
#include <getopt.h>
#include <sstream>

#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <fstream>
#include <chrono>

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamReader.h"
#include "Fractions.h"

const char* SIM_USAGE_MESSAGE = 
"\nProgram: sim - Simulate breakpoints and rearrangements on the genome\n"
"Contact: Jeremiah Wala <jwala@broadinstitute.org>\n"
"Usage: sim [covered.bed] <options> > out.bed\n"
"  General options\n"
"  -h, --help                Display this help and exit\n"
"  -L, --length              Length distribution factor (x^-L) for 2D simulation. [1]\n"
"  -N, --num-sims            Number of points / pairs to attempt simulation. [1e6]\n"
"  -m, --mode                Simulation mode. (1) one-dimensional (A) additive (M) multiplicative. [1]\n"
"  -r, --ratio               Ratio of intrachromosomal events to translocations. [4]\n"
"  -R, --random              Generate random BED file of 1D probabs, with variance -R. 10000 bp bins. Writes to random_bias.bed\n"
  "\n";

static const int MIN = 1000;
static const int MAX = 100000000;
static const int PRECISION = 1000000; // how carefully to parse random nums

static uint32_t GENOMESIZE = 3095693983;
static const uint32_t csum[24] = {0, 249250621, 492449994, 690472424, 881626700, 1062541960,
                             1233657027, 1392795690, 1539159712, 1680373143,
                             1815907890, 1950914406, 
			     2084766301, 2199936179, 2307285719, 2409817111, 2500171864, 2581367074,
			     2659444322, 2718573305, 2781598825, 2829728720, 2881033286, 3036303846};

static const uint32_t csize[24] = {249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
			      141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392,  90354753,
			      81195210,  78077248,  59128983,  63025520,  48129895,  51304566, 155270560,  59373566};

static const std::vector<std::string> CHR_NAME {"1", "2", "3", "4", "5", "6", "7", "8", "9",
      "10", "11", "12", "13", "14", "15", "16", "17", 
      "18", "19", "20", "21", "22", "X", "Y", "M"};

static SeqLib::HeaderSequenceVector hsv;
static SeqLib::BamHeader hdr;

static SeqLib::GRC points; 

namespace opt {

  std::string bias; // covered regions and their relative marginals
  uint32_t N = 1e6; // num to sample
  char mode = '1';
  double L = 1.001; // length distribution factor
  double R = 4;   // intra/inter ratio
  double variance = -1; 
};

static const char* shortopts = "hL:N:m:r:R:";
static const struct option longopts[] = {
  { "help",               no_argument, NULL, 'h' },
  { "length",             required_argument, NULL, 'L' },
  { "num-sims",           required_argument, NULL, 'N' },
  { "mode",               required_argument, NULL, 'm' },
  { "ratio",              required_argument, NULL, 'r' },
  { "random",              required_argument, NULL, 'R' },
  { NULL, 0, NULL, 0 }
};

void parseSimOptions(int argc, char** argv);

int runSim(int argc, char** argv) {

  parseSimOptions(argc, argv);

  // set up the header seq vec
  for (size_t i = 0; i < 23; ++i)
    hsv.push_back(SeqLib::HeaderSequence(CHR_NAME[i], csize[i]));

  // only sim from 1-X
  for (size_t i = 0; i < 23; ++i)
    GENOMESIZE += hsv[i].Length;

  // set up the random number generator
  std::default_random_engine generator;
  std::uniform_int_distribution<uint32_t> distribution(0,GENOMESIZE);

  // get a header file
  SeqLib::BamReader br;
  br.Open("/seq/picard_aggregation/G76270/NA12878/current/NA12878.bam");
  SeqLib::BamHeader h = br.Header();

  // start a timer
  std::chrono::time_point<std::chrono::system_clock> start, end, end2;

  // open the marginals
  Fractions fr;
  
  if (opt::variance == -1) {
    std::cerr << "...reading in " << opt::bias << std::endl;
    fr.readFromBed(opt::bias, h);
    std::cerr << "...read in " << SeqLib::AddCommas(fr.size()) << " regions " << std::endl;
  } else {
    std::cerr << "...generating randomed bias" << std::endl;
    fr = Fractions(10000, 0, hsv);
    double cfrac = 0.5; // current fraction
    for (auto& f : fr) {
      if (opt::variance == 0) { // uniform distribution
	f.frac = 1; 
      }
      else if (opt::variance == 1) { // totally random noise
	f.frac = ((double) (rand() % 10000)) / 10000; // random num between 0 1
      }
      else if (opt::variance == 2) { // slow random walk to get smooth curve
	double rfrac = ((double) (rand() % 10000)) / 10000; // random num between 0 1
	f.frac = cfrac + rfrac / 100 * (rand() % 2 ? -1 : 1);
	if (f.frac > 1)
	  f.frac = 1;
	else if (f.frac < 0)
	  f.frac = 0;
	cfrac = f.frac;
      } 
    }
    std::ofstream of("random_bias.bed");
    for (auto& f : fr)
      of << hsv[f.chr].Name << "\t" << f.pos1 << "\t" << f.pos2 << "\t" << f.frac << std::endl;
    of.close();
    std::cerr << "...generated bias of length " << SeqLib::AddCommas(fr.size()) << std::endl;    
  }
  

  // simulate points
  SeqLib::GRC grc;
  size_t initial_draw = 10000000;
  std::cerr << "...initially simulating " << SeqLib::AddCommas(initial_draw) << " breaks" << std::endl;
  for (size_t i = 0; i < initial_draw; ++i) 
    grc.add(regionFromNum(distribution(generator)));

  // do the overlaps
  std::cerr << "...creating interval tree on simulated points" << std::endl;
  start = std::chrono::system_clock::now();
  grc.CreateTreeMap();
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = end-start;
  std::cerr << "   interval tree construction time: " << elapsed.count() << "s\n";

  std::cerr << "...finding overlaps" << std::endl;
  start = std::chrono::system_clock::now();
  std::vector<int32_t> Query, Subject;
  fr.FindOverlaps(grc, Query, Subject, true); //good
  end = std::chrono::system_clock::now();
  elapsed = end-start;
  std::cerr << "   overlap time: " << elapsed.count() << "s\n";

  // loop through and subsample. Query is covered regions, subject is points
  std::cerr << "...sub-sampling " << SeqLib::AddCommas(Query.size()) << " simulated points to shape 1D marginals" << std::endl;
  start = std::chrono::system_clock::now();
  const int MOD = 100000;
  assert(Query.size() == Subject.size());
  for (size_t i = 0; i < Query.size(); ++i) {
    if (rand() % MOD < fr[Query[i]].frac * MOD) { // rand mod is [0,MOD], query*MOD is [0, MOD]
      int s = Subject[i];
      points.add(grc[s]);
    }
  }
  end = std::chrono::system_clock::now();
  elapsed = end - start;
  std::cerr << "   subsample time " << elapsed.count() << "s" << std::endl;

  start = std::chrono::system_clock::now();
  std::cerr << "...shuffling " << SeqLib::AddCommas(points.size()) << " points that survived marginal trimming" << std::endl;
  points.Shuffle();
  end = std::chrono::system_clock::now();
  elapsed = end - start;
  std::cerr << "   shuffle time " << elapsed.count() << "s" << std::endl;

  start = std::chrono::system_clock::now();
  if (opt::mode == '1') 
    output1D(points);
  else if (opt::mode == 'M')
    output2DMult(points);
  else if (opt::mode == 'A')
    output2DAdd(points);
  else 
    assert(false);
  end = std::chrono::system_clock::now();
  elapsed = end - start;
  std::cerr << "   point selection and output time " << elapsed.count() << "s" << std::endl;
    
  return 0;
}

void output1D(SeqLib::GRC& b) {

  std::cerr << "********* ONE DIMENSIONAL MODEL *********" << std::endl;  
  size_t c = 0;
  for (auto& i : b) {
    if (++c > opt::N)
      break;
    std::cout << hsv[i.chr].Name << "\t" << i.pos1 << "\t" << i.pos2 
	      << hsv[i.chr].Name << "\t" << i.pos1 << "\t" << i.pos2 << std::endl;
  }
  
  std::cerr << "...wrote " << SeqLib::AddCommas(b.size()) << " simulated 1D breaks" << std::endl;
  
}

void output2DMult(SeqLib::GRC& b) {
  
  const int FAIL_SAFE = 100;

  std::cerr << "********* MULTIPLICATIVE MODEL *********" << std::endl;

  size_t trans = 0;
  size_t intra = 0;

  // split into chromosomoes
  std::unordered_map<int, MRC> b_chr;
  size_t id =0;
  for (const auto& i : b) {
    MRegion rr;
    rr.chr = i.chr;
    rr.pos1 = i.pos1;
    rr.pos2 = i.pos2;
    rr.id = id++;
    b_chr[i.chr].add(rr);
  }
  assert(id == b.size());

  // interval tree them
  for (auto& i : b_chr)
    i.second.CreateTreeMap();

  // setup which ones are used
  std::vector<bool> used(id, false);
  size_t tried = 0;
  for (const auto& c : b_chr) { 
    for (size_t i = 0; i < c.second.size(); ++i) { 
      ++tried;
      const MRegion * draw = &c.second[i];

      if (trans + intra > opt::N)
	break;

      if (tried % 20000 == 0)
	std::cerr << "...drawing " << SeqLib::AddCommas(tried) << " of " << SeqLib::AddCommas(b.size()) << " trans " 
		  << SeqLib::AddCommas(trans) << " intra " << SeqLib::AddCommas(intra) << " ratio " << ((double)intra / trans) << std::endl;

      // if already used, skip
      if (used[draw->id]) {
	continue;
      }

      // reserve some for translocation
      double intra_ratio = (double)intra / trans;
      if (intra > 10 && intra_ratio > opt::R) {
	int p;
	bool pass;
	size_t fs = 0;
	do {
	  ++fs;
	  p = rand() % b.size();
	  pass = !used[p] && b[p].chr != draw->chr;
	} while(!pass && fs++ < FAIL_SAFE); // keep drawing if random partner is same chr or used
	
	if (!pass) 
	  continue;
	
	used[p] = true;
	++trans;
	const SeqLib::GenomicRegion * g = &b[p];
	std::cout << hsv[draw->chr].Name << "\t" << draw->pos1 << "\t" << draw->pos2 << "\t"
		  << hsv[g->chr].Name << "\t" << g->pos1 << "\t" << g->pos2 << std::endl;
	continue;
      }
      
      
      // draw a random distance as a cutoff threshold
      int p = drawFromPower(MIN, MAX, -opt::L);
      
      // find usable hits
      std::vector<int32_t> que, sub;
      MRegion padded = *draw;
      padded.Pad(p);
      MRC i_grc(padded);
      i_grc.CreateTreeMap();
      SeqLib::GRC out = i_grc.FindOverlaps(b_chr[draw->chr], que, sub, true);
      if (que.size() <= 1) { // should always have one (itself), but need more
	continue;
      }
      
      // of the hits, which are different from draw and not used
      std::vector<const MRegion*> avail;
      for (size_t z = 0; z < que.size(); ++z) {
	const MRegion * gg = &b_chr[draw->chr][sub[z]];
	if (gg != draw && !used[gg->id] && std::abs(gg->pos1 - draw->pos1) >= MIN)
	  avail.push_back(gg);
      }
      
      // if we have available hits, pick one. otherwise bail
      if (!avail.size()) {
	continue;
      }
      const MRegion * g = avail[rand() % avail.size()];
      used[g->id] = true;
      
      ++intra;
      std::cout << hsv[draw->chr].Name << "\t" << draw->pos1 << "\t" << draw->pos2 << "\t"
		<< hsv[g->chr].Name << "\t" << g->pos1 << "\t" << g->pos2 << std::endl;
      
    }

  }

  std::cerr << "...simulated " << SeqLib::AddCommas(intra) << " intrachr rearrangements and " << 
    SeqLib::AddCommas(trans) << " translocations" << std::endl;
  

  return;
}

void output2DAdd(SeqLib::GRC& b) {

  std::cerr << "********* ADDITIVE MODEL *********" << std::endl;  

  size_t intra = 0;
  size_t trans = 0;

  for (const auto& i : b) {
    
    if ((trans + intra) % 100000 == 0)
      std::cerr << "...drawing " << SeqLib::AddCommas(trans + intra) << " of " << SeqLib::AddCommas(b.size()) << " trans " 
		<< SeqLib::AddCommas(trans) << " intra " << SeqLib::AddCommas(intra) << " ratio " << ((double)intra / trans) << std::endl;


    if (trans + intra > opt::N)
      break;

    // randomly take out breaks to be used for translocations
    // actually need to sample at 2 * 1/ratio, bc we need to
    // remove two breaks for each break-pair. (for intra-chr we generate
    // the second break here
    if (rand() % PRECISION < (PRECISION / (1 + opt::R))) {
      SeqLib::GenomicRegion g;
      do {
	g = regionFromNum(rand() % GENOMESIZE);
      } while (g.chr == i.chr);
      std::cout << hsv[i.chr].Name << "\t" << i.pos1 << "\t" << i.pos2 << "\t"
		<< hsv[g.chr].Name << "\t" << g.pos1 << "\t" << g.pos2 << std::endl;
      ++trans;
      continue;
    }
      
    int32_t p = 0;
    
    do {
      p = i.pos1 + drawFromPower(MIN, MAX, -opt::L) * (rand() % 2 ? -1 : 1);
    } while (p < 0 || p > (int)hsv[i.chr].Length); // keep drawing until valid
   
    ++intra;
    std::cout << hsv[i.chr].Name << "\t" << i.pos1 << "\t" << i.pos2 << "\t"
	      << hsv[i.chr].Name << "\t" << p << "\t" << p << std::endl;
  }
    
  std::cerr << "...simulated " << SeqLib::AddCommas(intra) << " intrachr rearrangements and " << 
    SeqLib::AddCommas(trans) << " translocations" << std::endl;
}

void parseSimOptions(int argc, char** argv) {
  
  bool die = false;
  
  if (argc < 2) 
    die = true;
  else
    opt::bias = std::string(argv[1]);
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'N': arg >> opt::N; break;
    case 'L': arg >> opt::L; break;
    case 'r': arg >> opt::R; break;
    case 'm': arg >> opt::mode; break;
    case 'R': arg >> opt::variance; break;
    }
  }

  if (!(opt::mode == 'M' || opt::mode == '1' || opt::mode == 'A')) {
    std::cerr << "\nERROR: Simulation mode must be one of: M (multiplicative) A (additive) or 1 (one-dimensional)\n" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (die) {
    std::cerr << "\n" << SIM_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}

SeqLib::GenomicRegion regionFromNum(uint32_t number) {
  
  size_t j = 0; 
  for (;j < 24; ++j) 
    if (number - csum[j] <= csize[j])
      break;
  assert(j != 24);
  return SeqLib::GenomicRegion(j, number-csum[j], number-csum[j]);

  
}
