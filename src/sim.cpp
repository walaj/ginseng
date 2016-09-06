#include "sim.h"

#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <fstream>
#include <chrono>

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamReader.h"
#include "Fractions.h"

#define MAX_ENRICH 25

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: sim (simulate points on a genome)\n");
  //fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Jeremiah Wala <jwala@broadinstitute.org>\n\n");
  fprintf(stderr, "Usage:   sim covered.bed > out.bed\n\n");
  return -1;

}
int runSim(int argc, char** argv) {

  if (argc < 2) return usage();

  double max_enrich = -1;

  const std::string bias = std::string(argv[1]);

  const uint32_t GENOMESIZE = 3095693983;
  const size_t n = 1000000;
  const uint32_t csum[24] = {0, 249250621, 492449994, 690472424, 881626700, 1062541960,
                             1233657027, 1392795690, 1539159712, 1680373143,
                             1815907890, 1950914406, 
			     2084766301, 2199936179, 2307285719, 2409817111, 2500171864, 2581367074,
			     2659444322, 2718573305, 2781598825, 2829728720, 2881033286, 3036303846};

  const uint32_t csize[24] = {249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
			      141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392,  90354753,
			      81195210,  78077248,  59128983,  63025520,  48129895,  51304566, 155270560,  59373566};

  static const std::vector<std::string> CHR_NAME {"1", "2", "3", "4", "5", "6", "7", "8", "9",
      "10", "11", "12", "13", "14", "15", "16", "17", 
      "18", "19", "20", "21", "22", "X", "Y", "M"};

  SeqLib::HeaderSequenceVector hsv;
  for (size_t i = 0; i < 24; ++i)
    hsv.push_back(SeqLib::HeaderSequence(CHR_NAME[i], csize[i]));

  std::default_random_engine generator;
  std::uniform_int_distribution<uint32_t> distribution(0,GENOMESIZE);

  // get a header file
  SeqLib::BamReader br;
  br.Open("/seq/picard_aggregation/G76270/NA12878/current/NA12878.bam");
  SeqLib::BamHeader h = br.Header();

  // start a timer
  std::chrono::time_point<std::chrono::system_clock> start, end, end2;

  // open the thing
  std::cerr << "...reading in " << bias << std::endl;
  Fractions fr;
  fr.readFromBed(bias, h);
  std::cerr << "...read in " << SeqLib::AddCommas(fr.size()) << " regions " << std::endl;

  // find max enrichment
  for (auto& i : fr)
    if (i.frac > max_enrich) {
      max_enrich = i.frac;
      if (i.frac > MAX_ENRICH)
	i.frac = MAX_ENRICH;
    }
  std::cerr << "Max enrichment requested: " << max_enrich << std::endl;
  max_enrich = max_enrich > MAX_ENRICH ? MAX_ENRICH : max_enrich;
  
  // simulate 
  SeqLib::GRC grc;
  //uint32_t* p = (uint32_t*)malloc(n * sizeof(uint32_t));
  //uint8_t* c = (uint8_t*)malloc(n * sizeof(uint8_t));
  std::cerr << "...simulating" << std::endl;
  for (size_t i = 0; i < n * max_enrich; ++i) {
    uint32_t number = distribution(generator);

    size_t j = 0; 
    for (;j < 24; ++j) 
      if (number - csum[j] < csize[j])
	break;
    if (j == 24)
      continue;
    grc.add(SeqLib::GenomicRegion(j, number-csum[j], number-csum[j]));
  }
  std::cerr << "...simulated " << SeqLib::AddCommas(grc.size()) << " breaks " << std::endl;

  start = std::chrono::system_clock::now();
  // do the overlaps
  std::cerr << "...creating sim interval tree" << std::endl;
  grc.CreateTreeMap();
  end2 = std::chrono::system_clock::now();
  std::cerr << "...finding overlaps" << std::endl;
  std::vector<int32_t> Query, Subject;
  fr.FindOverlaps(grc, Query, Subject, true); //good

  end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds_create = end-start;
  std::chrono::duration<double> elapsed_seconds_query = end-end2;
  std::cerr << "tree construction time: " << elapsed_seconds_create.count() << "s\n";
  std::cerr << "tree query time: " << elapsed_seconds_query.count() << "s\n";
  std::cerr << "overlaps " << SeqLib::AddCommas(Query.size()) << std::endl;

  // loop through and subsample. Query is covered regions, subject is points
  std::cerr << "...writing" << std::endl;
  const int MOD = 100000;
  assert(Query.size() == Subject.size());
  size_t write_count = 0;
  for (size_t i = 0; i < Query.size(); ++i) {
    if (rand() % MOD < fr[Query[i]].frac / max_enrich * MOD) { // rand mod is [0,MOD], query*MOD is [0, MOD]
      int s = Subject[i];
      std::cout << CHR_NAME[grc[s].chr] << "\t" << grc[s].pos1 << "\t" << grc[s].pos2 << std::endl;
      ++write_count;
    }
  }

  std::cerr << "...wrote " << SeqLib::AddCommas(write_count) << " simulated breaks" << std::endl;
  return 0;
}
