#include "fishhook.h"

#include <getopt.h>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "FishHookInterval.h"
#include "Fractions.h"
#include "SeqLib/BamReader.h"

static SeqLib::BamHeader hdr;

static const char *FISHHOOK_USAGE_MESSAGE =
"Program: fishhook \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ] && Marcin Imielinski [ marcin.imielinski@gmail.com ]\n"
"Usage: fishhook breaks.bed -I <interval_track> -S <scored_track> \n\n"
"Commands:\n"
"  -I     Input a set of interval tracks (as a BED file)\n"
"  -S     Input a set of scored tracks (as a BED file)\n"
"  -w     Set the width of the tiles\n"
"  -s     Set the \"slop\" (overlap) of the tiles\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

namespace opt {

  std::vector<std::string> interval_files;
  std::vector<std::string> score_files;
  int width = 20000;
  int slop = 0;
  bool verbose = true;
}

static const char* shortopts = "hvI:S:w:s:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "interval",                required_argument, NULL, 'I' },
  { "scored",                  required_argument, NULL, 'S' },
  { "slop",                    required_argument, NULL, 's' },
  { "width",                   required_argument, NULL, 'w' },
  { NULL, 0, NULL, 0 }
};

int runFishhook(int argc, char** argv) {

  parseFishOptions(argc, argv);

  if (opt::verbose) {
    std::cerr << "FishHook Params: " << std::endl
	      << "\tWidth: " << SeqLib::AddCommas(opt::width) << std::endl
	      << "\tSlop: " << SeqLib::AddCommas(opt::slop) << std::endl
	      << "\tInterval Tracks: " << std::endl;
    for (auto& i : opt::interval_files)
      std::cerr << "\t-- " << i << std::endl;
    std::cerr << "\tScore Tracks: " << std::endl;
    for (auto& i : opt::score_files)
      std::cerr << "\t-- " << i << std::endl;
  }

  // read in the covariate tracks
  SeqHashMap<std::string, Fractions> intervals;
  SeqHashMap<std::string, Fractions> scores;
  
  // read a header for info
  SeqLib::BamReader rdr;
  rdr.Open("/seq/picard_aggregation/G76270/NA12878/current/NA12878.bam");
  hdr = rdr.Header();

  // read in the interval tracks
  for (auto& i : opt::interval_files)
    read_track(i, intervals);
  for (auto& i : opt::score_files)
    read_track(i, scores);
  
  // create the tiled regions
  FishHookTiles fish(opt::width, opt::slop, hdr.GetHeaderSequenceVector());
  if (opt::verbose)
    std::cerr << "...constructed " << SeqLib::AddCommas(fish.size()) << " fishhook intervals" << std::endl;
  fish.CreateTreeMap();
  if (opt::verbose)
    std::cerr << "\tcreated inteval tree map" << std::endl;
  
  
  // overlap the covariates with the tiles
  for (auto& i : intervals)
    fish.AddIntervalCovariate(i.first, i.second);
  for (auto& i : scores)
    fish.AddScoreCovariate(i.first, i.second);

  return 0;
}

void parseFishOptions(int argc, char** argv) {

  bool die = false;  
  bool help = false;

  if (argc <= 2) 
    die = true;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    case 'h': help = true; break;
    case 'I': arg >> tmp; opt::interval_files.push_back(tmp); break;
    case 'S': arg >> tmp; opt::score_files.push_back(tmp); break;
    case 'w': arg >> opt::width; break;
    case 's': arg >> opt::slop; break;
    }
  }  

  if (die || help) {
      std::cerr << "\n" << FISHHOOK_USAGE_MESSAGE;
      if (die) exit(EXIT_FAILURE);
      else exit(EXIT_SUCCESS);	
    }
}

void read_track(const std::string& track, SeqHashMap<std::string, Fractions>& frac) {
  
  if (opt::verbose) std::cerr << "...reading BED file: " << track << std::endl;
  frac[track].readFromBed(track, hdr);
  if (opt::verbose) std::cerr << "\tread in " << SeqLib::AddCommas(frac[track].size()) << " regions " << std::endl;
  if (opt::verbose) std::cerr << "\tcreating interval tree map" << std::endl;
  frac[track].CreateTreeMap();
  if (opt::verbose) std::cerr << "\tdone creating tree map" << std::endl;
  

}
