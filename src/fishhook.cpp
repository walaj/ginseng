#include "fishhook.h"

#include <getopt.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "FishHookInterval.h"
#include "FishModel.h"
#include "Events.h"
#include "Fractions.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/RefGenome.h"

static SeqLib::BamHeader hdr;

static const char *FISHHOOK_USAGE_MESSAGE =
"Program: fishhook \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ] && Marcin Imielinski [ marcin.imielinski@gmail.com ]\n"
"Usage: fishhook breaks.bed -I <interval_track> -S <scored_track> \n\n"
"Commands:\n"
"  -b     BAM file to use just to get header to define genome\n"
"  -G     Reference genome\n"
"  -I     Covariate as BED, considering only intervals\n"
"  -S     Covariate as BED. 4th column is \"score\" value (eg replication time)\n"
"  -F     Sequence feature to count (eg GC). Checks both strands\n"
"  -w     Set the width of the tiles\n"
"  -s     Set the \"slop\" (overlap) of the tiles\n"
"  -M     Optional coverage mask to define applicable regions (eg high mapq regions)\n"
"  -p     Number of threads to use (multi-threading only available in some parts)\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

namespace opt {

  std::string events;
  std::string refgenome;
  std::string bam; // = "/seq/picard_aggregation/G76270/NA12878/current/NA12878.bam"; // a bam to get reference from
  std::string coverage;
  std::unordered_set<std::string> interval_files;
  std::unordered_set<std::string> scored_files;
  std::unordered_set<std::string> seq_features;
  int width = 20000;
  int slop = 0;
  bool verbose = false;
  size_t num_threads = 1;
}

static const char* shortopts = "hvI:w:s:b:M:G:F:S:p:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "interval",                required_argument, NULL, 'I' },
  { "refgenome",               required_argument, NULL, 'G' },
  { "scored",                  required_argument, NULL, 'S' },
  { "coverage",                required_argument, NULL, 'M' },
  { "bam",                     required_argument, NULL, 'b' },
  { "slop",                    required_argument, NULL, 's' },
  { "width",                   required_argument, NULL, 'w' },
  { "seq-feature",             required_argument, NULL, 'F' },
  { "num-threads",             required_argument, NULL, 'p' },
  { NULL, 0, NULL, 0 }
};

// https://www.safaribooksonline.com/library/view/c-cookbook/0596007612/ch10s15.html
static std::string getFileName(const std::string& s) {

  char sep = '/';

#ifdef _WIN32
  sep = '\\';
#endif

  size_t i = s.rfind(sep, s.length());
  if (i != std::string::npos) {
    return(s.substr(i+1, s.length() - i));
  }

  return(s);
}

int runFishhook(int argc, char** argv) {

  parseFishOptions(argc, argv);

  if (opt::verbose) {
    std::cerr << "FishHook Params: " << std::endl 
	      << "\tWidth: " << SeqLib::AddCommas(opt::width) << std::endl
	      << "\tEvents: " << opt::events << std::endl
	      << "\tCoverage Mask: " << opt::coverage << std::endl
	      << "\tSlop: " << SeqLib::AddCommas(opt::slop) << std::endl
	      << "\tInterval Tracks: " << std::endl;
    for (auto& i : opt::interval_files)
      std::cerr << "\t-- " << i << std::endl;
    std::cerr << "\tScored Tracks: " << std::endl;
    for (auto& i : opt::scored_files)
      std::cerr << "\t-- " << i << std::endl;
    std::cerr << "\tSequence Features: " << std::endl;
    for (auto& i : opt::seq_features)
      std::cerr << "\t-- " << i << std::endl;

  }

  // read in the covariate tracks
  SeqHashMap<std::string, Fractions> intervals;
  
  // read a header for info
  SeqLib::BamReader rdr;
  if (!rdr.Open(opt::bam)) {
    std::cerr << "Error: Could not read BAM supplied by -b: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }
  hdr = rdr.Header();
  
  // read in the reference genome
  SeqLib::RefGenome ref;
  if (!ref.LoadIndex(opt::refgenome)) {
      if (opt::seq_features.size()) {
	std::cerr << "Error: Could not read referene genome supplied by -G: " << opt::refgenome << std::endl;
	exit(EXIT_FAILURE);
      }
  }

  // read in the events
  if (opt::verbose) std::cerr << "...reading events " << opt::events << std::endl;
  EventList events;
  if (!events.readFromBed(opt::events, hdr)) {
    std::cerr << "Error: Could not read events BED: " << opt::events << std::endl;
    exit(EXIT_FAILURE);
  }
  if (opt::verbose) std::cerr << "...read in " << SeqLib::AddCommas(events.size()) << " events" << std::endl;
  events.CreateTreeMap();
  
  // create the tiled regions
  FishHookTiles fish(opt::width, opt::slop, hdr.GetHeaderSequenceVector());
  if (opt::verbose)
    std::cerr << "...constructed " << SeqLib::AddCommas(fish.size()) << " fishhook intervals" << std::endl;
  fish.CreateTreeMap();

  // read the coverage mask
  SeqLib::GRC cov;
  if (!opt::coverage.empty()) {
    if (opt::verbose) std::cerr << "...reading coverage mask " << opt::coverage << std::endl;
    cov.ReadBED(opt::coverage, hdr);
    if (opt::verbose) std::cerr << "...read in " << SeqLib::AddCommas(cov.size()) << " covered regions " << std::endl;
    if (!cov.size()) {
      std::cerr << "Non-empty coverage track read with 0 regions. Check that is non-empty BED" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (opt::verbose) std::cerr << "...creating interval tree map on covered regions and overlapping with tiles" << std::endl;
    cov.CreateTreeMap();

    // find covered amount per tile
    std::vector<int32_t> q, s;
    SeqLib::GRC ovlp;
    // fish is subject
    if (fish.size() > cov.size()) // do in most efficient order
      ovlp = cov.FindOverlaps(fish, q, s, false);
    else
      ovlp = fish.FindOverlaps(cov, s, q, false);
    if (opt::verbose) std::cerr << "..." << SeqLib::AddCommas(ovlp.size()) << " regions are covered" << std::endl;

    // set the amount covered by each
    for (size_t i = 0; i < ovlp.size(); ++i) {
      fish[s[i]].covered += (double)ovlp[i].Width() / fish[s[i]].Width();
    }

    // mask the events
    q.clear(); s.clear(); ovlp.clear();
    // events is subject
    if (events.size() > cov.size()) // do in most efficient order
      ovlp = cov.FindOverlaps(events, q, s, false);
    else
      ovlp = events.FindOverlaps(cov, s, q, false);

    EventList newe;
    // set the amount covered by each
    for (size_t i = 0; i < ovlp.size(); ++i) {
      newe.add(Event(ovlp[i], events.at(s[i]).id));
    }
    events = newe;
    events.CreateTreeMap();

    if (opt::verbose) std::cerr << "...kept " << SeqLib::AddCommas(events.size()) << " events after mask" << std::endl;
    
  } else {
    for (auto& i : fish)
      i.covered = 1; // the entire thing is covered if no mask provided
  }

  // read in the interval tracks
  for (auto& i : opt::interval_files)
    read_track(i, intervals, cov, false);
  for (auto& i : opt::scored_files)
    read_track(i, intervals, cov, true);

  // count events per tile (also de-dupes on patient per bin)
  fish.CountEvents(events);
  // overlap the covariates with the tiles
  for (auto& i : intervals) {
    fish.AddIntervalCovariate(i.first, i.second);
  }

  // make the matrix
  FishModel fm;
  fm.AddTiles(fish);
  
  fm.SetNumThreads(opt::num_threads);
  fm.EstimateOLS();
  fm.CooksDistance(fm.GetOLS());

  // write the covariates
  fish.PrintBEDHeader(std::cout);
  fish.PrintBED(std::cout, hdr);

  

  return 0;
}

void parseFishOptions(int argc, char** argv) {

  bool die = false;  
  bool help = false;

  if (argc < 2) 
    die = true;
  else
    opt::events = std::string(argv[1]);

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    case 'h': help = true; break;
    case 'I': 
      arg >> tmp; 
      opt::interval_files.insert(tmp); 
      break;
    case 'S': 
      arg >> tmp; 
      opt::scored_files.insert(tmp); 
      break;
    case 'F': 
      arg >> tmp; 
      opt::seq_features.insert(tmp); 
      break;
    case 'w': arg >> opt::width; break;
    case 'G': arg >> opt::refgenome; break;
    case 'p': arg >> opt::num_threads; break;
    case 'M': arg >> opt::coverage; break;
    case 'b': arg >> opt::bam; break;
    case 's': arg >> opt::slop; break;
    }
  }  

  if (die || help) {
      std::cerr << "\n" << FISHHOOK_USAGE_MESSAGE;
      if (die) exit(EXIT_FAILURE);
      else exit(EXIT_SUCCESS);	
    }
}

void read_track(const std::string& track, SeqHashMap<std::string, Fractions>& frac, const SeqLib::GRC& cov, 
		bool scored) {
  
  // check for dup file basenames
  std::unordered_set<std::string> fname;

  std::string track_basename = getFileName(track);

  // in case two files have same basename
  while (fname.count(track_basename)) 
    track_basename += "_1"; 
  fname.insert(track_basename);

  if (opt::verbose) std::cerr << "...reading BED file: " << track << std::endl;

  frac[track_basename].readFromBed(track, hdr, scored);
  Fractions * f = &frac[track_basename];
  if (opt::verbose) std::cerr << "\tcreating interval tree map" << std::endl;
  f->CreateTreeMap();

  if (cov.size()) { // should we mask the track

    if (opt::verbose) std::cerr << "\tmasking track" << std::endl;
    std::vector<int32_t> q, s;
    SeqLib::GRC ovlp;
    if (f->size() > cov.size()) // do most efficient overlap of track and coverage
      ovlp = cov.FindOverlaps<FracRegion>(*f, q, s, true);
    else 
      ovlp = f->FindOverlaps(cov, s, q, true);
    
    // rebuild track with only covered regions
    // f (track) is subject
    Fractions newf;
    for (size_t i = 0; i < ovlp.size(); ++i)
      newf.add(FracRegion(ovlp[i], f->at(s[i]).frac)); 
    frac[track_basename] = newf;
    f = &frac[track_basename];
    f->CreateTreeMap();
  }

  // restrict to covered

  if (opt::verbose) std::cerr << "\tread in " << SeqLib::AddCommas(f->size()) << " regions " << std::endl;

  if (opt::verbose) std::cerr << "\tdone creating tree map" << std::endl;
  

}

