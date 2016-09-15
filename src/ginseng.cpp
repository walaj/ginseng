#include "sim.h"
#include "swap.h"
#include "fishhook.h"

//debug
#include "Model.h"

static const char *GINSENG_USAGE_MESSAGE =
"Program: ginseng \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ] and Marcin Imielinski [ mimielinski@nygenome.org ]\n"
"Usage: ginseng <command> [options]\n\n"
"Commands:\n"
"           sim            Simulate breakpoints (1D) and rearrangements (2D) on the genome\n"
"           swap           Run MCMC swap method to identify enriched 2D connections\n"
"           fishhook       Run Gamma-Poisson (fishhook) regression for 1D breaks using covariates\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

int main(int argc, char** argv) {

  //debug
  //test_apo();
  //return 0;

  if (argc <= 1) {
    std::cerr << GINSENG_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << GINSENG_USAGE_MESSAGE;
      return 0;
    } else if (command == "swap") {
      runSwap(argc -1, argv + 1);
    } else if (command == "fishhook") {
      runFishhook(argc-1, argv+1);
    } else if (command == "sim") {
      runSim(argc-1, argv+1);
    }
    else {
      std::cerr << GINSENG_USAGE_MESSAGE;
      return 0;
    }
  } 

}
