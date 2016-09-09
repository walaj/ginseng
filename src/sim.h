#ifndef GINSENG_SIM_H__
#define GINSENG_SIM_H__

#include <cmath>
#include <vector>
#include <cassert>

#include "SeqLib/GenomicRegionCollection.h"


int runSim(int argc, char** argv);
void output1D(SeqLib::GRC& b);
void output2DMult(SeqLib::GRC& b);
void output2DAdd(SeqLib::GRC& b);
std::vector<int> drawFromPower(double x0, double x1, double power, int n_draws);
SeqLib::GenomicRegion regionFromNum(uint32_t number);

inline int drawFromPower(double x0, double x1, double power) {

  assert(power != -1);

  const int PRECISION = 1e6;

  double r = (double)(rand() % PRECISION)/(double)PRECISION;
  double t1 = std::pow(x1, power+1);
  double t2 = std::pow(x0, power+1);
  double tsum = (t1-t2) * r + t2;
  return std::floor(std::pow(tsum, 1 / (power + 1)));

}

class MRegion : public SeqLib::GenomicRegion {

 public: 
  int id;
  
};

typedef SeqLib::GenomicRegionCollection<MRegion> MRC;

#endif
