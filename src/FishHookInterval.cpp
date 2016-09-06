#include "FishHookInterval.h"


FishHookTiles::FishHookTiles(int width, int ovlp, const SeqLib::HeaderSequenceVector& h) {

  size_t chr = 0;
  for (auto& i : h) {
    size_t start = 0; 
    for (; start + width < i.Length; start += width) {
      add(FishHookInterval(chr, start, start + width));
    }
    ++chr;
  }

}

void FishHookTiles::AddIntervalCovariate(const std::string& name, const Fractions& f) {

}

void FishHookTiles::AddScoreCovariate(const std::string& name, const Fractions& f) {

}
