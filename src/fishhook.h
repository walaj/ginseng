#ifndef GINSENG_FISH_HOOK_H__
#define GINSENG_FISH_HOOK_H__

#include <string>
#include "Fractions.h"

void read_track(const std::string& track, SeqHashMap<std::string, Fractions>& frac);
void parseFishOptions(int argc, char** argv);
int runFishhook(int argc, char** argv);

#endif
