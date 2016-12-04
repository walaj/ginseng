#ifndef GINSENG_SWAP_H__
#define GINSENG_SWAP_H__

#include <cstdlib>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <ostream>
#include <pthread.h>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>
#include <cassert>
#include <fstream>

#include "pthread-lite.h"
#include "Matrix.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "BEDIntervals.h"

typedef std::unordered_map<std::string, BEDIntervals> BEDMap;
typedef std::vector<SeqLib::GRC> BEDVec;

// Just a number to define how many sig-figs to go out on for randomness in temperature draw
// CANNOT BE BIGGER THAN 2^16-1 (using uint16_t)
#define TRAND 65535

// temperature to start at
#define TMAX 1000

class Matrix;

void import_bed_files(const std::string& bed_list, BEDMap& all_bed);
void parseMatrixOptions(int argc, char** argv);
void readStoredMatrices(const std::string& file);
int runSwap(int argc, char** argv);
void PrecalculateTemps(Matrix* m);
void PrecalculateHistogramBins(Matrix* m);
void PrerandomizeChromosomes(Matrix* m);
void printinfo();

static std::ofstream of_anim; 
static std::ofstream of_anim_hist; 
static std::ofstream of_anim_hist_small; 

// store results for all matrices that hit one thread
struct SwapThreadItem {

  SwapThreadItem(size_t i) : id(i) {}
 
  size_t id; // id of this thread

  std::stringstream inter_results; // results for points that span intervals
  std::stringstream intra_results; // results for points that are contained in one interval
};

/** Class to hold work items for swapping. 
 * 
 * Each work item holds a pointer to a Matrix object. 
 * Running this work item on a thread executes all of the 
 * swaps for that Matrix, changing the value of the matrix
 * object.
 */
class SwapWorkItem {

private:

  Matrix * m_orig_mat; // store pointer to orig data to copy to this swap run
  Matrix * m_this_mat; // the matrix for this work item
  size_t m_id;         // some unique id (passed in)
  
  const BEDMap * m_bed_mat; // list of bed file

public:

  /** Create a new swap run to do all steps on one matrix
   * @param mat Matrix to be swapped
   */
  SwapWorkItem(Matrix* mat, size_t id, 
	       const BEDMap *mB)
    : m_orig_mat(mat), m_id(id), m_bed_mat(mB) {}

  /** Kick off all of the swapping and annealing
   * @return true if allSwaps completed
   */
  bool run(SwapThreadItem* thread_data) { 

    m_this_mat = new Matrix(*m_orig_mat); // copy orig data to this mat
    m_this_mat->m_orig = m_orig_mat;       // let new one see original
    m_this_mat->id = m_id;                // set unique id
#ifdef DEBUG_SWAP
    pthread_mutex_lock(m_lock);  
    std::cerr << "...starting swaps on " << m_id << std::endl;
    pthread_mutex_unlock(m_lock);
#endif
    m_this_mat->allSwaps();               // do the actual swaps

    bool exclusive = true;
    thread_data->inter_results << m_this_mat->OutputOverlapsInter(exclusive);
    thread_data->intra_results << m_this_mat->OutputOverlapsIntraExclusive();
    
    /*
    // get the overlaps
    std::unordered_map<std::string, OverlapResult> all_overlaps;
    std::unordered_map<std::string, bool> ovl_ovl_seen;

    for (const auto& i : *m_bed_mat) {
      for (const auto& j : *m_bed_mat) {
	if (i.first < j.first || (!ovl_ovl_seen.count(i.first) && i.first==j.first) ) { // don't need to do each one twice
	  OverlapResult ovl = m_this_mat->checkOverlaps(&i.second, &j.second);
	  std::string ovl_name = i.first + "," + j.first;
	  all_overlaps[ovl_name] = ovl;
#ifdef DEBUG_SWAP
	  std::cerr << " Overlap for "  << ovl_name << " is "  << ovl.first << std::endl;      
#endif
	  thread_data->inter_results << ovl_name << "," << ovl.first << "," << ovl.second << "," << m_this_mat->id  << std::endl; 
	  if (i.first==j.first)
	    ovl_ovl_seen[i.first] = true;
	}
      }
    }

    // check intra unit overlaps
#ifdef DEBUG_SWAP
    std::cerr << "...checking intra-overlaps on " << m_id << std::endl;
#endif
    for (const auto& i : *m_bed_mat) {
      if (do_intra_overlap(i.first)) {
	OverlapResult ovl = m_this_mat->checkIntraUnitOverlaps(&i.second);
	thread_data->intra_results << i.first << "," << ovl.first << "," << ovl.second << "," << m_this_mat->id << std::endl; 
      }
    }
    */

    delete m_this_mat;
    return true; 
  }

};

#endif
