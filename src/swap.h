#ifndef GINSENG_SWAP_H__
#define GINSENG_SWAP_H__

#include <cstdlib>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <ostream>
#include <pthread.h>
#include <list>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>
#include <cassert>
#include <fstream>

#include "gzstream.h"
#include "Matrix.h"
#include "SeqLib/GenomicRegionCollection.h"

typedef std::unordered_map<std::string, SeqLib::GRC> BEDMap;
typedef std::vector<SeqLib::GRC> BEDVec;

// Just a number to define how many sig-figs to go out on for randomness in temperature draw
// CANNOT BE BIGGER THAN 2^16-1 (using uint16_t)
#define TRAND 65535

// temperature to start at
#define TMAX 100

class Matrix;

void import_bed_files(const std::string& bed_list, BEDMap& all_bed);
void parseMatrixOptions(int argc, char** argv);
void readStoredMatrices(const std::string& file);
int runSwap(int argc, char** argv);

static std::ofstream of_anim; 
static std::ofstream of_anim_hist; 
static std::ofstream of_anim_hist_small; 

template <typename T> class wqueue
{ 

  public:
  wqueue() {
    pthread_mutex_init(&m_mutex, NULL);
    pthread_cond_init(&m_condv, NULL);
  }

  ~wqueue() {
    pthread_mutex_destroy(&m_mutex);
    pthread_cond_destroy(&m_condv);
  }

  void add(T item) {
    pthread_mutex_lock(&m_mutex);
    m_queue.push_back(item);
    pthread_cond_signal(&m_condv);
    pthread_mutex_unlock(&m_mutex);
  }

  T remove() {
    pthread_mutex_lock(&m_mutex);
    while (m_queue.size() == 0) {
      pthread_cond_wait(&m_condv, &m_mutex);
    }
    T item = m_queue.front();
    m_queue.pop_front();
    pthread_mutex_unlock(&m_mutex);
    return item;
  }

  int size() {
    pthread_mutex_lock(&m_mutex);
    int size = m_queue.size();
    pthread_mutex_unlock(&m_mutex);
    return size;
  }

  std::list<T>   m_queue;
  pthread_mutex_t m_mutex;
  pthread_cond_t  m_condv;

};

class SwapThread {

  public:
  SwapThread() : m_tid(0), m_running(0), m_detached(0) {}
  
  ~SwapThread()
  {
    if (m_running == 1 && m_detached == 0) {
      pthread_detach(m_tid);
    }
    if (m_running == 1) {
      pthread_cancel(m_tid);
    }
  }

  static void* runThread(void* arg)
  {
    return ((SwapThread*)arg)->run();
  }
 
  int start()
  {
    int result = pthread_create(&m_tid, NULL, runThread, this);
    if (result == 0) {
      m_running = 1;
    }
    return result;
  }

  int join()
  {
    int result = -1;
    if (m_running == 1) {
      result = pthread_join(m_tid, NULL);
      if (result == 0) {
	m_detached = 1;
      }
    }
    return result;
  }
 
  int detach()
  {
    int result = -1;
    if (m_running == 1 && m_detached == 0) {
      result = pthread_detach(m_tid);
      if (result == 0) {
	m_detached = 1;
      }
    }
    return result;
  }

  pthread_t self() {
    return m_tid;
  }

  virtual void* run() = 0;
 
private:
  pthread_t  m_tid;
  int        m_running;
  int        m_detached;
};

template <class T>
  class ConsumerThread : public SwapThread {
 
public:

  ConsumerThread(wqueue<T*>& queue, bool verbose) : m_queue(queue), m_verbose(verbose) {}
 
  void* run() {
    // Remove 1 item at a time and process it. Blocks if no items are 
    // available to process.
    for (int i = 0;; i++) {
      //if (m_verbose)
	//printf("thread %lu, loop %d - waiting for item...\n", 
	//     (long unsigned int)self(), i);
      T* item = (T*)m_queue.remove();
      item->run();
      //m_output->push_back(item->output());
      delete item;
      if (m_queue.size() == 0)
        return NULL;
    }
    return NULL;
  }

  //std::vector<O*> getThreadOutput() {return m_output;}

 private: 
  wqueue<T*>& m_queue;
  bool m_verbose;

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

  Matrix * m_orig_mat;
  std::vector<Matrix*> * m_all_mats;
  pthread_mutex_t * m_lock;
  Matrix * m_this_mat;
  size_t m_id;
  SeqHashMap<std::string, SeqLib::GRC> *m_bed_mat;

  std::ofstream * m_results;
  std::ofstream * m_results2;
  ogzstream * m_oz_matrix;
  Matrix * summed_results;

public:

  /** Create a new swap run to do all steps on one matrix
   * @param mat Matrix to be swapped
   */
  SwapWorkItem(Matrix* mat, std::vector<Matrix*> *all_mats, pthread_mutex_t * lock, size_t id, 
	       std::unordered_map<std::string, SeqLib::GRC> *mB, std::ofstream * ro, std::ofstream * ro2, ogzstream * oz, 
	       Matrix* sr)  
    : m_orig_mat(mat), m_all_mats(all_mats), m_lock(lock), m_id(id), m_bed_mat(mB), m_results(ro), m_results2(ro2), 
    m_oz_matrix(oz), summed_results(sr)
  {
    
  }

  /// Destroy this work item
    ~SwapWorkItem() { }


  /** Kick off all of the swapping and annealing
   * @return true if allSwaps completed
   */
  bool run() 
  { 
    m_this_mat = new Matrix(*m_orig_mat);
    m_this_mat->id = m_id;

    pthread_mutex_lock(m_lock);  
    //std::cerr << "...starting swaps on " << m_id << std::endl;
    pthread_mutex_unlock(m_lock);
    m_this_mat->allSwaps(); 

    // get the overlaps
    std::unordered_map<std::string, OverlapResult> all_overlaps;
    std::unordered_map<std::string, bool> ovl_ovl_seen;
    std::stringstream ss, ss2;
    //std::cerr << "...checking overlaps on " << m_id << std::endl;
    for (auto& i : *m_bed_mat) {
      for (auto& j : *m_bed_mat) {
	if (i.first < j.first || (!ovl_ovl_seen.count(i.first) && i.first==j.first) ) { // don't need to do each one twice
	  OverlapResult ovl = m_this_mat->checkOverlaps(&i.second, &j.second);
	  std::string ovl_name = i.first + "," + j.first;
	  all_overlaps[ovl_name] = ovl;
	  //std::cout << " Overlap for "  << ovl_name << " is "  << ovl.first << std::endl;      
	  ss << ovl_name << "," << ovl.first << "," << ovl.second << "," << m_this_mat->id  << std::endl; 
	  if (i.first==j.first)
	    ovl_ovl_seen[i.first] = true;
	}
      }
    }

    // check intra unit overlaps
    //std::cerr << "...checking intra-overlaps on " << m_id << std::endl;
    for (auto& i : *m_bed_mat) {
      if (do_intra_overlap(i.first)) {
	OverlapResult ovl = m_this_mat->checkIntraUnitOverlaps(&i.second);
	ss2 << i.first << "," << ovl.first << "," << ovl.second << "," << m_this_mat->id << std::endl; 
      }
    }

    pthread_mutex_lock(m_lock);  
    //summed_results->add(*m_this_mat);
    (*m_results)  << ss.str();
    (*m_results2) << ss2.str();
    //m_this_mat->writeGzip(m_oz_matrix);
    pthread_mutex_unlock(m_lock);
    
    delete m_this_mat;
    return true; 
  }

};

#endif
