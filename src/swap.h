/*
 * By downloading the PROGRAM you agree to the following terms of use:
 *
 * BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
 * FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
 *
 * This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
 * WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
 * WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
 * NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
 *
 * 1. DEFINITIONS
 * 1.1"PROGRAM" shall mean copyright in the object code and source code known as MuTect and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.or/cancer/cga/mutect on the EFFECTIVE DATE.
 *
 * 2. LICENSE
 * 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
 * LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
 * The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
 * 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
 * 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
 *
 * 3. OWNERSHIP OF INTELLECTUAL PROPERTY
 * LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
 *
 * Copyright 2012 Broad Institute, Inc.
 * Notice of attribution:  The MuTect program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc. and is published at doi: 10.1038/nbt.2514.
 *
 * LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
 *
 * 4. INDEMNIFICATION
 * LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
 *
 * 5. NO REPRESENTATIONS OR WARRANTIES
 * THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
 * IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 *
 * 6. ASSIGNMENT
 * This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
 *
 * 7. MISCELLANEOUS
 * 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
 * 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
 * 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
 * 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
 * 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
 * 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
 * 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
 */


/*! \mainpage Swap 0.1
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

#ifndef SWAP_SWAP_H__
#define SWAP_SWAP_H__

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

#include "SnowTools/gzstream.h"

#include "Matrix.h"


// timer
#ifndef __APPLE__
static struct timespec start;
#endif

// Just a number to define how many sig-figs to go out on for randomness in temperature draw
// CANNOT BE BIGGER THAN 2^16-1 (using uint16_t)
#define TRAND 65535

// temperature to start at
#define TMAX 100

class Matrix;

void parseMatrixOptions(int argc, char** argv);
void readStoredMatrices(const std::string& file);

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
  std::unordered_map<std::string, SnowTools::GRC> *m_bed_mat;

  std::ofstream * m_results;
  std::ofstream * m_results2;
  ogzstream * m_oz_matrix;
  Matrix * summed_results;

public:

  /** Create a new swap run to do all steps on one matrix
   * @param mat Matrix to be swapped
   */
  SwapWorkItem(Matrix* mat, std::vector<Matrix*> *all_mats, pthread_mutex_t * lock, size_t id, 
	       std::unordered_map<std::string, SnowTools::GRC> *mB, std::ofstream * ro, std::ofstream * ro2, ogzstream * oz, 
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
    m_this_mat->m_id = m_id;

    pthread_mutex_lock(m_lock);  
    std::cerr << "...starting swaps on " << m_id << std::endl;
    pthread_mutex_unlock(m_lock);
    m_this_mat->allSwaps(); 

    // get the overlaps
    std::unordered_map<std::string, OverlapResult> all_overlaps;
    std::unordered_map<std::string, bool> ovl_ovl_seen;
    std::stringstream ss, ss2;
    std::cerr << "...checking overlaps on " << m_id << std::endl;
    for (auto& i : *m_bed_mat) {
      for (auto& j : *m_bed_mat) {
	if (i.first < j.first || (!ovl_ovl_seen.count(i.first) && i.first==j.first) ) { // don't need to do each one twice
	  OverlapResult ovl = m_this_mat->checkOverlaps(&i.second, &j.second);
	  std::string ovl_name = i.first + "," + j.first;
	  all_overlaps[ovl_name] = ovl;
	  //std::cout << " Overlap for "  << ovl_name << " is "  << ovl.first << std::endl;      
	  ss << ovl_name << "," << ovl.first << "," << ovl.second << "," << m_this_mat->m_id  << std::endl; 
	  if (i.first==j.first)
	    ovl_ovl_seen[i.first] = true;
	}
      }
    }

    // check intra unit overlaps
    std::cerr << "...checking intra-overlaps on " << m_id << std::endl;
    for (auto& i : *m_bed_mat) {
      if (do_intra_overlap(i.first)) {
	OverlapResult ovl = m_this_mat->checkIntraUnitOverlaps(&i.second);
	ss2 << i.first << "," << ovl.first << "," << ovl.second << "," << m_this_mat->m_id << std::endl; 
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
