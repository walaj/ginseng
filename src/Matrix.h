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

#ifndef SWAP_MATRIX_H__
#define SWAP_MATRIX_H__

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <fstream>

#include "SnowTools/Histogram.h"
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "Squares.h"

typedef uint32_t S;    
typedef std::pair<size_t, size_t> OverlapResult;

// Should we use lite MatrixValue code or GenomicRegion code?
//#define MV_LITE 1

using SnowTools::GenomicRegion;

/** Store basic information about the MCMC run
 */
struct MCMCData {

  size_t swap_tried = 0;
  size_t swap_tried_valid= 0;
  size_t old_accepted = 0;
  size_t accepted = 0;

  friend std::ostream& operator<<(std::ostream &out, const MCMCData &m);

};

/** Stores a single entry of a Matrix
 */
class MatrixValue {

  friend class Matrix;
  friend class Squares;
  
 public:
  
  MatrixValue() {}

  MatrixValue(int chr1, int pos1, int chr2, int pos2);

#ifndef MV_LITE
  MatrixValue(const GenomicRegion &gr1, const GenomicRegion &gr2) : r(gr1), c(gr2) {}
#endif
    /** Construct a new MatrixValue by inputting strings directly from VCF (eg. "chr1", "1" or "X", "chrX")
     */
  MatrixValue(const std::string &chr1, const std::string &pos1, const std::string &chr2, const std::string &pos2);

#ifdef MV_LITE
  MatrixValue(uint8_t rchr, uint32_t rpos, uint8_t cchr, uint32_t cpos) : r_chr(rchr), r(rpos), c_chr(cchr), c(cpos) {}
#endif

  /** Fill in a pair of new MatrixValue as the row/col swap of two other values
   * @param m1 Current MatrixValue 1 to get row from for n1 and col from for n2
   * @param m2 Current MatrixValue 2 to get row from for n2 and col from for n1
   * @param n1 New matrix value 1
   * @param n2 New matrix value 2
   */
  static void swapInter(const MatrixValue &m1, const MatrixValue &m2, MatrixValue &n1, MatrixValue &n2);
  static void swapIntra(const MatrixValue &m1, const MatrixValue &m2, MatrixValue &n1, MatrixValue &n2);

  /** Order MatrixValue objects by their distance (span)
   */
  bool operator< (const MatrixValue &m) const 
  {
    S dist1 = this->distance();
    S dist2 = m.distance();
    return (dist1 < dist2 || (dist1==dist2 && r < m.r));
  }
      
  /** Return a string key in the form row,col 
   */
  //std::string key() const;
  
  /** Print the matrix value for debuging in form: (row,col) - dist Intra: 1/0
   */
  friend std::ostream& operator<<(std::ostream &out, const MatrixValue &m);

  /** Return the distance between the row and column (eg span)
   */
  inline int32_t distance() const
    {
#ifdef MV_LITE
      if (r_chr == c_chr)
	return (r > c) ? (r-c) : (c-r);
      else 
	return INTERCHR;
#else
      int dist = r.distance(c);
      if (dist < 0) {
	assert(r.chr != c.chr);
	return INTERCHR;
      }
      assert(r.chr == c.chr);
      return dist;
#endif
    }

  /** Is this point an intra-chromosomal point?
   */
  inline bool isIntra() const { return distance() >= 0; }

  //private:
#ifdef MV_LITE
  uint8_t r_chr;
  int32_t r;
  uint8_t c_chr;
  int32_t c;
#else
  GenomicRegion r;
  GenomicRegion c;
#endif
  size_t count = 0;
};

typedef std::vector<MatrixValue> MVec;

/** Container for a sparse matrix 
 */
class Matrix {

  friend class SwapWorkItem;
  friend class Squares;

 public:

  // histograms
  SnowTools::Histogram hist; // original histogram
  SnowTools::Histogram hist_swap; // histogram after N swaps 

  size_t id = 0;

  /** Construct a Matrix with simlulated data.
   * @param nr Number of rows
   * @param nc Number of columns
   * @param ne Number of events
   * @param nb Number of bins for distance histogram
   * @param ns Number of steps (swaps) to attempt
   */
  //Matrix(S nr, S nc, S ne, size_t nb, size_t ns);

  /** Make a matrix from a list of VCF files
   * @param file_list Text file containing list of VCFs
   * @param nb Number of bins to divide distance histogram
   * @param ns Number of steps (swaps) to attempt
   * @param mk Mask regions such that events in this region are removed
   * @param inter_only only load interchromosomal events and only do inter chr
   */
  Matrix(const std::string &file_list, size_t nb, size_t ns, SnowTools::GRC &mk, bool inter_only);
  Matrix() {}

  /** Delete this Matrix and free memory.
   */
  ~Matrix() { 
    free(rand_rows); 
    free(rand_cols); 
    //free(rand_tval); 
  }
  
  void add(const Matrix& m);

  /** Add a new MatrixValue to this Matrix
   */
  void addMatrixValue(const MatrixValue &mv);

  /** Print out the MCMC statistics
   */
  std::string printMCMC() const;

  /** Attempt a single swap in the Matrix.
   *
   * Each swap is an attempt to mix-up the matrix while exactly 
   * preserving the row and column marginals. This function will 
   * also evaluate the current temperature of the Matrix and 
   * decide whether to keep swaps that move the length histogram
   * away from the desired distribution. 
   */
  void doSwap();

  /** Generate all of the random numbers for row and col swaps
   * 
   * Note that this function takes as seed the matrix ID, and 
   * mallocs an array of uint32_t and length = number of swaps.
   * The destructor is responsible for freeing the allocs.
   */
  void generateRandomVals();

  /** Run all of the swaps for this Matrix.
   * @param lock Mutex to lock threads while add to allm
   * @param allm Vector storing all Matrix objects across threads
   * @param num_steps Number of swaps to attempt
   */
  void allSwaps(); //pthread_mutex_t * lock, std::vector<Matrix*> * allm);

  /** Print brief information about this Matrix for debugging.
   */
  friend std::ostream& operator<<(std::ostream& out, const Matrix &m);

  /** Output the matrix to a CSV file
   *
   * The output matrix has format:
   * - row,col,step,temp,accept_percent,unmoved,shared
   * 
   * Note that if step or temp are not provided, will default to zero
   * @param fs File stream to output the matrix data to
   * @param fh File stream to output the histogram data to
   * @param step The step number (of attempted swaps) that this Matrix is on
   */
  void toCSV(std::ofstream &fs, std::ofstream &fh, std::ofstream &fh_small, size_t step = 0); 

  void toSimpleCSV(const std::string& file); 

  /** Get the spans of all events in the matrix
      @param pspanv pointer to an initialized std::vector of spans
   */
  void getSpans(std::vector<S>* pspanv);

  /** Return the total number of shared event pairs 
   * between this Matrix object at swap N and original data
   * @return Total number of elements that are identical between pr
   */
  size_t shared();
  
  /** Remove duplicate breakpoint pairs
   */
  void dedupe();

  /** Loop through all of the Matrix elements and fill the quantile histogram
   *
   * A quantile histogram is one where, upon filling, each histogram bin
   * contains the same number of events (i.e. some fraction of the total
   * number of events). This means that the bin widths have variable sizes.
   * @param num_bins Number of bins to distribute Matrix events into
   */
  void fillQuantileHistograms(size_t num_bins);

  /** Remove a distance from the swapped histogram object stored in this object.
   * @param span Distance to remove from the histogram
   */
  void removeFromHistogram(S span);

  /** Add a distance to the swapped histogram object stored in this object.
   * @param span Distance to add to the histogram
   */
  void addToHistogram(S span);

  /** Calculate the energy shift that would result from the proposed swap
   */
  int energyShift(S odist1, S odist2, S pdist1, S pdist2);

  /** Set the step size to take a snapshot for animation. 
   *
   * Note that by passing a positive value, this matrix will
   * output an animation
   */
  void setAnimationStep(size_t anim) {
    std::cout << "Setting Matrix " << id << " to output animation csv." << std::endl;
    m_anim_step = anim;
  }

  /** Get the Intra/Inter chromosomal ratio
   * @return intra/inter if inter > 0, 0 otherwise
   */
  double intraRatio() const {
    if (m_inter == 0)
      return 100;
    else
      return (double)m_intra/(double)m_inter;
  }

  /** Check for events that span two GenomicRegionVector objects
   */
  OverlapResult checkOverlaps(SnowTools::GRC* grvA, SnowTools::GRC * grvB);

  //double * temps; // pre-compute all probabilites at each temp, since same for every matrix
  uint16_t** probs; // an array holding the arrays of probabilites for each shift 
  uint8_t* bin_table;

  bool inter_only = false;
  
  size_t m_id = 0;

  /** Store the matrix in binary form
   *
   * Format is 4 uint32_t numbers per element, storing r.chr, r.pos, c.chr, c.pos
   */
  void writeBinary() const;

  void addPoint() const;
  
  void writeGzip(ogzstream * out) const;

  uint8_t* rand_chr;
  //private:

    size_t m_num_bins = 100; 
    size_t m_anim_step = 0;
    size_t m_verbose = 1;
    size_t m_num_steps = 0;
    
    size_t m_inter = 0;
    size_t m_intra = 0;

    uint32_t * rand_rows;
    uint32_t * rand_cols;
    //uint16_t * rand_tval; // random numbers for temperatue monte carlo
    
    Matrix * m_orig; // pointer to original data Matrix
    std::unordered_map<std::string, bool> m_orig_map;
    
    MCMCData m_mcmc;

    // map to store sparse matrix entries, per chrom
    std::vector<MVec> m_vec;

    SnowTools::Histogram m_hist_smallbins;
    
};

#endif
