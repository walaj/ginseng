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
	//assert(r.chr != c.chr);
	return INTERCHR;
      }
      //assert(r.chr == c.chr);
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
  uint16_t count = 0;
  uint16_t id = 0;
};

typedef std::vector<MatrixValue> MVec;

/** Container for a sparse matrix 
 */
class Matrix {

  friend class SwapWorkItem;
  friend class Squares;

 public:

  MatrixValue ms1, ms2; // these are holders for when vals get swapped in. Try 
  // declaring here to see if this is faster than re-declaring in every loop in 
  // doSwap

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
  Matrix(size_t ne, size_t nb, size_t nsteps, double pl, double frac_inter, const std::string& tid);

  /** Make a matrix from a list of VCF files
   * @param file_list Text file containing list of VCFs
   * @param nb Number of bins to divide distance histogram
   * @param ns Number of steps (swaps) to attempt
   * @param mk Mask regions such that events in this region are removed
   * @param inter_only only load interchromosomal events and only do inter chr
   */
  Matrix(const std::string &file_list, size_t nb, size_t ns, 
	 SnowTools::GRC &mk, bool inter_only, const std::vector<std::string>& identifiers, const std::string& tid);
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

  void __initialize_mvec();

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

    std::string analysis_id;
    // map to store sparse matrix entries, per chrom
    std::vector<MVec> m_vec;

    SnowTools::Histogram m_hist_smallbins;
    
};

#endif
