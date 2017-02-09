#ifndef SWAP_MATRIX_H__
#define SWAP_MATRIX_H__

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <fstream>
#include <memory>

#include "gzstream.h"
#include "Histogram.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "Squares.h"
#include "BEDIntervals.h"

struct GifWriter;
typedef uint32_t S;    
typedef std::pair<size_t, size_t> OverlapResult;

class GR : public SeqLib::GenomicRegion {
 public:
 GR() : SeqLib::GenomicRegion() {}
 GR(const SeqLib::GenomicRegion& gr, int ii) : SeqLib::GenomicRegion(gr), id(ii) {}
  int id;
};

inline bool do_intra_overlap(const std::string& name) {

  std::vector<std::string> intra_unit_todo = {"TAD", "GENE", "FRAG","ALL", "LINE"};  
  for (auto& i : intra_unit_todo) {
    if (name.find(i) != std::string::npos) {
      return true;
    }
  }
  return false;
}

// Should we use lite MatrixValue code or GenomicRegion code?

using SeqLib::GenomicRegion;

/** Store basic information about the MCMC run
 */
struct MCMCData {

  size_t swap_tried = 0;
  size_t swap_tried_valid= 0;
  size_t old_accepted = 0;
  size_t accepted = 0;

  friend std::ostream& operator<<(std::ostream &out, const MCMCData &m);

};

// store the name of a BED track, and the element in the BED that the
// 1D point overlaps (-1 for none)
typedef std::unordered_map<std::string, int> BedAndHitIDMap;
typedef std::pair<std::string, int> BedAndHitID;

class MatrixPoint : public SeqLib::GenomicRegion {

 public: 

 MatrixPoint(int32_t c, int32_t p1, int32_t p2) : SeqLib::GenomicRegion(c, p1, p2) { strand = '*'; }

 MatrixPoint(const std::string& c, const std::string& p1, const std::string& p2, const SeqLib::BamHeader& h) : SeqLib::GenomicRegion(c, p1, p2, h) { strand = '*'; }

  // add a BED track and check if this point overlaps with it
  void AddBED(const std::string& b, BEDIntervals& bi) {
    std::vector<int> out = bi.grc.FindOverlappedIntervals(*this, true);
    if (!out.size()) {
      olap_element.insert(BedAndHitID(b, -1));
      return;
    }
    assert(out.size() == 1); // should not have overlapping intervals
    olap_element.insert(BedAndHitID(b, out[0]));
  }

  // for each track, store what element this overlaps
  BedAndHitIDMap olap_element; 

};

/** Stores a single entry of a Matrix
 */
class MatrixValue {

  friend class Matrix;
  friend class Squares;
  
 public:
  
  MatrixValue() {}

  MatrixValue(int chr1, int pos1, int chr2, int pos2);

  void Flip() { 
    std::swap(r,c);
  }

    /** Construct a new MatrixValue by inputting strings directly from VCF (eg. "chr1", "1" or "X", "chrX")
     */
  MatrixValue(const std::string &chr1, const std::string &pos1, const std::string &chr2, const std::string &pos2);

  /** Fill in a pair of new MatrixValue as the row/col swap of two other values
   * @param m1 Current MatrixValue 1 to get row from for n1 and col from for n2
   * @param m2 Current MatrixValue 2 to get row from for n2 and col from for n1
   * @param n1 New matrix value 1
   * @param n2 New matrix value 2
   */
  static void swap(const MatrixValue * m1, const MatrixValue * m2, MatrixValue &n1, MatrixValue &n2);

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
      if (r->chr != c->chr)
	return INTERCHR;
      return std::abs(r->pos1 - c->pos1);
    }

  /** Is this point an intra-chromosomal point?
   */
  inline bool isIntra() const { return distance() >= 0; }

  std::shared_ptr<MatrixPoint> r;
  std::shared_ptr<MatrixPoint> c;
  //uint16_t count = 0;

  //uint16_t id = 0;
};

typedef std::vector<MatrixValue> MVec;

/** Container for a sparse matrix 
 */
class Matrix {

  friend class SwapWorkItem;
  friend class Squares;

 public:

  /** Load a new BEDPE file into the matrix */
  bool LoadBEDPE(const std::string& file);

  /** Load a new VCF file in BEND format into the matrix */
  bool LoadVCF(const std::string& file);
  
  Matrix();

  /** Delete this Matrix and free memory. */
  ~Matrix(); 
  
  /** Add a new MatrixValue to this Matrix */
  void addMatrixValue(const MatrixValue &mv);

  /** Print out the MCMC statistics */
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

  void doTransSwap();
  void doIntraSwap(size_t chr);

  /** Check all of the 1D overlaps with this BED */
  void AddBedElements(const std::string& b, BEDIntervals& bi);

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
  bool allSwaps(); //pthread_mutex_t * lock, std::vector<Matrix*> * allm);

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
  size_t shared() const;
  
  /** Scramble the matrix around the diagonal to make non-symmetric */
  void ScrambleAroundDiagonal();

  /** Return the total number of rearrangements stored here */
  size_t size() const;

  /** Retrieve as a binned image 
   * step is the amount to increment red for each hit (up to 255)
   */
  uint8_t* AsRGB8Image(uint32_t width, int step) const;

  /** Retrieve as a binned image 
   * step is the amount to increment red for each hit (up to 255)
   */
  uint8_t* AsBMPImage(uint32_t width, int step) const;

  /** Loop through all of the Matrix elements and fill the quantile histogram
   *
   * A quantile histogram is one where, upon filling, each histogram bin
   * contains the same number of events (i.e. some fraction of the total
   * number of events). This means that the bin widths have variable sizes.
   * @param num_bins Number of bins to distribute Matrix events into
   */
  void FillQuantileHistograms(size_t num_bins);

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
    m_anim_step = anim;
  }

  /** Set the number of swaps to attempt per matrix */
  void SetNumSwapSteps(int v) { m_num_steps = v; }

  /** Set whether to run translocations only */
  void SetInterOnly(bool v) { m_inter_only = v; }

  /** Set whether to run translocations only */
  void SetIntraOnly(bool v) { m_intra_only = v; }

  /** Set the minimum allowed rearrangement distance */
  void SetMinSize(int s) { m_min_size = s >= 0 ? s : 0; }

  /** Set the maximum allowed rearrangement distance */
  void SetMaxSize(int s) { m_max_size = s >= 0 ? s : 0; }

  //double * temps; // pre-compute all probabilites at each temp, since same for every matrix
  uint16_t** probs = nullptr; // an array holding the arrays of probabilites for each shift 
  uint8_t* bin_table;

  /** Store the matrix in binary form
   *
   * Format is 4 uint32_t numbers per element, storing r.chr, r.pos, c.chr, c.pos
   */
  void writeBinary() const;

  /** Return the fraction of inter-chromosomal rearrangements */
  double GetFractionInterChromosomal() const {
    return (double)m_inter / (double)(m_inter + m_intra);
  }

  /** Return the number intra chromosomal rearrangements */
  size_t GetNumIntraChromosomal() const {
    return m_intra;
  }

  /** Write the matrix to a GZipped tsv file */
  void writeGzip(ogzstream * out) const;

  uint8_t* rand_chr = nullptr;
  
  // map to store sparse matrix entries, per chrom
  std::vector<MVec> m_vec;
  
  Histogram m_hist_smallbins;
  
  // histograms
  Histogram hist; // original histogram
  Histogram hist_swap; // histogram after N swaps 

  size_t id = 0; ///< ID for this matrix

  // at the end of swapping, output a formatted
  // string with the output results (for inter-interval exclusive)
  // if exclusive, then don't count events that are within the exact same
  // element (e.g. the exact same SINE element)
  std::string OutputOverlapsInter(bool exclusive) const;

  // at the end of swapping, output a formatted
  // string with the output results (for rars with both ends in 
  // the exact same interval)
  std::string OutputOverlapsIntraExclusive() const;

 private:

  std::vector<std::string> bed_list; // which bed files are we tracking

  Matrix * m_orig = nullptr; // true data matrix

  GifWriter * gw = nullptr; // for writing the animation

  std::unordered_map<std::string, bool> m_orig_map;

  size_t m_anim_step = 0;
  size_t m_num_steps = 0;
  
  size_t m_inter = 0; // hold a count of number of inter/intra rars
  size_t m_intra = 0;

  //std::shared_ptr<uint32_t> rand_rows;  // hold a bunch of random numbers
  //std::shared_ptr<uint32_t> rand_cols;
  uint32_t * rand_rows;  // hold a bunch of random numbers
  uint32_t * rand_cols;
  
  // store data about run
  MCMCData m_mcmc;

  // these are holders for when vals get swapped in
  MatrixValue ms1, ms2; 
  
  size_t m_num_bins = 100; 
  size_t m_min_size = 1000;
  size_t m_max_size = 250e6;
    
  // inter or intra chromosomal only
  bool m_inter_only = false;
  bool m_intra_only = false;

  // check that a file line is not header or weird
  bool ValidateLine(const std::string& l) const;

  // Add a new matrix value from strings
  bool AddNewMatrixValue(const std::string& c1, const std::string& c2,
			 const std::string& p1, const std::string& p2);

  void __initialize_mvec();

  void Animate();
    
};

#endif
