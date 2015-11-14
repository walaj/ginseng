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


#include "Matrix.h"
#include "swap.h"

#include <regex>
#include <sstream>
#include <random>
#include <cassert>
#include <fstream>

#include "SnowTools/gzstream.h"
#include "SnowTools/SnowUtils.h"

#include <prng_engine.hpp>
#include <set>

// define a mask so we can store two chr in one uint16_t
#define CHR1_MASK = 0xFF;
#define CHR2_MASK = 0xFF00;

// min distance to be valid. Gets rid of Sanger high FP rate at small events
#define SANGER_DIST_LIM 8000
// dont include VCFs files with more than this many non-comment lines
#define SANGER_PER_FILE_LIMIT 5000

#define SITMO_RNG 1

#define MIN_BIN_WIDTH 10000

#define ANIMATION 1

static const size_t INTER = 25;

/*Matrix::Matrix(S nr, S nc, S ne, size_t nb, size_t ns) : m_num_bins(nb), m_num_steps(ns) { 
  
  std::cout << "Creating random sparse matrix of dims " << nr << "x" << nc << " with " << ne << " events (1's)" <<  std::endl;

  // initialize the m_vec
  assert(m_vec.size() == 0);
  for (size_t i = 0; i < 26; i++)
    m_vec.push_back({});

  // fill with random values
  for (S i = 0; i < ne; i++) {
    if (i % 10000 == 0 && i > 0) 
      std::cout << "...generating event " << i << " of " << ne <<  std::endl;

    std::string key;
    size_t dcount = 0;
    GenomicRegion gr1, gr2;

    int dist;
    bool retry = true; // keep trying until get acceptable simluated event
    do {

      gr1.random();
      gr2.random();

      if (dcount++ > 1000) {
	std::cerr << "Weird loop on gen" << std::endl;
	exit(EXIT_FAILURE);
      }

      dist = gr1.distance(gr2);

      if ( (dist >= 0 && intraRatio() <= 4) || (dist >= 0 && (m_intra+m_inter) < 10)) { // it's intra-chrom and we need one
	// reject sample for 1/L re-creation
	S rv;
	SnowTools::genRandomValue(rv, 1000);
	
	assert(gr1.chr == gr2.chr);

	// reject some
	double one_tenth = nc / 50; // 10Mb - 1/10
	double tmp = (one_tenth/(double)dist)*1000;
	
	if (dist < one_tenth || tmp > rv)
	  retry = false;
	//if (dist > one_tenth && tmp < rv)  
	//  key = "x";
	//else
	//  key = gr1.toString() + "--" + gr2.toString(); 
      } else if ( (dist < 0 && intraRatio() > 4) || (dist < 0 && (m_intra+m_inter) < 10)) { // is inter and we need one
	retry = false;
	//key = gr1.toString() + "--" + gr2.toString(); 
      }

    } while (retry);
    
    if (gr1 < gr2)
      addMatrixValue(MatrixValue(gr1, gr2));
    else
      addMatrixValue(MatrixValue(gr2, gr1));
  }

  if ((m_intra+m_inter) == 0) {
    std::cerr << "Matrix is empty. Try increasing non-zero fraction (-f)" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (m_verbose)
    std::cout << "Matrix created with " <<  m_intra << " intra- and " << m_inter << " inter- events (ratio " << intraRatio() << ")" << std::endl;

  // fill event histograms
  if (m_verbose) 
    std::cout << "filling histogram, m_num_bins=" << m_num_bins << std::endl;

  fillQuantileHistograms(m_num_bins);
  std::cout << "Histogram created " << std::endl;

  m_orig = this;
}
*/

void Matrix::add(const Matrix& m) {

  std::set<std::string> hash;

  //std::cout << "ORIGINAL MATRIX SIZE IS " << m_vec.size() << std::endl;

  // hash the existing matrix
  for (auto& i : m_vec) 
    for (auto& j : i)
      hash.insert(std::to_string(j.c.chr) + ":" + std::to_string(j.c.pos1) + 
		  std::to_string(j.r.chr) + ":" + std::to_string(j.r.pos1));
  

  //std::cout << "HASHED MATRIX IS SIZE " << hash.size() << std::endl;

  // add the element if not there, otherwise update counter
  for (auto& i : m_vec) {
    for (auto& j : i) {
      std::string thash = std::to_string(j.c.chr) + ":" + std::to_string(j.c.pos1) + 
	std::to_string(j.r.chr) + ":" + std::to_string(j.r.pos1);
      if (hash.count(thash))
	++j.count;
      else 
	addMatrixValue(j);
    }
  }

  //std::cout << "NEW MATRIX SIZE IS " << m_vec.size() << std::endl;
}
  
void Matrix::fillQuantileHistograms(size_t num_bins) {

  std::vector<S>* pspanv = new std::vector<S>();

  // fill a vector of the spans
  getSpans(pspanv);

  // fill the quantile histogram
  hist.initialSpans(num_bins, pspanv, MIN_BIN_WIDTH);
  m_hist_smallbins.initialSpans(num_bins*10, pspanv, 0);
  //hist.toCSV("hist.csv");

  delete pspanv;

  // copy initial histogram to swap histogram
  hist_swap = hist;  

}

void Matrix::getSpans(std::vector<S>* pspanv) {

  pspanv->clear();
  for (auto& i : m_vec) 
    for (auto& j : i)
      pspanv->push_back(j.distance());

}

void Matrix::allSwaps() { //pthread_mutex_t * lock, std::vector<Matrix*> * allm) {
  
  if (!(m_intra+m_inter))
    return;

#ifdef ANIMATION
  // open the animation file if need to
  if (m_anim_step > 0 && m_id == 0) {
    of_anim.open(anim_file.c_str());
    of_anim_hist.open(anim_hist_file.c_str());
    of_anim_hist_small.open(anim_hist_small_file.c_str());
  }
#endif

#ifdef MV_LITE
    // hash the original matrix, for later comparison
    for (auto& i : m_vec) 
      for (auto& j : i)
        m_orig_map[std::to_string(j.c_chr) + ":" + std::to_string(j.c) + 
		   std::to_string(j.r_chr) + ":" + std::to_string(j.r)] = true;
#else
    // hash the original matrix, for later comparison
    for (auto& i : m_vec) 
      for (auto& j : i)
	m_orig_map[j.c.toString() + j.r.toString()] = true;
#endif

  // make the random values
  generateRandomVals();

  // do the swaps
  for (S i = 0; i < m_num_steps; i++) {
    doSwap();
    //if (m_verbose > 1/* && i % 100 == 0*/) {
    //  std::cout << "...step " << i << " of " << m_num_steps << std::endl;
    //}
  }
  // save the matrix
  //writeBinary();
  
  // add to the final list
  //pthread_mutex_lock(lock);  

  //dummy
  //allm->push_back(this);

  // print it out if need be
  //size_t mcount = allm->size();
  //if ((mcount % 1 == 0 && m_verbose) || (mcount == 1 && m_verbose)) {
  std::cout << "Swapped " << m_id << " matrices -- ";
  SnowTools::displayRuntime(start);
  std::cout << std::endl;
  std::cout << printMCMC() << std::endl;

#ifdef ANIMATION
  if (m_anim_step > 0 && m_id == 0) {
    of_anim.close();
    of_anim_hist.close();
    of_anim_hist_small.close();
  }
#endif

  //pthread_mutex_unlock(lock);
  
}

void Matrix::doSwap() {

  // new and faster
  size_t chr = INTER; 

  if (!inter_only)
    chr = rand_chr[m_mcmc.swap_tried];
  size_t i1 = rand_rows[m_mcmc.swap_tried] % m_vec[chr].size();
  size_t i2 = rand_cols[m_mcmc.swap_tried] % m_vec[chr].size();
 
  MatrixValue mo1 = m_vec[chr][i1];
  MatrixValue mo2 = m_vec[chr][i2];

  //std::cerr << "chr " << chr << std::endl;
  // make the swapped vals
  MatrixValue ms1, ms2;
#ifdef MV_LITE
  if (chr != INTER)
    MatrixValue::swapIntra(mo1, mo2, ms1, ms2);
  else
    MatrixValue::swapInter(mo1, mo2, ms1, ms2);
#else
  MatrixValue::swapIntra(mo1, mo2, ms1, ms2);
#endif

  S d3 = ms1.distance();
  S d4 = ms2.distance();

  // inter switching to intra. NOT VALUD
  if (chr == INTER && (d3 != INTERCHR || d4 != INTERCHR)) {
    ++m_mcmc.swap_tried;
    return;
  }

  bool valid = (d3 >= SANGER_DIST_LIM && d4 >= SANGER_DIST_LIM);
  if (!valid) {
    ++m_mcmc.swap_tried;
    return;
  }
  S d1 = mo1.distance();
  S d2 = mo2.distance();

  // positive energy shift is move AWAY from optimal
  int es = 0;
  if (chr != INTER) 
    es = energyShift(d1, d2, d3, d4); 

  valid = (es <= 0 || probs[0][m_mcmc.swap_tried] == TRAND) ? true: false;
  if (!valid && probs[0][m_mcmc.swap_tried] > 10) { // must be more than 10/65535 to even try
    assert(es > 0 && es <= 4);
    size_t randt = rand_cols[m_mcmc.swap_tried] & TRAND;
    // proceed if probability greater than random tval number
    if (probs[es-1][m_mcmc.swap_tried] > randt)
      valid = true;
  }

  ++m_mcmc.swap_tried;  

  if (valid) {
    
    //assert(temps[m_mcmc.swap_tried-1] >= 0.001 || es <= 0);
    ++m_mcmc.accepted;
    
    // update the swap histogram
    if (chr != INTER) {

      ++hist_swap.m_bins[bin_table[d3]];
      ++hist_swap.m_bins[bin_table[d4]];
      --hist_swap.m_bins[bin_table[d1]];
      --hist_swap.m_bins[bin_table[d2]];
      
#ifdef ANIMATION
      if (id == 0 && m_anim_step > 0) {
	m_hist_smallbins.addElem(d3);
	m_hist_smallbins.addElem(d4);
	m_hist_smallbins.removeElem(d1);
	m_hist_smallbins.removeElem(d2);
      }
#endif
    }

    // add the new ones
    m_vec[chr][i1] = ms1;
    m_vec[chr][i2] = ms2;
  }
  
#ifdef ANIMATION
  // output the animation
  if (m_anim_step > 0 && m_id == 0)
    if ( (m_mcmc.swap_tried % m_anim_step) == 0 || m_mcmc.swap_tried == 1) {
      std::cout << "...writing animation for step " << m_mcmc.swap_tried << " matrix id " << id << std::endl;
      this->toCSV(of_anim, of_anim_hist, of_anim_hist_small, m_mcmc.swap_tried);
    }
#endif
}


void Matrix::toSimpleCSV(const std::string& file) {

  std::ofstream fs;
  fs.open(file);
  
  for (auto& i : m_vec) 
    for (auto& j : i)
      fs << j.r.chr << "," << j.r.pos1 << "," << j.c.chr << "," << j.c.pos1 << "," << j.count << std::endl;
}

void Matrix::toCSV(std::ofstream &fs, std::ofstream &fh, std::ofstream &fh_small, size_t step) {
  
  // get the number of shared object
  //int shared_count = -1;
  //if (m_orig)
  //  shared_count = shared();

  for (auto& i : m_vec) 
    for (auto& j : i)
      fs << j.r.chr << "," << j.r.pos1 << "," << j.c.chr << "," << j.c.pos1 << "," << j.count << "," << step << std::endl;

  m_mcmc.old_accepted = m_mcmc.accepted;
  
  // print out the histogram
  hist_swap.toCSV(fh);

  // print out the small histogram
  m_hist_smallbins.toCSV(fh_small);

}

void Matrix::writeGzip(ogzstream * out) const {

  char sep = '\t';
  for (auto& i : m_vec) {
    for (auto& j : i) 
#ifdef MV_LITE
      (*out) << j.r_chr << sep << j.r << sep 
			 << j.c_chr << sep << j.c << sep
			 << m_id << std::endl;
#else
      (*out) << j.r.chr << sep << j.r.pos1 << sep 
			 << j.c.chr << sep << j.c.pos1 << sep
			 << m_id << std::endl;
#endif
  }
}

void Matrix::writeBinary() const {
  /*
  FILE* binout;
  std::string name = "matrix" + std::to_string(m_id) + ".bin";
  binout = fopen(name.c_str(),"wb");
  for (auto& i : m_vec) {
    for (auto& j : i) {
      uint32_t buffer[] = {j.r.chr, j.r.pos1, j.c.chr, j.c.pos1};
      fwrite(buffer, sizeof(uint32_t), sizeof(buffer), binout);
    }
  }
  fclose(binout);
  */
}

/*std::ostream& operator<<(std::ostream& out, const Matrix &m) {
  
  for (auto& i : m.m) 
    out << i << std::endl;
  return (out);
  }*/

int Matrix::energyShift(S odist1, S odist2, S pdist1, S pdist2) {

  int shift = 0;

  // calculate energy shift
  int weight1 = abs((int)hist_swap.binCount(bin_table[odist1]) - (int)hist.binCount(bin_table[odist1]));
  shift += (hist_swap.binCount(bin_table[odist1]) >  hist.binCount(bin_table[odist1]) ? -1 : 1)*weight1;
  int weight2 = abs((int)hist_swap.binCount(bin_table[odist2]) - (int)hist.binCount(bin_table[odist2]));
  shift += (hist_swap.binCount(bin_table[odist2]) >  hist.binCount(bin_table[odist2]) ? -1 : 1)*weight2;
  int weight3 = abs((int)hist_swap.binCount(bin_table[pdist1]) - (int)hist.binCount(bin_table[pdist1]));
  shift += (hist_swap.binCount(bin_table[pdist1]) >= hist.binCount(bin_table[pdist1]) ? 1 : -1)*weight3;
  int weight4 = abs((int)hist_swap.binCount(bin_table[pdist2]) - (int)hist.binCount(bin_table[pdist2]));
  shift += (hist_swap.binCount(bin_table[pdist2]) >= hist.binCount(bin_table[pdist2]) ? 1 : -1)*weight4;

  if (shift < -4)
    shift = -4;
  else if (shift > 4)
    shift = 4;
      
  return shift;

}

Matrix::Matrix(const std::string &file_list, size_t nb, size_t ns, SnowTools::GRC &mk, bool inter_only) : m_num_bins(nb), m_num_steps(ns) {

  this->inter_only = inter_only;
  
  //open the file
  igzstream inFile(file_list.c_str());
  if (!inFile) {
    std::cerr << "Can't read file " << file_list << " for parsing VCF" << std::endl;
    return;
  }

  // get count of files
  size_t num_files = SnowTools::countLines(file_list, "#", "vcf");

  // counter for number of masked events
  size_t masked = 0;

  // initialize the m_vec
  assert(m_vec.size() == 0);
  for (size_t i = 0; i < 26; i++)
    m_vec.push_back({});

  // loop through the files and add the events
  std::string file;
  size_t file_counter = 0;
  while (std::getline(inFile, file)) {
    file_counter++;
    size_t new_events = 0;
    if (!SnowTools::read_access_test(file)) {
      std::cerr << "VCF File: "  << file << " not readable or does not exist" << std::endl;
    } else {
      std::string val;
      igzstream this_file(file.c_str());
      std::string event_line;
      
      size_t num_events_in_file = SnowTools::countLines(file,"#");
      if (num_events_in_file < SANGER_PER_FILE_LIMIT) // block if too many lines
      while (std::getline(this_file, event_line)) {
	if (event_line.length() > 0 && event_line.at(0) != '#') {
	  size_t count = 0;
	  std::string chr1 = "-1", pos1 = "-1", chr2 = "-1", pos2 = "-1";
	  std::istringstream iss(event_line);

	  // regex to get mate information
	  while (std::getline(iss, val, '\t')) {
	    switch(++count) {
	    case 1: chr1 = val; break;
	    case 2: pos1 = val; break;
	    case 5: 
	      std::regex reg(".*?(\\]|\\[)([0-9A-Z]+):([0-9]+).*?$");	  
	      std::smatch match;
	      if (std::regex_search(val, match, reg) ) {
		chr2 = match[2];
		pos2 = match[3];
	      } else {
		std::cerr << "No match on line "  << val << std::endl;
	      }
	      break;
	    } // end switch
	  } // end intra-line loop

	  // check that we have the data
	  if (chr1 == "-1" || pos1 == "-1" || chr2 == "-1" || pos2 == "-1")
	    std::cerr << "Failed to parse VCF on line " << event_line << std::endl;

	  // convert to numbers
	  try { 

	    MatrixValue mv(chr1, pos1, chr2, pos2);
	    
	    // check the mask
	    bool keep = true;
	    if (mk.size()) {
	      size_t ovl = 0;
#ifdef MV_LITE
	      ovl += mk.findOverlapping(GenomicRegion(mv.r_chr, mv.r, mv.r));	      
#else
	      ovl += mk.findOverlapping(mv.r);
#endif
	      if (!ovl) {
#ifdef MV_LITE
		ovl += mk.findOverlapping(GenomicRegion(mv.c_chr, mv.c, mv.c));	      
#else
		ovl += mk.findOverlapping(mv.c);
#endif
	      }
	      if (ovl) {
		keep = false;
		++masked;
	      }
	    }
	    
	    // sanger conditional on distance
	    if ( mv.distance() >= SANGER_DIST_LIM && keep && (mv.r < mv.c) && (!inter_only || mv.r.chr != mv.c.chr)) {

	      // add the event
	      addMatrixValue(mv);
	      
	      // increment the event counter
	      new_events++;
	    }

	  } catch(...) {
	    std::cerr << "Error converting VCF line from std::string to number. Line " << event_line << std::endl;
	    std::cerr << "chr1 " << chr1 << " chr2 " << chr2 << std::endl;
	  }
	} // end ## conditional
	
      } // end intra-file loop
      else 
	std::cerr << "!!! Removed file " << file << " with too many events: " << num_events_in_file<< std::endl;
      
      // remove events if too many
      //if (new_events > MAX_EVENTS && false) {
      //	std::cout << "!!! Too many events: " << new_events << " -- ignoring all in this file !!!" << std::endl;
      //	m.erase(m.end() - new_events, m.end());
      //}
      
    } // end file OK conditional

    if (file_counter % 100 == 0)
      std::cout << "...imported VCF file " << file_counter << " of " << num_files << " with " << new_events << " events " << std::endl;

  } // end VCF file list

  if (m_verbose) 
    std::cout << "Imported " << (m_inter+m_intra) << " breakpoints with " << masked << " masked events " << std::endl;

  if ((m_intra+m_inter) == 0) {
    std::cerr << "ERROR: Didn't import any events" << std::endl;
    exit(EXIT_FAILURE);
  }

  //dedupe();

  // scramble to make even around the diagonal
  srand(1);
  for (auto& i : m_vec) 
    for (auto& j : i) 
      if (rand() % 2)
#ifdef MV_LITE
	j = MatrixValue(j.c_chr, j.c, j.r_chr, j.r);
#else	
	j = MatrixValue(j.c, j.r);
#endif

  // fill event histograms
  fillQuantileHistograms(m_num_bins);

  m_orig = this;

}


/*S Matrix::energy() {
  S en = 0;
  for (size_t i = 0; i < hist.numBins(); i++)
    en += abs(hist_swap.binCount(i) - hist.binCount(i));
  return en;
}*/

//std::string MatrixValue::key() const {
//  return (std::to_string(r.chr) + ":" + std::to_string(r.pos1) + "-" + std::to_string(c.chr) + ":" + std::to_string(c.pos1));
//}

/*std::ostream& operator <<(std::ostream &out, const MatrixValue &m) { 
  out << m.key() << " - " << m.distance();
  return out;
  }*/

void MatrixValue::swapIntra(const MatrixValue &m1, const MatrixValue &m2, MatrixValue &n1, MatrixValue &n2) {

  n1.r = m1.r;
  n1.c = m2.c;
  n2.r = m2.r;
  n2.c = m1.c;

}

void MatrixValue::swapInter(const MatrixValue &m1, const MatrixValue &m2, MatrixValue &n1, MatrixValue &n2) {

  n1.r = m1.r;
  n1.c = m2.c;
  n2.r = m2.r;
  n2.c = m1.c;

#ifdef MV_LITE
  n1.r_chr = m1.r_chr;
  n1.c_chr = m2.c_chr;
  n2.r_chr = m2.r_chr;
  n2.c_chr = m1.c_chr;
#endif

}

MatrixValue::MatrixValue(int chr1, int pos1, int chr2, int pos2) {
  r = GenomicRegion(chr1, pos1, pos1);
  c = GenomicRegion(chr2, pos2, pos2);
}

MatrixValue::MatrixValue(const std::string &chr1, const std::string &pos1, const std::string &chr2, const std::string &pos2) {
  
#ifdef MV_LITE
  r_chr = GenomicRegion::chrToNumber(chr1);
  c_chr = GenomicRegion::chrToNumber(chr2);
  try {
    r = std::stoi(pos1);
    c = std::stoi(pos2);
  } catch (...) {
    std::cerr << "stoi failed on MatrixValue constructor on one of " << pos1 << " or " << pos2 << std::endl;
  }
#else
  r = GenomicRegion(chr1, pos1, pos1);
  c = GenomicRegion(chr2, pos2, pos2);
#endif
  
}

void Matrix::addMatrixValue(const MatrixValue &mv) {
  
  MatrixValue mr = mv;

#ifdef MV_LITE
  if ( (mv.c_chr < mv.r_chr) || ( (mv.c_chr == mv.r_chr)  && (mv.c < mv.r) ) )
    {
      mr.c = mv.r;
      mr.c_chr = mv.r_chr;
      mr.r = mv.c;
      mr.r_chr = mv.c_chr;
    }
  
  // intra
  if (mv.r_chr == mv.c_chr) {
    ++m_intra;
    m_vec[mr.r_chr].push_back(mr);
  } else {
    ++m_inter;
    m_vec[INTER].push_back(mr);
    assert(m_vec[INTER].size() == m_inter);
  }
  

#else

  // put in order with row smaller
  if (mv.c < mv.r) {
    mr.c = mv.r;
    mr.r = mv.c;
  }
    
  // intra chr
  if (mv.r.chr == mv.c.chr) {
    ++m_intra;
    m_vec[mr.r.chr].push_back(mr);

  // inter chr
  } else {
    ++m_inter;
    m_vec[INTER].push_back(mr);
    assert(m_vec[INTER].size() == m_inter);
  }
#endif
}

void Matrix::generateRandomVals() {

  //clock_t t = clock();

  // destructor is responsible for freeing
  rand_rows = (uint32_t*) malloc(m_num_steps * sizeof(uint32_t));
  rand_cols = (uint32_t*) malloc(m_num_steps * sizeof(uint32_t));
  //rand_tval = (uint16_t*) malloc(m_num_steps * sizeof(uint32_t));

#ifdef SITMO_RNG
  // sitmo
  sitmo::prng_engine eng1;
  sitmo::prng_engine eng2;
  eng1.seed((m_id+1)*1000);
  eng2.seed((m_id+1)*10000);
  size_t sumr = m_intra+m_inter;
  for (size_t i = 0; i < m_num_steps; ++i){
    //size_t r1 = eng1();
    rand_rows[i] = eng1() % sumr;
    rand_cols[i] = eng2() % sumr;
    //rand_tval[i] = r1 % TRAND;
  }
#else

  // fill the values
  //srand(m_id);
  unsigned short seed[3] = {155,0,155};
  for (size_t i = 0; i < m_num_steps; ++i){
    rand_rows[i] = nrand48(&seed[0]) % (m_intra+m_inter);
    rand_cols[i] = nrand48(&seed[0]) % (m_intra+m_inter);
    //rand_tval[i] = rand_rows[i] % TRAND;
  }
#endif

}

size_t Matrix::shared() {

  size_t count = 0;
  
#ifdef MV_LITE
  for (auto& i : m_vec)
    for (auto& j : i)
      count += m_orig_map.count(std::to_string(j.c_chr) + ":" + std::to_string(j.c) + 
			       std::to_string(j.r_chr) + ":" + std::to_string(j.r));
#else
  // run through swapped matrix and query hash from original
  for (auto& i : m_vec)
    for (auto& j : i)
      count += m_orig_map.count(j.c.toString() + j.r.toString());
#endif

  return count;

}

std::ostream& operator<<(std::ostream &out, const MCMCData &m) 
{
    out << "   Swap tried: " << m.swap_tried << " accepted " << m.accepted << " (" << SnowTools::percentCalc<size_t>(m.accepted, m.swap_tried) << "%)";
    return out;
}

std::string Matrix::printMCMC() const {
  std::stringstream ss;
  ss << m_mcmc;
  return ss.str();
}

void Matrix::dedupe() {

  size_t intra = 0;
  size_t inter = 0;
  //std::unordered_map<int, MVec> new_map;
  std::vector<MVec> new_vec;

  

  for (auto& i : m_vec) {

    std::sort(i.begin(), i.end());
    new_vec.push_back({});

    for (size_t j = 0; j < (i.size() - 1); j++)
      if (!(i[j].r == i[j+1].r) && !(i[j].c == i[j+1].c)) {
	new_vec.back().push_back(i[j]);
	//if (i == -1)
	//	  ++inter;
	//else
	//  ++intra;
      }

  }
  
  std::cout << "Events before dedupe: " << (m_intra+m_inter) << " after " << (inter+intra) << std::endl;
  m_vec = new_vec;
  m_intra = intra;
  m_inter = inter;
}

OverlapResult Matrix::checkOverlaps(SnowTools::GRC * grvA, SnowTools::GRC * grvB) {

  if (!grvA->size() || !grvB->size())
    return OverlapResult(0,m_intra+m_inter);

  size_t overlap = 0;
  size_t no_overlap = 0;

  // check overlaps

  for (auto& i : m_vec)
    for (auto& j : i) {
      if (j.r.chr != j.c.chr) {
	size_t ovl1 = grvA->findOverlapping(j.r);
	size_t ovl2 = grvB->findOverlapping(j.c);

	size_t ovl3 = grvA->findOverlapping(j.c);
	size_t ovl4 = grvB->findOverlapping(j.r);

	if ( (ovl1 > 0 && ovl2 > 0) || (ovl3 > 0 && ovl4 > 0) )
	  ++overlap;
	else
	  ++no_overlap;
      }
    }

  return OverlapResult(overlap, no_overlap);
 }
