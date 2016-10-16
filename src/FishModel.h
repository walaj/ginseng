#ifndef GINSENG_FISH_MODEL_H
#define GINSENG_FISH_MODEL_H

#include "apop.h"
#include "FishHookInterval.h"

class FishModel {

 public:

 FishModel() : mat(NULL), ols(NULL) {}

  // fill the data
  void AddTiles(const FishHookTiles& f);

  // fill the ols model
  void EstimateOLS();

  // return the Cook's distance for each tile
  gsl_vector* CooksDistance(apop_model* model);

  /** Set the number of threads to use throughout
   * @param n Number of threads
   */
  void SetNumThreads(size_t n);

  /** Return the linear OLS model as an apop_model
   */
  apop_model * GetOLS() const { return ols; }

  /** Return the data as an apop_data struct
   */
  apop_data * GetApopData() const { return mat; }
  
 private:

  apop_data *mat; 

  apop_model * ols;
  
};


#endif
