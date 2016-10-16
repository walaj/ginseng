#include "FishModel.h"
#include <omp.h>
#include <sstream>

#include <chrono>
#include <thread>

typedef double (*math_fn)(apop_data *);

// produce matrix with one row missing, for each possible row
// Do some math on the new sub-matrix
// Ben Klemens - MWD 4.6
static gsl_vector *jack_iteration(gsl_matrix *m, math_fn do_math){
  int height = m->size1;
  gsl_vector *out = gsl_vector_alloc(height);
  apop_data *reduced = apop_data_alloc(0, (size_t)height - 1, (int)m->size2);
  Apop_submatrix(m, 1, 0, height - 1, m->size2, mv);
  gsl_matrix_memcpy(reduced->matrix, mv);
  for (int i=0; i< height; i++){
    if (i % 100 == 0)
      std::cerr << "...jacknife at " << SeqLib::AddCommas(i) << " of " << SeqLib::AddCommas(height) << std::endl;
    gsl_vector_set(out, i, do_math(reduced)); // returns scalar output of do_math
    if (i < height - 1){ // create a new submatrix with new row ommited
      Apop_matrix_row(m, i, onerow);
      gsl_matrix_set_row(reduced->matrix, i, onerow);
    }
  }
  return out;
}

static apop_data *ols_data;
static gsl_vector * predicted;
static double p_dot_mse;
static apop_model * testmodel; //debug

static double sum_squared_diff(gsl_vector *left, const gsl_vector *right){
  gsl_vector_sub(left, right); //destroys the left vector
  return apop_vector_map_sum(left, gsl_pow_2);
}

static gsl_vector *project(const apop_data *d, const apop_model *m){
  return apop_dot(d, m->parameters, 0, 'v')->vector;
}

static double cook_math(apop_data *reduced){
  //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  //apop_model *r = apop_estimate(reduced, apop_ols); //debug
  gsl_vector* new_predicted = project(ols_data, testmodel); // debug r
  double out = sum_squared_diff(new_predicted, predicted)/p_dot_mse;
  free(new_predicted);
  //apop_model_free(r);
  return out;
}

static gsl_vector *cooks_distance(apop_model *in){
  testmodel = in; // debug
  apop_data  *c = apop_data_copy(in->data);
  apop_ols->prep(in->data, in);
  ols_data = in->data;
  predicted = project(in->data, in);
  p_dot_mse  = c->matrix->size2 * sum_squared_diff(in->data->vector, predicted);
  return jack_iteration(c->matrix, cook_math);
}

void FishModel::SetNumThreads(size_t n) { 
  omp_set_num_threads(n);
}

void FishModel::AddTiles(const FishHookTiles& f) {

  size_t nrows = f.size();
  size_t ncols = f.NumCovariates();
  // one for data, one for intercept
  mat = apop_data_alloc(0, nrows, (int)ncols+2); //vec, nrow, ncol
 
  // add the observed counts
  size_t i = 0; // row id (tile)
  for (i = 0; i < nrows; ++i) {
    float scaled_events = f.at(i).events;
    //scaled_events = f.at(i).covered > 0 ? scaled_events / f.at(i).covered : 0;
    apop_data_set(mat, i, (int)(ncols+1), scaled_events, NULL, NULL, NULL);
  }

  // set the covariates
  i = 0; // row id (tile)
  int j = 1; // column id (covariate) (0 is obs counts)
  for (const auto& r : f) { // loop the bins
    if (i == nrows) break;//debug, to break early
    j = 1;
    for (const auto& c : r) { // loop the covariates
      apop_data_set(mat, i, j++, c.second, NULL, NULL, NULL);
    }
    ++i; // update row interator
  }

  // add the intercept term
  for (i = 0; i < nrows; ++i)
    apop_data_set(mat, i, 0, 1.0, NULL, NULL, NULL);

  // add column names
  apop_name_add(mat->names, "events", 'c');
  for (const auto& c : f[0])
    apop_name_add(mat->names, c.first.c_str(), 'c');
  apop_name_add(mat->names, "intercept", 'c');

  // add row names
  std::stringstream ss;
  for (const auto& r : f) {
    ss << r.ChrName(SeqLib::BamHeader()) << ":" << r.pos1 << "-" << r.pos2;
    apop_name_add(mat->names, ss.str().c_str(), 'r');
    ss.str(std::string());
  }

}

void FishModel::EstimateOLS() {

  if (!mat) {
    std::cerr << "Need to fill the data with AddTiles first" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // run the model
  ols = apop_estimate(mat, apop_ols);

}

gsl_vector* FishModel::CooksDistance(apop_model* model) {
  
  if (!model) {
    std::cerr << "Need to supply CooksDistance with non-NULL model" << std::endl;
    exit(EXIT_FAILURE);
  }
  
 return cooks_distance(model);

}
