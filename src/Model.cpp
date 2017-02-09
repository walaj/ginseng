#include "Model.h"

#include <iostream>
#include <random>

typedef double (*math_fn)(apop_data *);

gsl_vector *jack_iteration(gsl_matrix *m, math_fn do_math){
  int height = m->size1;
  gsl_vector *out = gsl_vector_alloc(height);
  apop_data *reduced = apop_data_alloc(0, (size_t)(height - 1), (int)m->size2);
  Apop_submatrix(m, 1, 0, height - 1, m->size2, mv);
  gsl_matrix_memcpy(reduced->matrix, mv);
  for (int i=0; i< height; i++){
    gsl_vector_set(out, i, do_math(reduced));
    if (i < height - 1){
      Apop_matrix_row(m, i, onerow);
      gsl_matrix_set_row(reduced->matrix, i, onerow);
    }
  }
  return out;
}

apop_data *ols_data;
gsl_vector * predicted;
double p_dot_mse;

double sum_squared_diff(gsl_vector *left, gsl_vector *right){
  gsl_vector_sub(left, right); //destroys the left vector
  return apop_vector_map_sum(left, gsl_pow_2);
}

gsl_vector *project(apop_data *d, apop_model *m){
  return apop_dot(d, m->parameters, 0, 'v')->vector;
}

double cook_math(apop_data *reduced){
  apop_model *r = apop_estimate(reduced, apop_ols);
  double out = sum_squared_diff(project(ols_data, r), predicted)/p_dot_mse;
  apop_model_free(r);
  return out;
}

gsl_vector *cooks_distance(apop_model *in){
  apop_data  *c = apop_data_copy(in->data);
  apop_ols->prep(in->data, in);
  ols_data = in->data;
  predicted = project(in->data, in);
  p_dot_mse  = c->matrix->size2 * sum_squared_diff(in->data->vector, predicted);
  return jack_iteration(c->matrix, cook_math);
}

void test_apo() {
  
  // apop_data_set(apop, row, col, val, *, *, *)

  /*
  apop_data *scalar = apop_data_alloc(1, 0, 0);
  apop_data_set(scalar, 0,0,3,NULL,NULL,NULL);
  double three = apop_data_get(scalar,0,0,NULL,NULL,NULL);
  std::cerr << " THRE " << three << std::endl;

  //vector only, so the functions take only one coordinate
  apop_data *v = apop_data_alloc(3,0,0); // vec size 3, matrix 0x0
  apop_data_set(v, 0,-1,21,NULL,NULL,NULL);  // row, col, val
  apop_data_set(v, 1,-1,33,NULL,NULL,NULL);
  apop_data_set(v, 2,-1,43,NULL,NULL,NULL); 
  double two = apop_data_get(v, 0, -1, NULL, NULL, NULL); // row, col
  double two2 = apop_data_get(v, 1, -1, NULL, NULL, NULL); // row, col
  double two3 = apop_data_get(v, 2, -1, NULL, NULL, NULL); // row, col
  std::cerr << " TWO " << two << " " << two2 << " " << two3 << std::endl;
  */
  
  // make a big ol matrix
  /*apop_data *mat = apop_data_alloc(10, 10, 2);
  for (size_t i = 0; i < 10; ++i) { // loop each variable
    for (int j = 0; j < 2; ++j) { // loop each observation
      double r = i*3;
      if (j == 1)
	r = 1;
      apop_data_set(mat, i, j, r, NULL, NULL, NULL);
    }
    }*/

  // allocate the outcome
  apop_data *mat = apop_data_alloc(33, 0, 0);
  // mu = 100, var = 1000, p = .1  r = 111.1111
  double nb[33] = {93,109,96,85,102,100,111,105,92,89,112,118,96,103,102,108,103,94,85,107,103,125,102,91,98,106,108,97,101,89,97,115,107};
  for (size_t i = 0; i < 33; ++i) 
    apop_data_set(mat, i, -1, nb[i], NULL, NULL, NULL);
  //apop_data_show(mat);

  //apop_model *est = apop_estimate(mat, apop_negativebinomial);
  //apop_model_show(est);

  const size_t nd = 100000;
  gsl_rng * rrr = gsl_rng_alloc (gsl_rng_taus);
  apop_data *mat2 = apop_data_alloc(nd, 0, 0);
  // mu = 100, var = 1000, p = .1  r = 111.1111
  double nb2[nd];
  double mu = 500;
  double var = 10000;
  double p = mu / var;
  double r = mu*mu / (var - mu);
  std::cout << "Formulation 1 -- mu: " << mu << " variance: " << var << std::endl 
	    << "Formulation 2 (Apop form): probability p: " << p 
	    << " draws r: " << r << std::endl 
	    << "... estimating from " << nd << " points " << std::endl;
  for (size_t i = 0; i < nd; ++i) {
    nb2[i] = gsl_ran_negative_binomial(rrr, p, r);
    apop_data_set(mat2, i, -1, nb2[i], NULL, NULL, NULL);
  }

  //apop_data_show(mat2);

  std::cerr << "...estimating model" << std::endl;
  apop_model *est2 = apop_estimate(mat2, apop_negativebinomial);
  apop_model_show(est2);
  

    //std::cerr << " POT " << std::endl;
    //strcpy(apop_opts.output_delimiter, "\n");
    //apop_vector_show(cooks_distance(est));

  //  apop_data_add_names(v, 'r', "top row", "mid row", "low row");
  //   double two_again = apop_data_get(v, .rowname="low row");

  /*
  //matrix, so functions take two coordinates
  apop_data *m = apop_data_alloc(12, 12);
  apop_data_set(m, 3, 2, .val = 3.2);
  double threetwo = apop_data_get(v, 3, 2);

  //given vector and matrix, the vector is column -1.
  apop_data *mv = apop_data_alloc(12, 12, 12);
  apop_data_set(m, 3, -1, .val = 3.2);
  double vthree = apop_data_get(v, 3, -1);

  std::cerr << "VTHREE " << vthree << std::endl;
  */
}
