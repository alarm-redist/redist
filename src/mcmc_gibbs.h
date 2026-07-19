#ifndef MCMC_GIBBS_H
#define MCMC_GIBBS_H

#include "make_swaps_helper.h"
#include "map_calc.h"
#include "redist_types.h"
#include "scoring.h"
#include <RcppArmadillo.h>

double add_constraint(const std::string &name, Rcpp::List constraints, std::vector<int> districts,
                      Rcpp::NumericVector &psi_vec, std::function<double(Rcpp::List, int)> fn_constr);

double calc_gibbs_tgt(const arma::subview_col<arma::uword> &plan, int n_distr, int V,
                      std::vector<int> districts, Rcpp::NumericVector &psi_vec, const arma::uvec &pop,
                      double parity, const Graph &g, Rcpp::List constraints);

#endif
