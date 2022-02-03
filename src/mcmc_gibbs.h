#ifndef MCMC_GIBBS_H
#define MCMC_GIBBS_H

#include <RcppArmadillo.h>
#include "redist_types.h"
#include "make_swaps_helper.h"
#include "map_calc.h"

double add_constraint(const std::string& name, List constraints,
                      std::vector<int> districts, NumericVector &psi_vec,
                      std::function<double(List, int)> fn_constr);

double calc_gibbs_tgt(const subview_col<uword> &plan, int n_distr, int V,
                      std::vector<int> districts, NumericVector &psi_vec, const uvec &pop,
                      double parity, const Graph &g, List constraints);

#endif
