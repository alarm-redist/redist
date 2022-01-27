///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for sw_mh_helper.cpp, called by sw_mh_alg.cpp
///////////////////////////////////////

#ifndef SW_MH_HELPER_H
#define SW_MH_HELPER_H

#include <RcppArmadillo.h>
#include "redist_types.h"
#include "make_swaps_helper.h"
#include "constraint_calc_helper.h"
#include "map_calc.h"
#include "mcmc_gibbs.h"

using namespace Rcpp;

Rcpp::NumericVector init_pop(Rcpp::NumericVector popvec, arma::vec cds);
Rcpp::List add_ties(Rcpp::List aList);
Rcpp::List cut_edges(Rcpp::List aList_con,
		     double eprob);
Rcpp::List bsearch_boundary(Rcpp::List aList,
			    arma::vec boundary);
int count_valid(Rcpp::List aList, Rcpp::List boundarypart, Rcpp::NumericVector cdvec);
int draw_p(int lambda);
Rcpp::List make_swaps(Rcpp::List boundary_cc,
		      Rcpp::List aList,
		      Rcpp::NumericVector cds_old,
		      Rcpp::NumericVector pop_vec,
		      Rcpp::NumericVector cd_pop_vec,
		      Rcpp::List constraints,
		      CharacterVector psi_names,
		      double minparity,
		      double maxparity,
		      double parity,
		      int p,
		      double eprob,
		      double beta,
		      const Graph &g);
int mh_decision(double mh_prob);
Rcpp::List changeBeta(arma::vec betavec,
		      double beta,
		      double constraint,
		      Rcpp::NumericVector weights,
		      int adjswap);

#endif
