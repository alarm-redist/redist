///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for constraint_calc_helper.cpp
///////////////////////////////////////

#ifndef CONSTRAINT_CALC_HELPER_H
#define CONSTRAINT_CALC_HELPER_H

#include "redist_types.h"

Rcpp::NumericVector findBoundary(Rcpp::List fullList,
				 Rcpp::List conList);
arma::uvec getIn(arma::ivec vec1, arma::ivec vec2);
Rcpp::List genAlConn(Rcpp::List aList,
		     Rcpp::NumericVector cds);
Rcpp::List calc_psipop(arma::vec current_dists,
		       arma::vec new_dists,
		       Rcpp::NumericVector pops,
		       Rcpp::NumericVector distswitch);
Rcpp::List calc_psicompact(arma::vec current_dists,
			   arma::vec new_dists,
			   Rcpp::NumericVector distswitch,
			   std::string measure,
			   Rcpp::List aList,
			   Rcpp::NumericVector areas_vec,
			   arma::mat borderlength_mat,
			   bool discrete,
			   Rcpp::NumericVector pops,
			   Rcpp::NumericMatrix ssdmat,
			   int ndists,
			   const Graph &g,
			   arma::vec counties,
			   double denominator
			   );
Rcpp::List calc_psisegregation(arma::vec current_dists,
                               arma::vec new_dists,
                               Rcpp::NumericVector pops,
                               Rcpp::NumericVector distswitch,
                               Rcpp::NumericVector grouppop);
Rcpp::List calc_psivra(arma::vec current_dists,
			       arma::vec new_dists,
			       Rcpp::NumericVector pops,
			       Rcpp::NumericVector distswitch,
			       Rcpp::NumericVector grouppop,
			       double tgt_min,
			       double tgt_other);
Rcpp::List calc_psisimilar(arma::vec current_dists,
			   arma::vec new_dists,
			   arma::vec orig_dists,
			   Rcpp::NumericVector distswitch);
Rcpp::List calc_psicounty(arma::vec current_dists,
			  arma::vec new_dists,
			  arma::vec county_assignments,
			  arma::vec popvec);
Rcpp::List calc_psipartisan(arma::vec current_dists,
                            arma::vec new_dists,
                            Rcpp::IntegerVector rvote,
                            Rcpp::IntegerVector dvote,
                            std::string measure,
                            int ndists);
Rcpp::List calc_psiminority(arma::vec current_dists,
                      arma::vec new_dists,
                      Rcpp::NumericVector pops,
                      Rcpp::NumericVector grouppop,
                      int ndists,
                      Rcpp::NumericVector minorityprop);
Rcpp::List calc_psihinge(arma::vec current_dists,
                            arma::vec new_dists,
                            Rcpp::NumericVector pops,
                            Rcpp::NumericVector grouppop,
                            int ndists,
                            Rcpp::NumericVector minorityprop);

#endif

