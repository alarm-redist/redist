///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for constraint_calc_helper.cpp
///////////////////////////////////////

#ifndef CONSTRAINT_CALC_HELPER_H
#define CONSTRAINT_CALC_HELPER_H

Rcpp::List calc_psipop(arma::vec current_dists,
		       arma::vec new_dists,
		       Rcpp::NumericVector pops);
Rcpp::List calc_psicompact(arma::vec current_dists,
			   arma::vec new_dists,
			   Rcpp::NumericVector pops,
			   Rcpp::NumericMatrix ssdmat,
			   double denominator);
Rcpp::List calc_psisegregation(arma::vec current_dists,
			       arma::vec new_dists,
			       Rcpp::NumericVector pops,
			       Rcpp::NumericVector grouppop);
Rcpp::List calc_psisimilar(arma::vec current_dists,
			   arma::vec new_dists,
			   arma::vec orig_dists);
Rcpp::List calc_psicounty(arma::vec current_dists,
			  arma::vec new_dists,
			  arma::vec county_assignments);

#endif

