///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for constraint_calc_helper.cpp
///////////////////////////////////////

#ifndef CONSTRAINT_CALC_HELPER_H
#define CONSTRAINT_CALC_HELPER_H

Rcpp::NumericVector findBoundary(Rcpp::List fullList,
				 Rcpp::List conList);
Rcpp::List calc_psipop(arma::vec current_dists,
		       arma::vec new_dists,
		       Rcpp::NumericVector pops,
		       Rcpp::NumericVector distswitch);
Rcpp::List calc_psicompact(arma::vec current_dists,
			   arma::vec new_dists,
			   NumericVector distswitch,
			   std::string measure,
			   List aList,
			   NumericVector areas_vec,
			   List boundarylength_list,
			   NumericVector pops,
			   NumericMatrix ssdmat,
			   double denominator);
Rcpp::List calc_psisegregation(arma::vec current_dists,
			       arma::vec new_dists,
			       Rcpp::NumericVector pops,
			       Rcpp::NumericVector distswitch,
			       Rcpp::NumericVector grouppop);
Rcpp::List calc_psisimilar(arma::vec current_dists,
			   arma::vec new_dists,
			   arma::vec orig_dists,
			   Rcpp::NumericVector distswitch);
Rcpp::List calc_psicounty(arma::vec current_dists,
			  arma::vec new_dists,
			  arma::vec county_assignments);

#endif

