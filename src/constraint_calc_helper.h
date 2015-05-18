///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for constraint_calc_helper.cpp
///////////////////////////////////////

#ifndef CONSTRAINT_CALC_HELPER_H
#define CONSTRAINT_CALC_HELPER_H

Rcpp::List calc_betapop(arma::vec current_dists,
			arma::vec new_dists,
			Rcpp::NumericVector pops,
			double beta_population,
			Rcpp::NumericVector distswitch);
Rcpp::List calc_betacompact(arma::vec current_dists,
			    arma::vec new_dists,
			    Rcpp::NumericVector pops,
			    double beta_compact,
			    Rcpp::NumericVector distswitch,
			    Rcpp::NumericMatrix ssdmat,
			    double denominator);
Rcpp::List calc_betasegregation(arma::vec current_dists,
				arma::vec new_dists,
				Rcpp::NumericVector pops,
				double beta_segregation,
				Rcpp::NumericVector distswitch,
				Rcpp::NumericVector grouppop);
Rcpp::List calc_betasimilar(arma::vec current_dists,
			    arma::vec new_dists,
			    arma::vec orig_dists,
			    double beta_similar,
			    Rcpp::NumericVector distswitch);

#endif

