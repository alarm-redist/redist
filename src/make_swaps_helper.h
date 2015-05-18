///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for make_swaps_helper.cpp,
//          called by sw_mh_helper.cpp and sw_mh_alg.cpp
///////////////////////////////////////

#ifndef MAKE_SWAPS_HELPER_H
#define MAKE_SWAPS_HELPER_H

Rcpp::List adjcheck_propcd(Rcpp::List aList,
			   Rcpp::NumericVector prop_partitions,
			   Rcpp::NumericVector accepted_partitions,
			   Rcpp::NumericVector cds);
int elim_check(Rcpp::NumericVector prop_partition,
	       Rcpp::NumericVector cds);
int countpartitions(Rcpp::List aList);
Rcpp::NumericVector update_distpop(Rcpp::NumericVector prop_partition,
				   Rcpp::NumericVector unitpop_vec,
				   int prop_cd,
				   int curr_cd,
				   Rcpp::NumericVector distpop_vec);
double update_mhprob(Rcpp::NumericVector prop_partition,
		     Rcpp::List aList,
		     arma::vec cds,
		     int prop_cd,
		     double eprob,
		     double mh_prob);
Rcpp::NumericMatrix calcPWDh(Rcpp::NumericMatrix x);
Rcpp::NumericVector distParity(Rcpp::NumericMatrix mat,
			       Rcpp::NumericVector popvec);

#endif

