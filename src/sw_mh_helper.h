///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for sw_mh_helper.cpp, called by sw_mh_alg.cpp
///////////////////////////////////////

#ifndef sw_mh_helper_h
#define sw_mh_helper_h

Rcpp::NumericVector init_pop(Rcpp::NumericVector popvec,
			     arma::vec cds);
Rcpp::List genAlConn(Rcpp::List aList,
		     Rcpp::NumericVector cds);
Rcpp::NumericVector findBoundary(Rcpp::List fullList,
				 Rcpp::List conList);
Rcpp::List add_ties(Rcpp::List aList);
Rcpp::List cut_edges(Rcpp::List aList_con,
		     double eprob);
Rcpp::List bsearch_boundary(Rcpp::List aList,
			    arma::vec boundary);
int draw_p(int lambda);
Rcpp::List make_swaps(Rcpp::List boundary_cc,
		      Rcpp::List aList,
		      Rcpp::NumericVector cds_old,
		      Rcpp::NumericVector cds_orig,
		      Rcpp::NumericVector pop_vec,
		      Rcpp::NumericVector cd_pop_vec,
		      Rcpp::NumericVector group_pop_vec,
		      Rcpp::NumericMatrix ssdmat,
		      double minparity,
		      double maxparity,
		      int p,
		      double eprob,
		      double beta_population,
		      double beta_compact,
		      double beta_segregation,
		      double beta_similar,
		      double ssd_denominator);
int mh_decision(double mh_prob);
Rcpp::List changeBeta(arma::vec betavec,
		      double beta,
		      double constraint,
		      Rcpp::NumericVector weights,
		      int adjswap = 1);

#endif

