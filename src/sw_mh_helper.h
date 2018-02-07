///////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/01/05
// Date Modified: 2015/02/26
// Purpose: header file for sw_mh_helper.cpp, called by sw_mh_alg.cpp
///////////////////////////////////////

#ifndef SW_MH_HELPER_H
#define SW_MH_HELPER_H


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
		      Rcpp::NumericVector cds_orig,
		      Rcpp::NumericVector pop_vec,
		      Rcpp::NumericVector cd_pop_vec,
		      Rcpp::NumericVector group_pop_vec,
		      Rcpp::NumericVector areas_vec,
		      arma::mat borderlength_mat,
		      Rcpp::NumericMatrix ssdmat,
		      Rcpp::NumericVector county_membership,
		      double minparity,
		      double maxparity,
		      int p,
		      double eprob,
		      double beta,
		      double weight_population,
		      double weight_compact,
		      double weight_segregation,
		      double weight_similar,
		      double weight_countysplit,
		      double ssd_denominator,
		      std::string compactness_measure);
int mh_decision(double mh_prob);
Rcpp::List changeBeta(arma::vec betavec,
		      double beta,
		      double constraint,
		      Rcpp::NumericVector weights,
		      int adjswap = 1);

#endif

