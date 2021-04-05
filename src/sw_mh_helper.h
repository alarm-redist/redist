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
		      Rcpp::IntegerVector county_membership,
		      double minparity,
		      double maxparity,
		      int p,
		      double eprob,
		      double beta,
		      double weight_population,
		      double weight_compact,
		      double weight_segregation,
		      double weight_vra,
		      double weight_similar,
		      double weight_countysplit,
		      double weight_partisan,
		      double weight_minority,
		      double weight_hinge,
		      double ssd_denominator,
		      double tgt_min,
		      double tgt_other,
		      Rcpp::IntegerVector rvote,
		      Rcpp::IntegerVector dvote,
		      Rcpp::NumericVector minorityprop,
		      std::string compactness_measure,
		      std::string partisan_measure,
		      const Graph &g);
int mh_decision(double mh_prob);
Rcpp::List changeBeta(arma::vec betavec,
		      double beta,
		      double constraint,
		      Rcpp::NumericVector weights,
		      int adjswap = 1);

#endif

