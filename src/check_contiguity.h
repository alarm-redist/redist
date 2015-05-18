#ifndef CHECK_CONTIG_H
#define CHECK_CONTIG_H

#include <RcppArmadillo.h>

int check_contiguity(Rcpp::List adj_list,
					Rcpp::IntegerVector p_neighbors,
					int p_neighbors_size,
					Rcpp::IntegerVector d_neighbors,
					int i_dist,
					Rcpp::IntegerVector member_dvec
                    );


#endif
