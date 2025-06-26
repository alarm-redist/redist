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

#endif

