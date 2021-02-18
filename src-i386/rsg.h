#ifndef RSG_H
#define RSG_H

#include <RcppArmadillo.h>

Rcpp::List rsg (List adj_list,
               arma::ivec adj_length,
               arma::vec population,
               int Ndistrict,
               double target_pop,
               double thresh,
               int maxiter
               );

#endif
