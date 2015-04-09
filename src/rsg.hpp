// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef RSG_HPP
#define RSG_HPP

#include <RcppArmadillo.h>

Rcpp::List rsg (List adj_list,
               arma::ivec adj_length,
               arma::vec population,
               int Ndistrict,
               double target_pop,
               double thresh
               );

#endif
