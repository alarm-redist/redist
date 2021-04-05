#ifndef CLOSEST_ADJ_POP_H
#define CLOSEST_ADJ_POP_H

#include <Rcpp.h>
using namespace Rcpp;

int closest_adj_pop(IntegerVector adj,
                int i_dist,
                NumericVector g_prop);

#endif