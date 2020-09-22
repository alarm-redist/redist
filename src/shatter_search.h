#ifndef SHATTER_SEARCH_H
#define SHATTER_SEARCH_H
#include <Rcpp.h>
using namespace Rcpp;

int shatter_search(List adj_list,
                   int p,
                   int i_dist,
                   IntegerVector member_plist_i,
                   IntegerVector dist_assignment);


#endif


#ifndef CAN_SWAP_H
#define CAN_SWAP_H
#include <Rcpp.h>
using namespace Rcpp;

bool can_swap(List adj_list, 
                   int p, 
                   int j_dist, 
                   IntegerVector dist_assignment);


#endif
