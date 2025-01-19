#include "forest_plan_type.h"


ForestPlan::ForestPlan(arma::subview_col<arma::uword> region_ids_col, 
               arma::subview_col<arma::uword> region_sizes_col, 
               int ndists, int num_regions, const arma::uvec &pop, 
               bool split_district_only,
              const Rcpp::List &initial_forest_adj_list):
              Plan(region_ids_col, region_sizes_col, ndists, num_regions, pop, split_district_only){

    if(num_regions > 1){
        forest_graph = list_to_graph(initial_forest_adj_list);
    }
    
              }