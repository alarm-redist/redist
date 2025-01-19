#include "base_plan_type.h"


class ForestPlan : public Plan {


public:
    // constructor
    ForestPlan(arma::subview_col<arma::uword> region_ids_col, 
               arma::subview_col<arma::uword> region_sizes_col, 
               int ndists, int num_regions, const arma::uvec &pop, 
               bool split_district_only,
              const Rcpp::List &initial_forest_adj_list = {});

    // We now need to keep track of trees as undirected graphs
    Graph forest_graph;
    std::vector<int> tree_roots;


};