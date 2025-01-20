#include "base_plan_type.h"


class ForestPlan : public Plan {


public:
    // constructor
    ForestPlan(arma::subview_col<arma::uword> region_ids_col, 
               arma::subview_col<arma::uword> region_sizes_col, 
               int ndists, int num_regions, const arma::uvec &pop, 
               bool split_district_only,
              const Rcpp::List &initial_forest_adj_list = {});
    

    // Clone method to create a copy of the ForestPlan object
    std::unique_ptr<Plan> deep_clone() const override {
        return std::make_unique<ForestPlan>(*this); // Copy the entire object
    }

    // We now need to keep track of trees as undirected graphs

    Graph get_forest_adj();

    void update_vertex_info_from_cut(
            Tree &ust, EdgeCut cut_edge, 
            const int split_region1_id, const int split_region2_id,
            bool split_district_only
    );

};