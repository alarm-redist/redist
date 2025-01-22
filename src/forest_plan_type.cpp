#include "forest_plan_type.h"


ForestPlan::ForestPlan(arma::subview_col<arma::uword> region_ids_col, 
               arma::subview_col<arma::uword> region_sizes_col, 
               int ndists, int num_regions, const arma::uvec &pop, 
               bool split_district_only,
              const Rcpp::List &initial_forest_adj_list):
              Plan(region_ids_col, region_sizes_col, ndists, num_regions, pop, split_district_only){

    if(num_regions == 1){
        forest_graph.resize(V);
    }
    if(num_regions > 1){
        forest_graph = list_to_graph(initial_forest_adj_list);
    }    
}


Graph ForestPlan::get_forest_adj(){
    return forest_graph;
}


void ForestPlan::update_vertex_info_from_cut(
        Tree &ust, EdgeCut cut_edge, 
        const int split_region1_id, const int split_region2_id,
        bool split_district_only
){


    // Get the root of the tree associated with region 1 and 2
    int split_region1_tree_root, split_region2_tree_root;
    int split_region1_size, split_region2_size;
    int split_region1_pop, split_region2_pop;

    cut_edge.get_split_regions_info(
        split_region1_tree_root, split_region1_size, split_region1_pop,
        split_region2_tree_root, split_region2_size, split_region2_pop
    );

    // update the vertex labels and the tree
    assign_region_id_and_forest_from_tree(
        ust, region_ids, forest_graph,
        split_region1_tree_root, split_region1_id);

    // Rprintf("Root is %d and its size is %d\n", 
    //     split_region1_tree_root, (int) forest_graph.at(split_region1_tree_root).size());

    assign_region_id_and_forest_from_tree(
        ust, region_ids, forest_graph,
        split_region2_tree_root, split_region2_id);


    // Rprintf("Root is %d and its size is %d\n", 
    //     split_region2_tree_root, (int) forest_graph.at(split_region2_tree_root).size());

    return;
}
