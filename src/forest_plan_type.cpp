#include "forest_plan_type.h"


ForestPlan::ForestPlan(arma::subview_col<arma::uword> region_ids_col, 
               arma::subview_col<arma::uword> region_sizes_col, 
               int ndists, int num_regions, const arma::uvec &pop, 
               bool split_district_only,
              const Rcpp::List &initial_forest_adj_list):
              Plan(region_ids_col, region_sizes_col, ndists, num_regions, pop, split_district_only){

    if(num_regions == 1){
        forest_graph.resize(region_ids.n_elem);
    }
    if(num_regions > 1){
        forest_graph = list_to_graph(initial_forest_adj_list);
    }    
}


Graph ForestPlan::get_forest_adj(){
    return forest_graph;
}





void ForestPlan::update_vertex_info_from_cut(
        Tree const &ust, EdgeCut const cut_edge, 
        const int split_region1_id, const int split_region2_id
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



/*  
 * @title Get Valid Adjacent Regions to Effective Tree Boundary Lengths Map
 * 
 * Returns a unordered_map mapping valid pairs of adjacent regions to 
 * the length of the effect tree boundary between them in the forest
 * 
 * Finds all valid pairs of adjacent regions (meaning their sizes are valid 
 * to merge) and returns a hash map mapping the pairs of valid adjacent 
 * regions to the length of the effective tree boundary between them in `g`.
 * Recall this is the probability that we chose the actual split edge among
 * the tree made by merging each boundary edge
 * 
 * 
 * @param g A graph (adjacency list) passed by reference
 * @param plan A plan object
 * @param check_adj_to_regions A vector tracking whether or not we should 
 * check for edges adjacent to vertices in a region of a particular size. For
 * example, `check_adj_to_regions[i] == true` means we will attempt to find 
 * edges adjacent to any vertex in a region of size i. This vector is 1-indexed
 * meaning we don't need to subtract region size by 1 when accessing.
 * @param valid_merge_pairs A 2D `ndists+1` by `ndists+1` boolean matrix that
 * uses region size to check whether or not two regions are considered a valid
 * merge that can be counted in the map. For example `valid_merge_pairs[i][j]`
 * being true means that any regions where the sizes are (i,j) are considered
 * valid to merge. 2D matrix is 1 indexed (don't need to subtract region size)
 * 
 * @details No modifications to inputs made
 * @return A hash map mapping (std::pair<int,int> to double) that maps a pair of
 * region ids in the form (smaller id, bigger id) to the effective tree boundary length 
 * between them in `g`.
 */
std::unordered_map<std::pair<int, int>, double, bounded_hash> NEW_get_valid_pairs_and_tree_eff_boundary_len_map(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    const Plan &plan, const TreeSplitter &edge_splitter,
    std::unordered_map<std::pair<int, int>, double, bounded_hash> existing_pair_map = {}
){
    std::vector<bool> visited(map_params.V);
    std::vector<int> pops_below_vertex(map_params.V, 0);
    
    // Initialize unordered_map with num_region * 2.5 buckets
    // Hueristic. Bc we know planar graph has at most 3|V| - 6 edges
    int const init_bucket_size = std::ceil(2.5*plan.num_regions);

    // create the hash map
    std::unordered_map<std::pair<int, int>, double, bounded_hash> region_pair_map(
        init_bucket_size, bounded_hash(plan.num_regions-1)
        );

    int const V = map_params.V;

    // boolean for whether or not we should 
    bool const check_existing_map = existing_pair_map.size() != 0;


    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = plan.region_ids(v);
        auto v_region_size = plan.region_sizes(v_region_num);

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!splitting_schedule.check_adj_to_regions.at(v_region_size)){
            continue;
        }

        // get neighbors
        std::vector<int> nbors = map_params.g.at(v);

        // now iterate over its neighbors
        for (int v_nbor : nbors) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = plan.region_ids(v_nbor);

            // ignore if they are in the same region
            if(v_region_num == v_nbor_region_num) continue;
            // ignore if pair can't be merged
            auto v_nbor_region_size = plan.region_sizes(v_nbor_region_num);
            if(!splitting_schedule.valid_merge_pair_sizes[v_region_size][v_nbor_region_size]) continue;

            // else they are different so regions are adj
            // check if we need to be aware of double counting. IE if both regions
            // are ones where check adjacent is true then its possible to double count edges
            bool double_counting_not_possible = !splitting_schedule.check_adj_to_regions[v_nbor_region_size];
            // Rprintf("Neighbor size is %d and Double counting is %d\n", 
            //     (int) v_nbor_region_size, (int) double_counting_not_possible);

            // Now add this edge if either we can't double count it 
            // or `v_region_num` is the smaller of the pair 
            if(double_counting_not_possible || v_region_num < v_nbor_region_num){
                // Rprintf("Counting (%d,%d)\n", v, v_nbor);
                int smaller_region_id  = std::min(v_nbor_region_num,v_region_num);
                int bigger_region_id  = std::max(v_nbor_region_num,v_region_num);

                // if we need to see if the pair is already in the hash and its not then ignore
                if(check_existing_map && 
                    existing_pair_map.find({smaller_region_id, bigger_region_id}) == existing_pair_map.end()){
                        Rprintf("Skipping!!\n");
                        continue;
                    }
                

                auto merged_region_size = v_region_size+v_nbor_region_size;

                auto cut_size_bounds = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];
                int min_possible_cut_size = cut_size_bounds.first;
                int max_possible_cut_size = cut_size_bounds.second;

                // get the log probability
                double log_edge_selection_prob = get_log_retroactive_splitting_prob_for_joined_tree(
                    map_params, plan, edge_splitter,
                    visited, pops_below_vertex,
                    v, v_nbor,
                    min_possible_cut_size, max_possible_cut_size,
                    splitting_schedule.all_regions_smaller_cut_sizes_to_try[merged_region_size]);

                // Rprintf("Adding (%d,%d) w/ %.4f\n", v, nbor, std::exp(log_edge_selection_prob));
                // add the non log thing 
                region_pair_map[{smaller_region_id, bigger_region_id}] += std::exp(log_edge_selection_prob);
            }
        }
    }

    // adj_pairs_vec
    return region_pair_map;
}




std::vector<std::tuple<int, int, double>> ForestPlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter,
    std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map
) const{
    // get pairs with graph theoretic boundary 
    auto valid_adj_region_pairs_to_boundary_map = NEW_get_valid_pairs_and_tree_eff_boundary_len_map(
        map_params, splitting_schedule, 
        *this,tree_splitter, existing_pair_map
    );

    // declare the vector of tuples
    std::vector<std::tuple<int, int, double>> region_pairs_tuple_vec;
    region_pairs_tuple_vec.reserve(
        valid_adj_region_pairs_to_boundary_map.size()
    );

    for(auto const &a_pair: valid_adj_region_pairs_to_boundary_map){
        // add the region ids and the log of the boundary 
        region_pairs_tuple_vec.emplace_back(
            a_pair.first.first, 
            a_pair.first.second,
            std::log(a_pair.second)
        );
    }

    return region_pairs_tuple_vec;

}




double ForestPlan::get_log_eff_boundary_len(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, 
    const int region1_id, int const region2_id
) const{

    std::vector<bool> visited(map_params.V);
    std::vector<int> pops_below_vertex(map_params.V);

    int const V = map_params.V;
    int const merged_region_size = region_sizes(region1_id)+ region_sizes(region2_id);

    auto cut_size_bounds = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    double tree_selection_probs = 0.0;
    for (int v = 0; v < V; v++) {
        if (region_ids(v) != region1_id) continue; // Only count if starting vertex in region 1
        for (int nbor : map_params.g[v]) {
            if (region_ids(nbor) != region2_id) continue;

            double log_edge_selection_prob = get_log_retroactive_splitting_prob_for_joined_tree(
                map_params, *this, tree_splitter,
                visited, pops_below_vertex,
                v, nbor,
                min_possible_cut_size,max_possible_cut_size,
                splitting_schedule.all_regions_smaller_cut_sizes_to_try[merged_region_size]);
            
            tree_selection_probs += std::exp(log_edge_selection_prob);
        }
    }

    return std::log(tree_selection_probs);
}