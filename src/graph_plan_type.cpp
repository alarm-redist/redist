#include "graph_plan_type.h"


constexpr bool DEBUG_GRAPH_PLANS_VERBOSE = false; // Compile-time constant

void GraphPlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter,
    USTSampler &ust_sampler, EdgeCut const cut_edge, 
    const int split_region1_id, const int split_region2_id,
    bool const add_region
){
    // Get the root of the tree associated with region 1 and 2
    int split_region1_tree_root, split_region2_tree_root;
    int split_region1_size, split_region2_size;
    int split_region1_pop, split_region2_pop;

    cut_edge.get_split_regions_info(
        split_region1_tree_root, split_region1_size, split_region1_pop,
        split_region2_tree_root, split_region2_size, split_region2_pop
    );

    // update the vertex labels
    assign_region_id_from_tree(ust_sampler.ust, region_ids,
        split_region1_tree_root, split_region1_id,
        ust_sampler.vertex_queue);

    assign_region_id_from_tree(ust_sampler.ust, region_ids,
        split_region2_tree_root, split_region2_id,
        ust_sampler.vertex_queue);

    return;
}





/*  
 * @title Get Valid Adjacent Regions to Boundary Lengths Map
 * 
 * Returns a unordered_map mapping valid pairs of adjacent regions to 
 * the length of the boundary between them in the graph
 * 
 * Finds all valid pairs of adjacent regions (meaning their sizes are valid 
 * to merge) and returns a hash map mapping the pairs of valid adjacent 
 * regions to the length of the boundary between them in `g`
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
 * @param existing_pair_map A hash map mapping pairs of regions to a double. 
 * This is an optional parameter and if its empty then its ignored. If its not
 * empty then pairs already in the hash won't be counted added to the output.
 * 
 * @details No modifications to inputs made
 * @return A hash map mapping (std::pair<int,int> to int) that maps a pair of
 * region ids in the form (smaller id, bigger id) to the length of the boundary
 * between them in `g`.
 */




std::vector<std::tuple<RegionID, RegionID, double>> GraphPlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function, bool const is_final_split,
    USTSampler &ust_sampler, TreeSplitter &tree_splitter
) const{

    // build the multigraph 
    plan_multigraph.build_plan_multigraph(region_ids, num_regions);
    // remove invalid hierarchical merges if needed
    plan_multigraph.remove_invalid_hierarchical_merge_pairs(*this);
    // remove invalid size pairs 
    plan_multigraph.remove_invalid_size_pairs(*this, splitting_schedule);
    // remove invalid hard constraint pairs 
    plan_multigraph.remove_invalid_hard_constraint_pairs(*this, scoring_function, is_final_split);
    // now make the output vector 
    std::vector<std::tuple<RegionID, RegionID, double>> region_pairs_tuple_vec;
    region_pairs_tuple_vec.reserve(plan_multigraph.pair_map.num_hashed_pairs);

    for(auto const key_val_pair: plan_multigraph.pair_map.get_all_values()){
        
        double log_boundary_len;
        if(key_val_pair.second.admin_adjacent){
            // if administratively adjacent only count within same county
            log_boundary_len = std::log(key_val_pair.second.within_county_edges);
        }else{
            // else only count across county 
            log_boundary_len = std::log(key_val_pair.second.across_county_edges);
        }

        region_pairs_tuple_vec.emplace_back(
            key_val_pair.first.first, 
            key_val_pair.first.second, 
            log_boundary_len
        );
    }

    return region_pairs_tuple_vec;
    
}


double GraphPlan::get_log_eff_boundary_len(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    USTSampler &ust_sampler, TreeSplitter &tree_splitter, 
    ScoringFunction const &scoring_function,
    const int region1_id, int const region2_id
) const{
    // Return the log of the graph theoretic boundary
    auto search_result = plan_multigraph.pair_map.get_value(region1_id, region2_id);

    if(search_result.second.admin_adjacent){
        // if administratively adjacent only count within same county
        return std::log(search_result.second.within_county_edges);
    }else{
        // else only count across county 
        return std::log(search_result.second.across_county_edges);
    }
}





/*
 * Calculate the deviations for cutting at every edge in a spanning tree.
 * and returns them ordered.
//'
//'
//' For each edge it returns the larger of the two deviations associated with 
//' the best region sizes assignment 
 */
std::vector<double> get_ordered_tree_cut_devs(Tree &ust, int root,
                             std::vector<int> const &cut_below_pop, double const target,
                             PlanVector const &region_ids,
                             RegionID const region_id1, RegionID const region_id2, 
                             int const region_size, int const region_pop,
                             int const min_potential_cut_size, int const max_potential_cut_size,
                             std::vector<int> const &smaller_cut_sizes_to_try
                             ) {
    int V = cut_below_pop.size();
    // compile a list of candidate edges to cut
    std::vector<double> devs; 
    devs.reserve(V); // reserve V which is overkill but ok
    // REprintf("Only looking for region id %d!\n", region_id);
    // REprintf("Starting at %d and %2.f and ", region_pop, static_cast<double>(region_pop)*2.0);
    for (int i = 0; i < V; i++) {
        // ignore vertices not in the region
        if (i == root || (region_ids[i] != region_id1 && region_ids[i] != region_id2)) continue;


        // start at total pop since deviance will never be more than total_pop/2
        double smallest_dev = static_cast<double>(region_pop)*2.0;
        
        int below_pop = cut_below_pop.at(i);
        int above_pop = region_pop - below_pop;

        // if one of the populations is zero just skip 
        if(below_pop == 0 || above_pop == 0){
            continue;
        } 


        // for each possible d value get the deviation above and below
        for(auto const cut_region1_size: smaller_cut_sizes_to_try){
            int cut_region2_size = region_size - cut_region1_size;
            double cut_region1_target = target*cut_region1_size; 
            double cut_region2_target = target*cut_region2_size;


            // REprintf("Size 1=%d-Target %.2f, Size 2=%d-Target %.2f", cut_region1_size, cut_region1_target, cut_region2_size, cut_region2_target);
            // find the larger deviation associated with assigning cut_region1_size to cutting below 
            double max_below_dev = std::max(
                std::fabs(below_pop - cut_region1_target) / cut_region1_target,
                std::fabs(above_pop - cut_region2_target) / cut_region2_target
            );
            // find the larger deviation associated with assigning cut_region1_size to cutting above
            double max_above_dev = std::max(
                std::fabs(above_pop - cut_region1_target) / cut_region1_target,
                std::fabs(below_pop - cut_region2_target) / cut_region2_target
            );
            // REprintf("pair dev %.2f, %.2f", max_below_dev, max_above_dev);

            // take this minimum of this d value and all previous ones 
            double min_dev = std::min(
                max_below_dev,
                max_above_dev
            ) ;
            smallest_dev = std::min(smallest_dev, min_dev);
        }

        devs.push_back(smallest_dev);
    }

    std::sort(devs.begin(), devs.end());
    // legacy code expects devs to be length V-1 so just pad out the remainder with junk
    devs.resize(V-1, 2.06);

    return devs;
}



/*
 * Choose k and multiplier for efficient, accurate sampling
 * Assumes plan multigraph is already built 
 */
int estimate_mergesplit_cut_k(
    Plan const &plan, PlanMultigraph const &plan_multigraph,  
    SplittingSchedule const &splitting_schedule,
    double const thresh, double const tol,
    RNGState &rng_state) {

    int k;
    // sample some spanning trees and compute deviances
    int V = plan_multigraph.map_params.g.size();
    // IN FUTURE USE MY OWN FUNCTION THIS HAS INEXING ERRORS
    // Graph dist_g = district_graph(g, region_ids, n_distr, true);
    int k_max = std::min(20 + ((int) std::sqrt(V)), V - 1); // heuristic
    int N_adapt = (int) std::floor(4000.0 / sqrt((double) V));
    

    double lower = plan_multigraph.map_params.target * (1 - tol);
    double upper = plan_multigraph.map_params.target * (1 + tol);
    
    std::vector<std::vector<double>> devs;
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    std::vector<bool> visited(V);
    std::vector<int> cut_below_pop(V,0);
    TreePopStack pop_stack(V+1);
    Tree county_tree = init_tree(plan_multigraph.map_params.num_counties);
    TreePopStack county_stack(plan_multigraph.map_params.num_counties);
    arma::uvec county_pop(plan_multigraph.map_params.num_counties, arma::fill::zeros);
    std::vector<std::vector<int>> county_members(plan_multigraph.map_params.num_counties, std::vector<int>{});
    std::vector<bool> c_visited(plan_multigraph.map_params.num_counties, true);
    std::vector<int> cty_pop_below(plan_multigraph.map_params.num_counties, 0);
    std::vector<std::array<int, 3>> county_path;
    std::vector<int> path;

    int max_V = 0;
    Tree ust = init_tree(V);
    for (int i = 0; i < N_adapt; i++) {
        double joint_pop = 0;
        auto random_pair_index = rng_state.r_int(plan_multigraph.pair_map.hashed_pairs.size());
        auto a_pair = plan_multigraph.pair_map.hashed_pairs[random_pair_index];
        auto merged_size = plan.region_sizes[a_pair.first] + plan.region_sizes[a_pair.second];
        auto merged_pop = plan.region_pops[a_pair.first] + plan.region_pops[a_pair.second];

        int n_vtx = 0;
        for (int j = 0; j < V; j++) {
            if (plan.region_ids[j] == a_pair.first || plan.region_ids[j] == a_pair.second) {
                joint_pop += plan_multigraph.map_params.pop(j);
                ignore[j] = false;
                n_vtx++;
            } else {
                ignore[j] = true;
            }
        }
        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(
            plan_multigraph.map_params, ust, root, 
            lower * merged_size, upper * merged_size,
            visited, ignore, county_tree, county_stack, county_pop, county_members, 
            c_visited, cty_pop_below, county_path, path,
            rng_state);
        if (result != 0) {
            i--;
            continue;
        }

        // reset the cut below pop to zero
        std::fill(cut_below_pop.begin(), cut_below_pop.end(), 0);
        get_tree_pops_below(ust, root, pop_stack, plan_multigraph.map_params.pop, cut_below_pop);

        std::pair<int, int> min_and_max_possible_cut_sizes = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_size];
        int min_possible_cut_size = min_and_max_possible_cut_sizes.first;
        int max_possible_cut_size = min_and_max_possible_cut_sizes.second;

        devs.push_back(
           get_ordered_tree_cut_devs(
                ust, root, cut_below_pop, plan_multigraph.map_params.target, 
                plan.region_ids,
                a_pair.first, a_pair.second, merged_size, 
                merged_pop,
                min_possible_cut_size, max_possible_cut_size,
                splitting_schedule.all_regions_smaller_cut_sizes_to_try[merged_size]
            )
        );
        
        int n_ok = 0;
        // NOTE very hacky, not sure if need to change because MMD now
        for (int j = 0; j < V-1; j++) {
            // REprintf("i = %d, j = %d\n", i, j);
            n_ok += devs.at(i).at(j) <= tol;
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max)
            max_ok = n_ok;
    }

    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    for (k = 1; k <= k_max; k++) {
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            int rand_index = rng_state.r_int(k);
            double dev = devs.at(i).at(rand_index);
            // REprintf("i = %d | Dev is %f", i, dev);
            if (dev > tol) continue;
            else n_ok++;
            for (int j = 0; j < N_adapt; j++) {
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k == k_max + 1) {
        Rcerr << "Warning: maximum hit; falling back to naive k estimator.\n";
        k = max_ok + 1;
    }

    k = std::min(k, max_V - 1);
    return(k);
}

/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void estimate_cut_k(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    RNGState &rng_state,
    int &k, int const last_k, 
    const arma::vec &unnormalized_weights, double thresh,
    double tol, std::vector<std::unique_ptr<Plan>> const &plan_ptrs_vec, 
    bool split_district_only,
    int const verbosity) {
    // sample some spanning trees and compute deviances
    int V = map_params.V;
    int k_max = std::min(10 + (int) (2.0 * V * tol), last_k + 4); // heuristic
    int N_max = plan_ptrs_vec.size();
    int N_adapt = std::min(60 + (int) std::floor(5000.0 / sqrt((double)V)), N_max);

    double target = map_params.target;
    double lower = target * (1 - tol);
    double upper = target * (1 + tol);
    int num_regions = plan_ptrs_vec[0]->num_regions;

    std::vector<std::vector<double>> devs;
    devs.reserve(N_adapt);
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    std::vector<bool> visited(V);
    std::vector<int> cut_below_pop(V,0);
    TreePopStack pop_stack(V+1);
    Tree county_tree = init_tree(map_params.num_counties);
    TreePopStack county_stack(map_params.num_counties);
    arma::uvec county_pop(map_params.num_counties, arma::fill::zeros);
    std::vector<std::vector<int>> county_members(map_params.num_counties, std::vector<int>{});
    std::vector<bool> c_visited(map_params.num_counties, true);
    std::vector<int> cty_pop_below(map_params.num_counties, 0);
    std::vector<std::array<int, 3>> county_path;
    std::vector<int> path;

    int idx = 0;
    int max_V = 0;
    Tree ust = init_tree(V);

    bool any_size_split = splitting_schedule.schedule_type == SplittingSizeScheduleType::AnyValidSizeSMD;
    
    for (int i = 0; i < N_max && idx < N_adapt; i++, idx++) {
        if (unnormalized_weights.at(i) == 0) { // skip if not valid
            idx--;
            continue;
        }

        int n_vtx = V;

        // Get the index of the region with the largest dval
        int biggest_region_id; int biggest_region_size; 

        // if split district only just do remainder 
        if(split_district_only){
            biggest_region_id = plan_ptrs_vec.at(i)->num_regions-1;
            biggest_region_size = plan_ptrs_vec.at(i)->region_sizes[biggest_region_id];
        }else if(any_size_split){
            // else get 
            auto max_it = std::max_element(
                plan_ptrs_vec.at(i)->region_sizes.begin(),
                plan_ptrs_vec.at(i)->region_sizes.begin() + num_regions
            );
            biggest_region_id = std::distance(
                plan_ptrs_vec.at(i)->region_sizes.begin(), 
                max_it
            );
            biggest_region_size = plan_ptrs_vec.at(i)->region_sizes[biggest_region_id];
        }else{// custom size 
            biggest_region_size = -1; biggest_region_id = num_regions -1;
            for (int j = 0; j < num_regions; j++)
            {
                int region_size = plan_ptrs_vec[i]->region_sizes[j];
                // check if valid and bigger
                if(splitting_schedule.valid_region_sizes_to_split[region_size] && 
                   region_size > biggest_region_size){
                    biggest_region_id = j;
                    biggest_region_size = region_size;
                }
            }
        }

        

        int biggest_size_region_pop = plan_ptrs_vec[i]->region_pops[biggest_region_id];
        std::pair<int, int> min_and_max_possible_cut_sizes = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[biggest_region_size];
        int min_possible_cut_size = min_and_max_possible_cut_sizes.first;
        int max_possible_cut_size = min_and_max_possible_cut_sizes.second;
        
        for (int j = 0; j < V; j++) {
            // if not the biggest region mark as ignore
            if (plan_ptrs_vec.at(i)->region_ids[j] != biggest_region_id) {
                ignore[j] = true;
                n_vtx--;
            }
        }
        

        // Rprintf("Tree on region %d of size %d has %d vertices and pop %d!\n",
        // biggest_region_id,
        //     biggest_region_size, n_vtx, plan_ptrs_vec.at(i)->region_pops.at(biggest_region_id));

        // plan_ptrs_vec.at(i)->Rprint();
        // Rprintf("\n\n");

        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(
            map_params, ust, root, 
            lower * min_possible_cut_size, upper * min_possible_cut_size,
            visited, ignore, county_tree, county_stack, county_pop, county_members, 
            c_visited, cty_pop_below, county_path, path,
            rng_state
        );


        


            
        if (result != 0) {
            idx--;
            continue;
        }

        // reset the cut below pop to zero
        std::fill(cut_below_pop.begin(), cut_below_pop.end(), 0);
        get_tree_pops_below(ust, root, pop_stack, map_params.pop, cut_below_pop);
        
        devs.push_back(
            get_ordered_tree_cut_devs(ust, root, cut_below_pop, target, 
                    plan_ptrs_vec.at(i)->region_ids,
                    biggest_region_id, biggest_region_id, biggest_region_size, 
                    biggest_size_region_pop,
                    min_possible_cut_size, max_possible_cut_size,
                    splitting_schedule.all_regions_smaller_cut_sizes_to_try[biggest_region_size]
                )
            );

        int n_ok = 0;
        for(const auto &a_dev : devs.at(idx)){
            if (a_dev <= tol) { // sorted
                n_ok++;
            } else {
                break;
            }
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max){
            max_ok = n_ok;
        }
            
        Rcpp::checkUserInterrupt();
    }

    if (idx < N_adapt) N_adapt = idx; // if rejected too many in last step
    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    for (k = 1; k <= k_max; k++) {
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            int rand_index = rng_state.r_int(k);
            double dev = devs.at(i).at(rand_index);
            if (dev > tol) continue;
            else n_ok++;
            for (int j = 0; j < N_adapt; j++) {
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k >= k_max) {
        if (verbosity >= 3) {
            Rcout << " [maximum hit; falling back to naive k estimator]";
        }
        k = max_ok;
    }

    if (last_k < k_max && k < last_k * 0.6) k = (int) (0.5*k + 0.5*last_k);

    k = std::min(std::max(max_ok + 1, k) + 1 - (distr_ok(k) > 0.99) + (thresh == 1),
                 max_V - 1);
}


