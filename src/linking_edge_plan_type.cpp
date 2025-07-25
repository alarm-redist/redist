/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/03
 * Purpose: Implementation of Plan type for Linking Edge space
 ********************************************************/

 #include "linking_edge_plan_type.h"

constexpr bool DEBUG_L_EDGE_PLANS_VERBOSE = false; // Compile-time constant


// randomly initalize a set of linking edges 
std::vector<std::pair<int, int>> get_intial_linking_edges(
    PlanMultigraph &plan_multigraph, PlanVector const &region_ids,
    int const num_regions, Graph &region_graph
){
    // First clear the region graph 
    for (size_t i = 0; i < region_graph.size(); i++)
    {
        region_graph[i].clear();
    }

    // Now build the plan multigraph 
    plan_multigraph.build_plan_multigraph(region_ids, num_regions);

    // Now iterate over the pairs and build the region graph 
    // Where two regions are adjacent iff they are a valid merge 
    for (auto const a_pair: plan_multigraph.pair_map.hashed_pairs){
        // Get the value 
        auto val = plan_multigraph.pair_map.get_value(a_pair.first, a_pair.second).second;
        // continue if not a valid hierarchical merge
        if(!val.merge_is_hier_valid) continue;
        // else add edge to the graph 
        region_graph[a_pair.first].push_back(a_pair.second);
        region_graph[a_pair.second].push_back(a_pair.first);
    }
    
    // Now we're going to make a random walk through the graph to make
    // a spanning tree. 
    std::vector<bool> visited(num_regions, false);
    std::vector<std::pair<int, int>> tree_edges;
    tree_edges.reserve(num_regions - 1);

    // DFS helper
    std::function<void(int)> dfs = [&](int v) {
        visited[v] = true;
        for (int u : region_graph[v]) {
            if (!visited[u]) {
                tree_edges.emplace_back(v, u);  // record tree edge
                dfs(u);
            }
        }
    };

    dfs(0);  // start from node 0

    std::vector<std::pair<int, int>> initial_edges;
    initial_edges.reserve(num_regions-1);

    // Now for each pair walk through the graph and find the first edge 
    // in those regions 
    for (auto const a_pair: tree_edges){
        // REprintf("(%d, %d)\n", a_pair.first, a_pair.second);
        bool pair_found = false;
        for (size_t v = 0; v < plan_multigraph.map_params.V; v++)
        {
            auto v_region = region_ids[v];
            // REprintf("%d in region %d\n", v, v_region);
            if(v_region != a_pair.first) continue;

            for(auto const u: plan_multigraph.map_params.g[v]){
                auto u_region = region_ids[u];
                if(u_region == a_pair.second){
                    // means we found edge between (region1, region2) so push back
                    initial_edges.push_back({v, u});
                    pair_found = true;
                    break;
                }
            }
            if(pair_found) break;
        }
        if(!pair_found) REprintf("ERROR! No edge found!\n");
    }


    return initial_edges;

}

VertexGraph LinkingEdgePlan::get_forest_adj(){
    return forest_graph;
}


LinkingEdgePlan::LinkingEdgePlan(
    int const total_seats, 
    int const total_pop,
    PlanVector &this_plan_region_ids, 
    RegionSizes &this_plan_region_sizes,
    IntPlanAttribute &this_plan_region_pops,
    IntPlanAttribute &this_plan_order_added
):
    Plan(total_seats, total_pop, 
        this_plan_region_ids, this_plan_region_sizes, this_plan_region_pops, this_plan_order_added
    )
{
    forest_graph.resize(region_ids.size());
    linking_edges.reserve(this_plan_region_sizes.size()-1);
};

LinkingEdgePlan::LinkingEdgePlan(
        int const ndists, int const num_regions,
        const arma::uvec &pop,
        PlanVector &this_plan_region_ids, 
        RegionSizes &this_plan_region_sizes,
        IntPlanAttribute &this_plan_region_pops,
        IntPlanAttribute &this_plan_order_added,
        TreeSplitter const &tree_splitter,
        USTSampler &ust_sampler,
        PlanMultigraph &plan_multigraph,
        Graph &region_graph,
        RNGState &rng_state, 
        const Rcpp::List &initial_forest_adj_list,
        const std::vector<std::array<double, 3>> &input_initial_linking_edges 
):
Plan(num_regions, pop, 
    this_plan_region_ids, this_plan_region_sizes, this_plan_region_pops, this_plan_order_added
){
    // resize the forest graph 
    forest_graph.resize(region_ids.size());


    if(initial_forest_adj_list.size() > 1){
        throw Rcpp::exception("Input Forest list not supported right now\n");
    }else{
        // else just build a forest at random 
        for (size_t region_id = 0; region_id < num_regions; region_id++)
        {
            int root;
            auto result = draw_tree_on_region(
                ust_sampler.map_params, region_id, 
                ust_sampler.ust, ust_sampler.visited, ust_sampler.ignore,
                root, rng_state, 1000000
            );

            if(!result.first){
                throw Rcpp::exception("Could not draw a tree on a region after 1000000 attempts");
            } 
        }
    
    }

    auto initial_linking_edges = get_intial_linking_edges(
        plan_multigraph, region_ids, num_regions, region_graph
    );

    // reserve space for ndists - 1 linking edges 
    linking_edges.reserve(ndists-1);

    for (auto const a_pair: initial_linking_edges){
        // REprintf("(%d, %d)\n", a_pair.first, a_pair.second);
        double log_prob = get_regions_log_splitting_prob(
            tree_splitter,  ust_sampler,
            a_pair.first, a_pair.second
        );
        linking_edges.push_back({a_pair.first, a_pair.second, log_prob});
    }
    
}





// NOTE: This really only belongs to Forest and linking edge class, not graph 
double LinkingEdgePlan::get_regions_log_splitting_prob(
    TreeSplitter const &tree_splitter, USTSampler &ust_sampler,
    const int region1_root, const int region2_root
) const{
    int region1_id = region_ids[region1_root]; int region2_id = region_ids[region2_root];
    int region1_size = region_sizes[region1_id]; int region2_size = region_sizes[region2_id];
    int merged_region_size = region1_size + region2_size;
    auto cut_size_bounds = ust_sampler.splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    // get the log probability
    return tree_splitter.get_log_retroactive_splitting_prob_for_joined_tree(
        ust_sampler.map_params, forest_graph, ust_sampler.stack,
        ust_sampler.visited, ust_sampler.pops_below_vertex,
        region1_root, region2_root,
        region_pops[region1_id], region_pops[region2_id],
        region1_size, region2_size,
        min_possible_cut_size, max_possible_cut_size,
        ust_sampler.splitting_schedule.all_regions_smaller_cut_sizes_to_try[merged_region_size]
    );
}


void LinkingEdgePlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter,
    USTSampler &ust_sampler, EdgeCut const cut_edge, 
    const int split_region1_id, const int split_region2_id,
    bool const add_region
){

    
    // If we're not adding a region we need to erase the old linking edge associated 
    // with these two regions if it exists 
    if(!add_region){
        if(DEBUG_L_EDGE_PLANS_VERBOSE){
        Rprintf("Adding Edge (%d, %d) - Regions (%d, %d)\n",
            cut_edge.cut_vertex, cut_edge.cut_vertex_parent,
            split_region1_id, split_region2_id);
        Rprintf("Current Linking Edges:\n");
        for (auto const &a_edge: linking_edges)
        {
            Rprintf("\tEdge(%d, %d) - Region (%d, %d)- Prob %f \n", 
                std::get<0>(a_edge), std::get<1>(a_edge), 
            region_ids[std::get<0>(a_edge)], region_ids[std::get<1>(a_edge)],
            std::exp(std::get<2>(a_edge)));
        }
        Rprintf("\n");
        }
        int erase_index;
        bool need_to_erase = false;
        std::set<int> merge_region_ids = {split_region1_id, split_region2_id};
        // Go through and find that pair 
        for (int i = 0; i < linking_edges.size(); i++)
        {
            // get the regions associated with the linking edge
            std::set<int> candidate_pairs = {
                static_cast<int>(region_ids[std::get<0>(linking_edges[i])]), 
                static_cast<int>(region_ids[std::get<1>(linking_edges[i])])
            };
            // if they match break
            if(candidate_pairs == merge_region_ids){
                erase_index = i;
                need_to_erase = true;
                break;
            }
        }
        // erase that pair of linking edges if needed
        if(need_to_erase){
            linking_edges.erase(linking_edges.begin()+erase_index);
        }else{
            Rprint(true);
            REprintf("Error could not find old edge!!\n");
            throw Rcpp::exception("NOOOO");
        }
    }


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
        ust_sampler.ust, region_ids, forest_graph,
        split_region1_tree_root, split_region1_id,
        ust_sampler.vertex_queue);

    assign_region_id_and_forest_from_tree(
        ust_sampler.ust, region_ids, forest_graph,
        split_region2_tree_root, split_region2_id,
        ust_sampler.vertex_queue);

    // TODO need to find the edge and update that stuff 
    // Now update the linking edge stuff 

    // First iterate through the old linking edges and update ones touching a split region
    for(auto &a_linking_edge: linking_edges){
        int v = std::get<0>(a_linking_edge); int u = std::get<1>(a_linking_edge);
        int edge_region1 = region_ids[v];
        int edge_region2 = region_ids[u];
        // If either of the vertices in the edge is now in a split region we need 
        // to update the probability that edge was chosen 
        if(edge_region1 == split_region1_id || edge_region2 == split_region1_id ||
            edge_region1 == split_region2_id || edge_region2 == split_region2_id){
            if(DEBUG_L_EDGE_PLANS_VERBOSE){
                Rprintf("(%d,%d) in Region (%d,%d)", v,u, edge_region1, edge_region2);
                Rprintf(" Linking edge prob was %f\n", std::get<2>(a_linking_edge));
            } 
            // update the log probability 
            std::get<2>(a_linking_edge) = get_regions_log_splitting_prob(
                tree_splitter, ust_sampler,
                u, v
            );
            if(DEBUG_L_EDGE_PLANS_VERBOSE) Rprintf("Linking edge prob is now %f\n", std::get<2>(a_linking_edge));
        }
    }

    // Append the new linking edge from the region we just split 
    linking_edges.push_back({cut_edge.cut_vertex, cut_edge.cut_vertex_parent, cut_edge.log_prob});
    if(DEBUG_L_EDGE_PLANS_VERBOSE){
        Rprintf("\nAdding New Edge (%d,%d) in Region (%d,%d), log prob %f\n",
            cut_edge.cut_vertex, cut_edge.cut_vertex_parent,
            region_ids[cut_edge.cut_vertex], region_ids[cut_edge.cut_vertex_parent],
            cut_edge.log_prob);
    } 

    return;
}



double LinkingEdgePlan::get_log_eff_boundary_len(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter, 
    const int region1_id, int const region2_id
) const {
    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        int edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        int edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        if(edge_region1 == region1_id && edge_region2 == region2_id){
            return std::get<2>(edge_pair);
        }
    }
    throw Rcpp::exception("Linking Pair not found!\n");
}


// - for linking edge sampling its the edge selection probability PLUS 
// the ratio log(merged plan linking edges) - log(plan linking edges)
std::vector<std::tuple<RegionID, RegionID, double>> LinkingEdgePlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function, 
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter
) const{
    // build the multigraph 
    plan_multigraph.build_plan_multigraph(region_ids, num_regions);
    // remove invalid hard constraint merges 
    plan_multigraph.remove_invalid_hard_constraint_pairs(*this, scoring_function);

    // compute log linking edges for the plan
    double const plan_log_linking_edge_term = plan_multigraph.compute_log_multigraph_tau(
        num_regions, scoring_function
    );


    std::vector<std::tuple<RegionID, RegionID, double>> valid_adj_region_pairs_to_boundary_map;
    valid_adj_region_pairs_to_boundary_map.reserve(linking_edges.size());

    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        int edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        int edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);

        int edge_region1_size = region_sizes[edge_region1]; int edge_region2_size = region_sizes[edge_region2];

        double log_linking_edge_ratio = -plan_log_linking_edge_term;
        
        // Add if its ok to merge 
        if(splitting_schedule.valid_merge_pair_sizes[edge_region1_size][edge_region2_size]){
            // sanity check that we never get invalid pair 
            auto search_result = plan_multigraph.pair_map.get_value(edge_region1, edge_region2);
            if(!search_result.first){
                REprintf("BIG BIG ERROR: A pair of regions with linking edge is somehow hierarchically invalid now!\n");
                throw Rcpp::exception("A pair of regions with linking edge is somehow hierarchically invalid now!\n");
            }else if(!search_result.second.merge_is_hier_valid){
                REprintf("ERROR: A pair of regions with linking edge is somehow hierarchically invalid now!\n");
                throw Rcpp::exception("A pair of regions with linking edge is somehow hierarchically invalid now!\n");
            }


            if(num_regions > 2){ 
            std::vector<int> merge_index_reshuffle(num_regions, 0);


            // TEMP just rebuild the multigraph 
            std::vector<RegionID> flattened_all_plans(plan_multigraph.map_params.V);
            PlanVector plan_region_ids(flattened_all_plans, 0, plan_multigraph.map_params.V);
            // REprintf("Size is %u!\n", plan_region_ids.size());

            // set merge reindex 
            int const merged_reindex = num_regions-2;
            for (int current_reindex = 0, i = 0; i < num_regions; i++){
                if(i == edge_region1 || i == edge_region2){
                    merge_index_reshuffle[i] = merged_reindex;
                }else{
                    merge_index_reshuffle[i] = current_reindex;
                    ++current_reindex;
                }
                // REprintf("Mapping %d to %d!\n", i, merge_index_reshuffle[i]);
            }

            // REprintf("Merging (%u, %u)\n", region1_id, region2_id);
            for (size_t i = 0; i < region_ids.size(); i++)
            {
                plan_region_ids[i] = merge_index_reshuffle[region_ids[i]];
                // REprintf("Plan %u | Merged %u \n", plan.region_ids[i], plan_region_ids[i]);
            }

            PlanMultigraph temp_multi(plan_multigraph.map_params, true);
            temp_multi.build_plan_multigraph(plan_region_ids, num_regions -1);

            double merged_tau = plan_multigraph.compute_merged_log_multigraph_tau(
                num_regions,
                edge_region1, edge_region2, scoring_function
            );
            double temp_tau = temp_multi.compute_log_multigraph_tau(num_regions-1, scoring_function);

            if(std::fabs(merged_tau - temp_tau) > 1e16){
                REprintf("Difference %f | Merging (%u, %u) | Merged Plan - Merged Code - %f, NO MERGE Code - %f and Equality Check = %s\n",
                    merged_tau - temp_tau,
                    edge_region1, edge_region2,
                merged_tau, temp_tau, 
                (merged_tau == temp_tau) ? "TRUE" : "FALSE" );
                Rprint(true);
                throw Rcpp::exception("DIFFERENT LOG TAU VALUES!\n");
            }
            // the ratio is merged/current so since its already current term we just add merged
            log_linking_edge_ratio += merged_tau;

            }

            // Now we add the edge probability and linking edge ratio 
            valid_adj_region_pairs_to_boundary_map.push_back(
                {edge_region1, edge_region2, std::get<2>(edge_pair) + log_linking_edge_ratio}
            );
        }
    }
    return(valid_adj_region_pairs_to_boundary_map);
}




std::vector<std::pair<RegionID,RegionID>> LinkingEdgePlan::get_valid_smc_merge_regions(
    PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule,
    ScoringFunction const &scoring_function
) const{
    // build the multigraph 
    plan_multigraph.build_plan_multigraph(region_ids, num_regions);

    // remove invalid hard constraint merges
    plan_multigraph.remove_invalid_hard_constraint_pairs(*this, scoring_function);

    std::vector<std::pair<RegionID,RegionID>> valid_adj_region;
    valid_adj_region.reserve(linking_edges.size());

    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        int edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        int edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);

        int edge_region1_size = region_sizes[edge_region1]; int edge_region2_size = region_sizes[edge_region2];

        
        // Add if its ok to merge 
        if(splitting_schedule.valid_merge_pair_sizes[edge_region1_size][edge_region2_size]){
            // sanity check that we never get invalid pair 
            auto search_result = plan_multigraph.pair_map.get_value(edge_region1, edge_region2);
            if(!search_result.first){
                REprintf("A pair of regions with linking edge is somehow hierarchically invalid now!\n");
                throw Rcpp::exception("A pair of regions with linking edge is somehow hierarchically invalid now!\n");
            }

            // add 
            valid_adj_region.push_back(
                    {edge_region1, edge_region2}
            );
        }
    }
    return(valid_adj_region);
}


std::pair<bool, std::vector<std::pair<RegionID,RegionID>>> LinkingEdgePlan::attempt_to_get_valid_mergesplit_pairs(
    PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule,
    ScoringFunction const &scoring_function
) const{
    // attempt to build valid multigraph 
    bool const result = plan_multigraph.build_plan_multigraph(region_ids, num_regions);
    // return false if not successful
    if(!result) return std::make_pair(false, std::vector<std::pair<RegionID,RegionID>>{}); 

    // else remove all the invalid hierarchical merge pairs for multigraph computation later
    // plan_multigraph.remove_invalid_hierarchical_merge_pairs(*this);
    plan_multigraph.remove_invalid_hard_constraint_pairs(*this, scoring_function);


    std::vector<std::pair<RegionID,RegionID>> valid_adj_region;
    valid_adj_region.reserve(linking_edges.size());

    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        auto edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        auto edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);

        int edge_region1_size = region_sizes[edge_region1]; int edge_region2_size = region_sizes[edge_region2];
        // REprintf("(%u, %u) \n", edge_region1, edge_region2);
        
        // Add if its ok to merge 
        if(splitting_schedule.valid_merge_pair_sizes[edge_region1_size][edge_region2_size]){
            // sanity check that we never get invalid pair 
            auto search_result = plan_multigraph.pair_map.get_value(edge_region1, edge_region2);
            if(!search_result.first){
                 plan_multigraph.Rprint();
                REprintf("Getting Mergesplie Pairs: A pair of regions with linking edge is somehow hierarchically invalid now!\n");
                throw Rcpp::exception("Getting Mergesplie Pairs: A pair of regions with linking edge is somehow hierarchically invalid now!\n");
            }

            valid_adj_region.push_back(
                {edge_region1, edge_region2}
            );
        }
    }

    return std::make_pair(true, valid_adj_region);
};


void LinkingEdgePlan::Rprint(bool verbose) const{
    Plan::Rprint(verbose);
    if(verbose){
        Rprintf("Current Linking Edges:\n");
        for (auto const &a_edge: linking_edges)
        {
            Rprintf("\tEdge(%d, %d) - Region (%d, %d)- Prob %f \n", 
                std::get<0>(a_edge), std::get<1>(a_edge), 
            region_ids[std::get<0>(a_edge)], region_ids[std::get<1>(a_edge)],
            std::exp(std::get<2>(a_edge)));
        }
        Rprintf("\n");
    }
}