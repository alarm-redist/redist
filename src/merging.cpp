/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: Functions for Merging Plans
********************************************************/

#include "merging.h"

double pick_adj_regions_to_merge(
    std::vector<std::pair<int, int>> const &adj_pairs_vec,
    std::pair<int, int> &adj_pair
){

    // Return 1 over the log of the size 
    return 3.4;
}


int run_merge_split_step( 
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    bool split_district_only,
    Tree &ust, int k_param,
    Plan &plan, int nsteps_to_run,
    double &lower, double upper, double target,
    std::vector<bool> &visited, std::vector<bool> &ignore
){
    // initialize number of successes to zero
    int num_successes = 0;

    // tracks 
    bool changed_plan = true;

    // initialize a region level graph
    Graph rg;

    // get all pairs of adjacent regions 
    std::set<std::pair<int, int>> adj_region_pairs_set;
    std::vector<std::pair<int, int>> adj_pairs_vec;
    // stores the pair of adjacent regions to merge
    std::pair<int, int> adj_pair;

    // Seed the random generator with the current time
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib;
        

    for (int j = 0; j < nsteps_to_run; j++){
        // Rprintf("Iter %d!\n", j);
        // If plan change from previous step then recreate the region graph and 
        // adjacent pairslist and update tree stuff 
        if(changed_plan){
            rg = get_region_graph(g, plan);
            adj_region_pairs_set.clear();
            for (int u = 0; u < rg.size(); u++){
                // iterate over neighbors
                for(auto v: rg.at(u)){
                    // store the edges such that smaller vertex is always first
                    if (u < v) {
                        adj_region_pairs_set.emplace(u, v);
                    } else {
                        adj_region_pairs_set.emplace(v, u);
                    }
                }
            }
            // Now convert to vec
            adj_pairs_vec.assign(adj_region_pairs_set.begin(), adj_region_pairs_set.end());
            // Define the range for random index
            distrib = std::uniform_int_distribution<size_t>(0, adj_pairs_vec.size() - 1);
        }

        // Pick two regions to merge
        // pick_regions_to_merge(adj_region_pairs)
        auto adj_pair = adj_pairs_vec.at(distrib(gen));
        // 
        int region1_id = adj_pair.first;
        int region2_id = adj_pair.second;

        // just make the first id the merged id
        int merged_id = region1_id;
        // create vector of vertex ids
        std::vector<int> region_ids = plan.region_ids;
        // Replace all occurance of region 2 with region 1
        for (int v = 0; v < plan.V; v++)
        {
            if(region_ids[v] == region2_id){
                region_ids[v] = region1_id;
            }
        }

        // now create the merged dval and population
        int merged_dval = plan.region_dvals.at(region1_id) + plan.region_dvals.at(region2_id);
        double merged_pop = plan.region_pops.at(region1_id) + plan.region_pops.at(region2_id);

        
        // Prepare to draw uniform spanning tree
        for (int i = 0; i < plan.V; i++){
            ignore[i] = region_ids.at(i) != merged_id;
        }
        int root;
        clear_tree(ust);

        // Get a uniform spanning tree drawn on that region
        int result = sample_sub_ust(g, ust, plan.V, root, visited, ignore, pop, lower, upper, counties, cg);
        // Return unsuccessful if tree not drawn
        if (result != 0){
            REprintf("SPANNING TREE NOT DRAWN SOMETHING WENT REALLY WRONG!!\n");
            changed_plan = false;
            continue;
        }

            // splitting related params
        int new_region1_tree_root, new_region2_tree_root;
        int new_region1_dval, new_region2_dval;
        double new_region1_pop, new_region2_pop;

        bool successful_edge_found;

        int max_potential_d;

        if(split_district_only){
            max_potential_d = 1;
        }else{
            max_potential_d = merged_dval - 1;
        }

        // try to get an edge to cut
        successful_edge_found = get_edge_to_cut(ust, root,
                        k_param, max_potential_d,
                        pop, region_ids, 
                        merged_id, merged_pop,
                        merged_dval,
                        lower, upper, target,
                        new_region1_tree_root, new_region1_dval, new_region1_pop,
                        new_region2_tree_root, new_region2_dval, new_region2_pop);

        if(successful_edge_found){
            Rprintf(
                "Yippeee! Succesfully merged region %d and %d to make "
                "a new region with %d dval and %.1f population. It was "
                "then split into (R%d, %d, %.1f) and (R%d, %d, %.1f) \n", 
                region1_id, region2_id,merged_dval,merged_pop,
                region1_id, new_region1_dval,  new_region1_pop, 
                region2_id, new_region2_dval,  new_region2_pop
                );
            num_successes++;
            // now update it using the old ids 
            update_plan_from_cut(
                    ust, plan, split_district_only,
                    new_region1_tree_root, new_region1_dval, new_region1_pop,
                    new_region2_tree_root, new_region2_dval, new_region2_pop,
                    region1_id, region2_id
                );
            Rprintf("out of cut update!\n");
            // changed plan successfully 
            changed_plan = true;
        }else{
            // plan didn't change 
            changed_plan = false;
        }

    }
    

    return num_successes;
}
