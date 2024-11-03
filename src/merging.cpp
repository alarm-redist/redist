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

int run_merge_split_step_on_a_plan(
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    bool const split_district_only,
    int const k_param,
    Plan &plan, int const nsteps_to_run,
    double const lower, double const upper, double const target)
{
    int V = plan.V;
    Tree ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V);

    // initialize number of successes to zero
    int num_successes = 0;

    // tracks 
    bool changed_plan = true;

    // get all pairs of adjacent regions 
    std::vector<std::pair<int, int>> adj_pairs_vec;

    // Seed the random generator with the current time
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib;

    // Create the Unif[0,1] value 
    std::random_device rd_u;  // Seed for randomness
    std::mt19937 gen_u(rd_u()); // Mersenne Twister generator
    std::uniform_real_distribution<> unif_rv(0.0, 1.0); // Uniform distribution on [0, 1]
    // get draw by doing unif_rv(rd_u);


    for (int j = 0; j < nsteps_to_run; j++){
        // Rprintf("Iter %d!\n", j);
        // If plan change from previous step then recreate the region graph and 
        // adjacent pairslist and update tree stuff 
        if(changed_plan){
            std::vector<bool> valid_regions;
            if(split_district_only){
                valid_regions.resize(plan.num_regions, false);
                valid_regions.at(plan.remainder_region) = true;
            }else{
                valid_regions.resize(plan.num_regions, true);
            }
            
            // Now convert to vec
            get_all_adj_pairs(
                g, adj_pairs_vec,
                plan.region_ids,
                valid_regions
            );

            // Define the range for random index
            distrib = std::uniform_int_distribution<size_t>(0, adj_pairs_vec.size() - 1);
        }

        // Pick two regions to merge
        // pick_regions_to_merge(adj_region_pairs)
        auto adj_pair = adj_pairs_vec.at(distrib(gen));
        // 
        int region1_id = adj_pair.first;
        int region2_id = adj_pair.second;

        int old_region1_dval = plan.region_dvals.at(region1_id);
        int old_region2_dval = plan.region_dvals.at(region2_id);

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
        int merged_dval = old_region1_dval + old_region2_dval;
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
            // REprintf("SPANNING TREE NOT DRAWN SOMETHING WENT REALLY WRONG!!\n");
            changed_plan = false;
            continue;
        }

            // splitting related params
        int new_region1_tree_root, new_region2_tree_root;
        int new_region1_dval, new_region2_dval;
        double new_region1_pop, new_region2_pop;

        bool successful_edge_found; 
        bool proposal_accepted = false;

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

        // If successful then do uniform thing 
        if(successful_edge_found){
            double u_draw = unif_rv(rd_u);
            double mh_ratio = .99;
            proposal_accepted = u_draw <= mh_ratio;
        }

        if(proposal_accepted){
            // Rprintf(
            //     "Yippeee! Succesfully merged region %d and %d to make "
            //     "a new region with %d dval and %.1f population. It was "
            //     "then split into (R%d, %d, %.1f) and (R%d, %d, %.1f) \n", 
            //     region1_id, region2_id,merged_dval,merged_pop,
            //     region1_id, new_region1_dval,  new_region1_pop, 
            //     region2_id, new_region2_dval,  new_region2_pop
            //     );
            num_successes++;
            // check if we need to change the number of districts or multidistricts
            
            // Check if region 1 changed from multidist to district
            if(old_region1_dval > 1 && new_region1_dval == 1){
                plan.num_multidistricts--;
                plan.num_districts++;
            }else if(old_region1_dval == 1 && new_region1_dval > 1){
                // check if region 1 changed from district to multidistrict
                plan.num_multidistricts++;
                plan.num_districts--;
            }

            // Check if region 2 changed from multidist to district
            if(old_region2_dval > 1 && new_region2_dval == 1){
                plan.num_multidistricts--;
                plan.num_districts++;
            }else if(old_region2_dval == 1 && new_region2_dval > 1){
                // check if region 2 changed from district to multidistrict
                plan.num_multidistricts++;
                plan.num_districts--;
            }

            // now update it using the old ids 
            update_plan_from_cut(
                    ust, plan, split_district_only,
                    new_region1_tree_root, new_region1_dval, new_region1_pop,
                    new_region2_tree_root, new_region2_dval, new_region2_pop,
                    region1_id, region2_id
                );
            
            

            // changed plan successfully 
            changed_plan = true;
        }else{
            // plan didn't change 
            changed_plan = false;
        }

    }
    

    return num_successes;
}

void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    std::vector<Plan> &plans_vec, 
    bool const split_district_only, int const k_param,
    int const nsteps_to_run,
    double const lower, double const upper, double const target,
    std::vector<int> &success_count_vec
){
    int M = (int) plans_vec.size();

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // Create variables needed for each 

        // store the number of succesful runs
        success_count_vec.at(i) = run_merge_split_step_on_a_plan( 
            g, counties, cg, pop,
            split_district_only,
            k_param,
            plans_vec.at(i), nsteps_to_run,
            lower, upper, target
        );

    });

    // Wait for all the threads to finish
    pool.wait();

    return;
}

