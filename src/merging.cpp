/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: Functions for Merging Plans
********************************************************/

#include "merging.h"



// Break into 


// Merges two regions together in a plan
// DOES NOT CHECK FOR ADJACENCY
// UNLESS SPLIT AFTER THIS BREAKS THE PLAN!!
void merge_regions(
    Plan &plan, 
    const int region1_id, const int region2_id, const int merged_id
){
    // check that the merged id is not greater than the number of regions -1
    // bc your new merged plan has num_regions-1 and since 0 indexing any new
    // merged id should be less than that
    if(merged_id > plan.num_regions-1 ){
        Rprintf("Merging id %d and size is %d \n", merged_id, plan.num_regions);
        throw Rcpp::exception("Merged id is too big");
    }

    // check the merged id is one of the ids of the region being merged
    if(merged_id != region1_id && merged_id != region2_id){
        throw Rcpp::exception("Merged id is not one of the two regions being merged");
    }
    int const V = plan.region_ids.n_elem;
    // update the region vertex ids
    for (size_t v = 0; v < V; v++)
    {
        if(plan.region_ids(v) == region2_id){
            plan.region_ids(v) = merged_id;
        } else if(plan.region_ids(v) == region1_id){
            plan.region_ids(v) = merged_id;
        }
    }

    // if either region is the remainder then make it the merged one
    // TODO make it so merged region is proper id. I guess that should be new num regions -1 
    // Since no remainder region anymore ok

    // update the number of regions 
    plan.num_regions--;


    // now update the region population and dval 
    int merged_dval = plan.region_sizes(region1_id)+ plan.region_sizes(region2_id);
    int merged_pop = plan.region_pops.at(region1_id) + plan.region_pops.at(region2_id);

    plan.region_sizes(merged_id) = merged_dval;
    plan.region_pops.at(merged_id) = merged_pop;

}


int run_merge_split_step_on_a_plan(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    bool split_district_only, std::string const merge_prob_type,
    Plan &plan, Plan &new_plan, 
    TreeSplitter &tree_splitter,
    int const nsteps_to_run)
{
    int V = plan.region_ids.n_elem;
    Tree ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V);

    // If its the final step then remainder doesn't matter so make split district only false 
    if(split_district_only && plan.num_regions == splitting_schedule.ndists){
        split_district_only = false;
    }

    // initialize number of successes to zero
    int num_successes = 0;


    // get all adjacent pairs in initial plan
    auto adj_pairs_and_boundary_len_vec = plan.get_valid_adj_regions_and_eff_log_boundary_lens(
        map_params, splitting_schedule, tree_splitter
    );

    

    // get sampler of adj pairs
    auto pair_index_sampler = get_adj_pair_sampler(
        plan,
        adj_pairs_and_boundary_len_vec,
        merge_prob_type
    );


    
    std::vector<std::tuple<int, int, double >> proposed_adj_pairs_and_boundary_len_vec;
    // get sampler of adj pairs in proposed plan
    std::discrete_distribution<> proposed_pair_index_sampler;


    // Seed the random generator with the current time
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create the Unif[0,1] value 
    std::random_device rd_u;  // Seed for randomness
    std::mt19937 gen_u(rd_u()); // Mersenne Twister generator
    std::uniform_real_distribution<> unif_rv(0.0, 1.0); // Uniform distribution on [0, 1]
    // get draw by doing unif_rv(rd_u);

    // Make sure the new plan has the same data as the old
    new_plan = plan;

    RNGState rng_state;

    for (int j = 0; j < nsteps_to_run; j++){


        // Now just pick one with the oldest index

        // Find the index of the smallest element so oldest
        // std::vector<int> nums(
        // plan.region_added_order.begin(), 
        // plan.region_added_order.begin() + plan.num_regions
        // );
        // auto min_iter = std::min_element(nums.begin(), nums.end());
        // size_t min_index = std::distance(nums.begin(), min_iter);

        // // Now go through the adjacency pairs and find the one 
        // // with the oldest region 
        // for (const auto& an_adj_pair : adj_pairs_vec) {
        //     if(an_adj_pair.first == min_index || an_adj_pair.second == min_index){
        //         adj_pair = an_adj_pair;
        //         break;
        //     }
        // }

        // Pick two regions to merge
        int pair_index = pair_index_sampler(gen);
        double index_selection_prob = pair_index_sampler.probabilities().at(pair_index);
        auto pair_triple = adj_pairs_and_boundary_len_vec.at(pair_index);

        // get info std::get<0>(adj_pairs_and_boundary_lens[i]);
        int region1_id = std::get<0>(pair_triple); int region2_id = std::get<1>(pair_triple);
        double old_region_boundary_length = std::get<2>(pair_triple);
        int old_region1_dval = plan.region_sizes(region1_id);
        int old_region2_dval = plan.region_sizes(region2_id);

        // make the merged region the one associated with the larger dval
        int merged_id; int new_split_region;

        if (old_region1_dval > old_region2_dval)
        {
            merged_id = region1_id; new_split_region = region2_id;
        }else{
            merged_id = region2_id; new_split_region = region1_id;
        }
        

        
        // REprintf("Pre-merged %d and %d into %d is \n", region1_id, region2_id, merged_id);
        // new_plan.Rprint();
        // REprintf("\n\n");
        
        // merge the regions in the new plan 
        merge_regions(new_plan, region1_id, region2_id, merged_id);

        // REprintf("Post-merged %d and %d into %d is \n", region1_id, region2_id, merged_id);
        // new_plan.Rprint();
        // REprintf("\n\n");

        // now attempt to split them
        // bool successful_split = attempt_region_split(map_params, ust,
        //     new_plan, merged_id,
        //     new_split_region,
        //     visited, ignore, 
        //     k_param, split_district_only);

        int min_region_cut_size, plan_specific_max_region_cut_size;
        throw Rcpp::exception("Fix cut size thing!");

        bool successful_split = new_plan.attempt_split(
                map_params, splitting_schedule,
                ust, tree_splitter, 
                visited, ignore, rng_state,
                min_region_cut_size, plan_specific_max_region_cut_size,
                splitting_schedule.all_regions_smaller_cut_sizes_to_try[new_plan.region_sizes(merged_id)],
                split_district_only,
                merged_id, new_split_region
            );


        bool proposal_accepted = false;

        // If successful then do MH adjustment
        if(successful_split){

            // get all valid adjacent pairs in the new proposed plan
            proposed_adj_pairs_and_boundary_len_vec = plan.get_valid_adj_regions_and_eff_log_boundary_lens(
                map_params, splitting_schedule, tree_splitter
            );

            // get sampler of adj pairs in proposed plan
            proposed_pair_index_sampler = get_adj_pair_sampler(
                new_plan,
                proposed_adj_pairs_and_boundary_len_vec,
                merge_prob_type
            );

            // index of (region1, region2) in proposed plan
            int region_pair_in_proposal_index;

            int smaller_region_id = std::min(region1_id, region2_id);
            int bigger_region_id = std::max(region1_id, region2_id);

            // Rprintf("Looking for Pair (%d, %d)!\n", smaller_region_id, bigger_region_id);

            // Find the index of the pair in the proposed plan
            for (size_t i = 0; i < proposed_adj_pairs_and_boundary_len_vec.size(); i++)
            {
                // Rprintf(" (%d, %d) ",
                //     proposed_adj_pairs_and_boundary_len_vec.at(i).at(0),
                //     proposed_adj_pairs_and_boundary_len_vec.at(i).at(1)
                //     );

                if(std::get<0>(proposed_adj_pairs_and_boundary_len_vec[i]) == smaller_region_id 
                && std::get<1>(proposed_adj_pairs_and_boundary_len_vec[i]) == bigger_region_id){
                    region_pair_in_proposal_index = i;
                    break;
                }

                // means disconnected thing glitch. If encountered just print and ignore
                if(i == (proposed_adj_pairs_and_boundary_len_vec.size()-1)){
                    REprintf("Merged and Split Regions %d and %d no longer adjacent!\n",
                        smaller_region_id, bigger_region_id);
                    smaller_region_id = -1;
                }
            }
            // check if merged and resplit regions no longer adj 
            if(smaller_region_id == -1){
                proposal_accepted = false;
            }else{
                // get length of boundary in proposed plan
                int new_region_boundary_length = std::get<2>(proposed_adj_pairs_and_boundary_len_vec[region_pair_in_proposal_index]);
                // get probability we would have picked that pair in the proposed plan
                double proposed_pair_index_selection_prob = proposed_pair_index_sampler.probabilities().at(region_pair_in_proposal_index);
                
                double log_mh_ratio = get_log_mh_ratio(
                    map_params.g, region1_id, region2_id,
                    old_region_boundary_length, new_region_boundary_length, 
                    index_selection_prob, proposed_pair_index_selection_prob, 
                    plan, new_plan
                );

                // draw a uniform and compare to log_mh_ratio
                double u_draw = unif_rv(rd_u);
                proposal_accepted = u_draw <= std::exp(log_mh_ratio);
            }
        }


        // if proposal accepted then make current plan the new one
        if(proposal_accepted){
            num_successes++;
            // update plan to new plan
            plan = new_plan;
            // if plan changed then the proposal adj pairs and sampler becomes 
            // the current one 
            adj_pairs_and_boundary_len_vec = proposed_adj_pairs_and_boundary_len_vec;
            pair_index_sampler = proposed_pair_index_sampler;
        }else{
            // plan didn't change so revert new plan back to old plan
            new_plan = plan;
        }

    }
    
    return num_successes;
}

void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &new_plan_ptrs_vec, 
    std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
    bool const split_district_only, std::string const merge_prob_type, 
    int const nsteps_to_run,
    Rcpp::IntegerMatrix::Column success_count_vec
){
    int M = (int) plan_ptrs_vec.size();

    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // Create variables needed for each 

        // store the number of succesful runs
        success_count_vec[i] = run_merge_split_step_on_a_plan(
            map_params, splitting_schedule,
            split_district_only, merge_prob_type,
            *plan_ptrs_vec.at(i), *new_plan_ptrs_vec.at(i), 
            *tree_splitters_ptr_vec.at(i),
            nsteps_to_run
        );

    });

    // Wait for all the threads to finish
    pool.wait();

    return;
}

