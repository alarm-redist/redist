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

    // update the region vertex ids
    for (size_t v = 0; v < plan.V; v++)
    {
        if(plan.region_ids(v) == region2_id){
            plan.region_ids(v) = merged_id;
        } else if(plan.region_ids(v) == region1_id){
            plan.region_ids(v) = merged_id;
        }
    }

    // if either region is the remainder then make it the merged one
    if(region1_id == plan.remainder_region || region2_id == plan.remainder_region){
        plan.remainder_region = merged_id;
    }

    // update the number of regions 
    plan.num_regions--;
    
    // update the num districts and multi-districts 

    // if two districts then -2 districts, +1 multi
    if(plan.region_dvals(region1_id) == 1 && plan.region_dvals(region2_id) == 1){
        plan.num_districts -= 2;
        plan.num_multidistricts++;
    }else if(plan.region_dvals(region1_id) == 1 || plan.region_dvals(region2_id) == 1){
        // else if only one is a district just decrease number of districts by 1
        plan.num_districts--;
    }else{
        plan.num_multidistricts--;
    }

    // now update the region population and dval 
    int merged_dval = plan.region_dvals(region1_id)+ plan.region_dvals(region2_id);
    double merged_pop = plan.region_pops.at(region1_id) + plan.region_pops.at(region2_id);

    plan.region_dvals(merged_id) = merged_dval;
    plan.region_pops.at(merged_id) = merged_pop;

}


int run_merge_split_step_on_a_plan(
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    bool split_district_only,
    int const k_param,
    Plan &plan, Plan &new_plan, int const nsteps_to_run,
    double const lower, double const upper, double const target)
{
    int V = plan.V;
    Tree ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V);

    // If its the final step then remainder doesn't matter so make split district only false 
    if(split_district_only && plan.num_regions == plan.N){
        split_district_only = false;
    }

    // initialize number of successes to zero
    int num_successes = 0;

    // tracks 
    bool changed_plan = true;

    // get all pairs of adjacent regions 
    std::vector<std::pair<int, int>> adj_pairs_vec;


    std::vector<std::pair<int, int>> proposal_adj_pairs_vec;


    // For the first run fill in the adjacency thing 
    std::vector<bool> valid_regions;
    if(split_district_only){
        // if splitting district only then only find adjacent to remainder
        valid_regions.resize(plan.num_regions, false);
        valid_regions.at(plan.remainder_region) = true;
    }else{
        valid_regions.resize(plan.num_regions, true);
    }

    // Now fill in proposal_adj_pairs_vec
    get_all_adj_pairs(
        g, proposal_adj_pairs_vec,
        plan.region_ids,
        valid_regions
    );

    // Seed the random generator with the current time
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib;

    // Create the Unif[0,1] value 
    std::random_device rd_u;  // Seed for randomness
    std::mt19937 gen_u(rd_u()); // Mersenne Twister generator
    std::uniform_real_distribution<> unif_rv(0.0, 1.0); // Uniform distribution on [0, 1]
    // get draw by doing unif_rv(rd_u);

    // Make sure the new plan has the same data as the old
    new_plan.shallow_copy(plan);

    for (int j = 0; j < nsteps_to_run; j++){
        // Rprintf("Iter %d!\n", j);
        // If plan change from previous step then recreate the region graph and 
        // adjacent pairslist and update tree stuff 
        if(changed_plan){
            // if plan changed then the proposal becomes the current one
            adj_pairs_vec = proposal_adj_pairs_vec;

            // TODO: add function which takes the adjacent pairs vec and returns 
            // a random sampler 

            // Define the range for random index
            distrib = std::uniform_int_distribution<size_t>(0, adj_pairs_vec.size() - 1);
        }

        // Pick two regions to merge
        // pick_regions_to_merge(adj_region_pairs)
        auto adj_pair = adj_pairs_vec.at(distrib(gen));

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
        

        // 
        int region1_id = adj_pair.first;
        int region2_id = adj_pair.second;

        int old_region1_dval = plan.region_dvals(region1_id);
        int old_region2_dval = plan.region_dvals(region2_id);

        // make the merged region the one associated with the larger dval
        int merged_id; int new_split_region;

        if (old_region1_dval > old_region2_dval)
        {
            merged_id = region1_id; new_split_region = region2_id;
        }else{
            merged_id = region2_id; new_split_region = region1_id;
        }
        

        if(split_district_only)
        {
            if(merged_id != plan.remainder_region){
                REprintf("ERROROROROR merged id is %d but remainder is %d",
                merged_id, plan.remainder_region);
            }
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
        bool successful_split = attempt_region_split(
            g, ust, counties, cg,
            new_plan, merged_id,
            new_split_region,
            visited, ignore, pop,
            lower, upper, target,
            k_param, split_district_only);

        bool proposal_accepted = false;

        // If successful then do uniform thing 
        if(successful_split){

            std::vector<bool> valid_regions;
            if(split_district_only){
                // if splitting district only then only find adjacent to remainder
                // which is region 2 because its the bigger one
                valid_regions.resize(plan.num_regions, false);
                valid_regions.at(new_plan.remainder_region) = true;
            }else{
                valid_regions.resize(plan.num_regions, true);
            }

            // get adjacent region pairs under the proposal
            get_all_adj_pairs(
                g, proposal_adj_pairs_vec,
                new_plan.region_ids,
                valid_regions
            );

            double u_draw = unif_rv(rd_u);
            double log_mh_ratio = get_log_mh_ratio(
                g, 
                region1_id, region2_id,
                plan.region_ids,
                new_plan.region_ids,
                static_cast<int>(adj_pairs_vec.size()), 
                static_cast<int>(proposal_adj_pairs_vec.size())
            );
            // Rprintf("(%d and %d and %.3f)\n", 
            //     static_cast<int>(adj_pairs_vec.size()),
            //     static_cast<int>(proposal_adj_pairs_vec.size()),
            //     std::exp(log_mh_ratio)
            //     );
            proposal_accepted = std::log(u_draw) <= log_mh_ratio;
        }


        // if proposal accepted then make current plan the new one
        if(proposal_accepted){
            num_successes++;
            // changed plan successfully 
            changed_plan = true;
            // update plan to new plan
            plan.shallow_copy(new_plan);
        }else{
            // plan didn't change 
            changed_plan = false;
            // revert new plan back to old plan
            new_plan.shallow_copy(plan);
        }

    }
    
    return num_successes;
}

void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    std::vector<Plan> &plans_vec, std::vector<Plan> &new_plans_vec, 
    bool const split_district_only, int const k_param,
    int const nsteps_to_run,
    double const lower, double const upper, double const target,
    Rcpp::IntegerMatrix::Column success_count_vec
){
    int M = (int) plans_vec.size();

    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // Create variables needed for each 

        // store the number of succesful runs
        success_count_vec[i] = run_merge_split_step_on_a_plan( 
            g, counties, cg, pop,
            split_district_only,
            k_param,
            plans_vec.at(i), new_plans_vec.at(i), nsteps_to_run,
            lower, upper, target
        );

    });

    // Wait for all the threads to finish
    pool.wait();

    return;
}

