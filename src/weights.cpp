/********************************************************
* Author: Philip O'Sullivan
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: SMC weight calculation related functions
********************************************************/

#include "weights.h"

//' Computes the effective sample size from log incremental weights
//'
//' Takes a vector of log incremental weights and computes the effective sample
//' size which is the sum of the weights squared divided by the sum of squared
//' weights
//'
//'
//' @title Compute Effective Sample Size
//'
//' @param log_wgt vector of log incremental weights
//'
//' @details No modifications to inputs made
//'
//' @return sum of weights squared over sum of squared weights (sum(wgt)^2 / sum(wgt^2))
//'
double compute_n_eff(const std::vector<double> &log_wgt) {
    double sum_wgt = 0.0;
    double sum_wgt_squared = 0.0;

    // compute sum of squares and square of sum
    for (const double& log_w : log_wgt) {
        double wgt = std::exp(log_w);
        sum_wgt += wgt;
        sum_wgt_squared += wgt * wgt;
    }


    return std::exp(
        (2 * std::log(sum_wgt)) - std::log(sum_wgt_squared)
    );
}

// only gets regions adjacent to indices where valid_regions is true
void get_all_adj_pairs(
    Graph const &g, std::vector<std::pair<int, int>> &adj_pairs_vec,
    std::vector<int> const &vertex_region_ids,
    std::vector<bool> const valid_regions
){
    int V = g.size();

    // Set to put adjacent pairs in 
    std::set<std::pair<int, int>> adj_region_pairs;

    // iterate over all vertices in g
    for (int i = 0; i < V; i++) {
        // Find out which region this vertex corresponds to
        int region_num_i = vertex_region_ids.at(i);

        // check if its a valid region and if not continue
        if(!valid_regions.at(region_num_i)){
            continue;
        }

        std::vector<int> nbors = g.at(i);

        // now iterate over its neighbors
        for (int nbor : nbors) {
            // find which region neighbor corresponds to
            int region_num_j = vertex_region_ids.at(nbor);

            // if they are different regions mark matrix true since region i
            // and region j are adjacent as they share an edge across
            if (region_num_i != region_num_j) {
                if (region_num_i < region_num_j) {
                    adj_region_pairs.emplace(region_num_i, region_num_j);
                } else {
                    adj_region_pairs.emplace(region_num_j, region_num_i);
                }
            }
        }
    }

    // Now convert to vec
    adj_pairs_vec.assign(adj_region_pairs.begin(), adj_region_pairs.end());

    return;

}


//' Returns the log of the count of the number of edges across two regions in
//' the underlying graph.
//'
//' Given a graph and two regions in the graph this function counts the number of
//' edges across the the two regions (ie one vertex in one and the other in
//' the other also know as the grap theoretic length of the boundary) and
//' returns the log of that count.
//'
//'
//' @title Log Region boundary count
//'
//' @param g A graph (adjacency list) passed by reference
//' @param vertex_region_ids A vector mapping vertex id to the region its in 
//' so `vertex_region_ids[i] = r` means vertex i is in region r
//' @param region1_id The id of the first region
//' @param region2_id The id of the second region
//'
//' @details No modifications to inputs made
//'
//' @return the log of the boundary count
//'
 double region_log_boundary(const Graph &g, 
                            const std::vector<int> &vertex_region_ids,
                            int const&region1_id,
                            int const&region2_id
) {
     int V = g.size(); // get number of vertices
     double count = 0; // count of number of edges across the two regions

     // Iterate over every vertex in the graph
     for (int i = 0; i < V; i++) {
         /* Check vertex if in first region, if not continue
          Since edges appear twice in adjacency list to avoid double
          counting we will only count those where first edge is in
          region 1
          */
         if (vertex_region_ids.at(i) != region1_id) continue;

         // Get vertice's neighbors
         std::vector<int> nbors = g[i];

         // Since edges show up twice to avoid double counting we will only count
         // edges where first one is region1_id and second is region2_id

         // Now check if neighbors are in second region
         for (int nbor : nbors) {
             if (vertex_region_ids.at(nbor) != region2_id)
                 continue;
             // if they are increase count by one
             count += 1.0;
         }
     }

     return std::log(count);
}

// computes log metropolis hastings ratio
double get_log_mh_ratio(
    const Graph &g, 
    const int region1_id, const int region2_id,
    const std::vector<int> &old_vertex_region_ids,
    const std::vector<int> &new_vertex_region_ids,
    const int num_old_adj_regions, const int num_new_adj_regions
){
    // get the log boundary lengths
    double old_log_boundary_length = region_log_boundary(
        g, old_vertex_region_ids,
        region1_id, region2_id
    );

    double new_log_boundary_length = region_log_boundary(
        g, new_vertex_region_ids,
        region1_id, region2_id
    );

    // Get the merge probability which for now is uniform
    double log_old_merge_prob = std::log(1) - std::log(
        static_cast<double>(num_old_adj_regions)
        );

    double log_new_merge_prob = std::log(1) - std::log(
        static_cast<double>(num_new_adj_regions)
        );

    double log_mh_ratio_numerator = old_log_boundary_length + log_old_merge_prob;
    double log_mh_ratio_denominator = new_log_boundary_length + log_new_merge_prob;

    double log_mh_ratio = log_mh_ratio_numerator - log_mh_ratio_denominator;

    return log_mh_ratio;
}


//' Get the probability the union of two regions was chosen to split
//'
//' Given a plan object and two regions in the plan this returns the probability
//' the union of the two regions was chosen to be split.
//'
//' @title Get Retroactive Split Selection Probability
//'
//' @param plan A plan object
//' @param region1_id The id of the first region to union
//' @param region2_id The id of the second region to union
//'
//' @details No modifications to inputs made
//'
//' @return the log of the probability the union of the two regions would be
//' chosen to split.
//'
double get_log_retroactive_splitting_prob(
        const Plan &plan,
        const int region1_id, const int region2_id
){

    // count total dvals of all the multidistricts
    int total_multi_ds = 0;


    // Iterate over all regions
    for(int region_id = 0; region_id < plan.num_regions; region_id++) {
        int d_val = plan.region_dvals.at(region_id);

        // collect info if multidistrict and not the two we started with
        if(d_val > 1 && region_id != region1_id && region_id != region2_id){
            // Add that regions d value to the total
            total_multi_ds += d_val;
        }
    }

    // Now get the sum of dnk values of two regions aka the unioned old region
    int unioned_region_dnk = plan.region_dvals.at(region1_id) + plan.region_dvals.at(region2_id);
    // update the total number of multi district dvals with the value the union region would have been
    total_multi_ds += unioned_region_dnk;

    // so prob of picking is the sum of the two regions dnk over the dnk of all
    // multidistricts

    double log_prob = std::log(
        static_cast<double>(unioned_region_dnk)
    ) - std::log(
            static_cast<double>(total_multi_ds)
    );

    return log_prob;

}

// TODO: DOCUMENTATION NEEDED
// Computes log population tempering term
double compute_log_pop_temper(
        const Plan &plan,
        const int region1_id, const int region2_id,
        double target, double pop_temper
){
    int N = plan.N;
    // Get pop and dval of two regions and their union
    double region1_pop = plan.region_pops[region1_id];
    double region2_pop = plan.region_pops[region2_id];
    double old_region_pop = region1_pop + region2_pop;

    int region1_dval = plan.region_dvals[region1_id];
    int region2_dval = plan.region_dvals[region2_id];
    int old_region_dval = region1_dval + region2_dval;

    // now compute the devations for each of them
    double dev1 = std::fabs(
        region1_pop/static_cast<double>(region1_dval) - target
    )/target;
    double dev2 = std::fabs(
        region2_pop/static_cast<double>(region2_dval) - target
    )/target;
    double old_dev = std::fabs(
        old_region_pop/static_cast<double>(old_region_dval) - target
    )/target;

    // now compute population penalties
    double pop_pen1 = std::sqrt(static_cast<double>(N) - 2) * std::log(1e-12 + dev1);
    double pop_pen2 = std::sqrt(static_cast<double>(N) - 2) * std::log(1e-12 + dev2);
    double old_pop_pen = std::sqrt(static_cast<double>(N) - 2) * std::log(1e-12 + old_dev);

    // now return the values for the old region minus the two new ones
    return old_pop_pen * pop_temper - (pop_pen1 + pop_pen2) * pop_temper;
}




//' Compute the optimal log incremental weight of a plan
//'
//' Given a plan object this computes the minimum variance weights as derived in
//' <PAPER NAME HERE>. This is equal to the inverse of a sum over all
//' adjacent regions in a plan if using generalized region split and 
//' sum over all districts adajacent to the remainder if using one district
//' split approach.
//'
//' @title Compute Optimal Incremental Weight of a plan
//'
//' @param g The underlying map graph
//' @param plan A plan object
//' @param split_district_only whether or not to compute the weights under 
//' the district only split scheme or not.
//' @param target The target population for a single district
//' @param pop_temper The population tempering parameter
//'
//' @details No modifications to inputs made
//'
//' @return the log of the incremental weight of the plan
//'
double compute_optimal_log_incremental_weight(
        const Graph &g, const Plan &plan,
        bool split_district_only,
        const double target, const double pop_temper){
    // get a region level graph
    Graph rg = get_region_graph(g, plan);

    // sanity check: make sure the region level graph's number of vertices is
    // the same as the number of regions in the plan
    if((int)rg.size() != plan.num_regions){
        Rprintf("SOMETHING WENT WRONG IN COMPPUTE LOG INCREMENT WEIGHJT\n");
    }

    // Now get all unique adjacent region pairs and store them such that
    // the first vertex is always the smaller numbered of the edge
    // We use a set because it doesn't keep duplicates
    std::set<std::pair<int, int>> adj_region_pairs;


    // If n < N and split_district_only true then just get adjacent to remainder
    if(plan.num_regions < plan.N && split_district_only){
        // get remainder region id
        int u = plan.remainder_region;

        // Iterate through the adjacency list of remainder region
        for(auto v: rg.at(u)){
            // store the edges such that smaller vertex is always first
            if (u < v) {
                adj_region_pairs.emplace(u, v);
            } else {
                adj_region_pairs.emplace(v, u);
            }
        }
    }else{ // else get all adjacent districts
        // Iterate through the adjacency list
        // Now get all unique adjacent region pairs and store them such that
        // the first vertex is always the smaller numbered of the edge
        for (int u = 0; u < rg.size(); u++){
            // iterate over neighbors
            for(auto v: rg.at(u)){
                // store the edges such that smaller vertex is always first
                if (u < v) {
                    adj_region_pairs.emplace(u, v);
                } else {
                    adj_region_pairs.emplace(v, u);
                }
            }
        }
    }


    double incremental_weight = 0.0;

    bool do_pop_temper = (plan.num_regions < plan.N) && pop_temper > 0;


    // Now iterate over all adjacent pairs
    for (const auto& edge : adj_region_pairs) {
        // get log of boundary length between regions
        double log_boundary =  region_log_boundary(g, plan.region_ids, edge.first, edge.second);
        // get prob of selecting union of adj regions
        double log_splitting_prob;

        // for one district split the probability that region was chosen to be split is always 1
        if(split_district_only){
            log_splitting_prob = 0;
        }else{
            // in generalized region split find probability you would have 
            // picked to split the union of the the two regions 
            log_splitting_prob = get_log_retroactive_splitting_prob(plan, edge.first, edge.second);
        }


        // Do population tempering term if not final
        double log_temper;
        if(do_pop_temper){
            log_temper = compute_log_pop_temper(
                plan,
                edge.first, edge.second,
                target, pop_temper
            );
        }else{
            log_temper = 0;
        }

        // multiply the boundary length and selection probability by adding the logs
        // now exponentiate and add to the sum
        incremental_weight += std::exp(log_boundary
                                           + log_splitting_prob
                                           + log_temper);
    }


    // Check its not infinity
    if(incremental_weight == -std::numeric_limits<double>::infinity()){
        Rprintf("Error! weight is negative infinity for some reason \n");
    }


    // now return the log of the inverse of the sum
    return -std::log(incremental_weight);

}


//' Computes log unnormalized weights for vector of plans
//'
//' Using the procedure outlined in <PAPER HERE> this function computes the log
//' incremental weights and the unnormalized weights for a vector of plans (which
//' may or may not be the same depending on the parameters).
//'
//' @title Compute Log Unnormalized Weights
//'
//' @param pool A threadpool for multithreading
//' @param g A graph (adjacency list) passed by reference
//' @param plans_vec A vector of plans to compute the log unnormalized weights
//' of
//' @param split_district_only whether or not to compute the weights under 
//' the district only split scheme or not. If `split_district_only` is true
//' then uses optimal weights from one-district split scheme.
//' @param log_incremental_weights A vector of the log incremental weights
//' computed for the plans. The value of `log_incremental_weights[i]` is
//' the log incremental weight for `plans_vec[i]`
//' @param unnormalized_sampling_weights A vector of the unnormalized sampling
//' weights to be used with sampling the `plans_vec` in the next iteration of the
//' algorithm. Depending on the other hyperparameters this may or may not be the
//' same as `exp(log_incremental_weights)`
//' @param target Target population of a single district
//' @param pop_temper <DETAILS NEEDED>
//'
//' @details Modifications
//'    - The `log_incremental_weights` is updated to contain the incremental
//'    weights of the plans
//'    - The `unnormalized_sampling_weights` is updated to contain the unnormalized
//'    sampling weights of the plans for the next round
void get_all_plans_log_gsmc_weights(
        RcppThread::ThreadPool &pool,
        const Graph &g, std::vector<Plan> &plans_vec,
        bool split_district_only,
        std::vector<double> &log_incremental_weights,
        std::vector<double> &unnormalized_sampling_weights,
        double target, double pop_temper
){
    int M = (int) plans_vec.size();

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        double log_incr_weight = compute_optimal_log_incremental_weight(
            g, plans_vec.at(i), split_district_only,
            target, pop_temper);
        log_incremental_weights[i] = log_incr_weight;
        unnormalized_sampling_weights[i] = std::exp(log_incr_weight);
    });

    // Wait for all the threads to finish
    pool.wait();


    return;
}