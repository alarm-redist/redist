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
double compute_n_eff(const arma::subview_col<double> log_wgt) {
    double sum_wgt = 0.0;
    double sum_wgt_squared = 0.0;

    // compute sum of squares and square of sum
    for (const double& log_w : log_wgt) {
        double wgt = std::exp(log_w);
        sum_wgt += wgt;
        sum_wgt_squared += std::exp(2*log_w);
    }


    return std::exp(
        (2 * std::log(sum_wgt)) - std::log(sum_wgt_squared)
    );
}



// only gets regions adjacent to indices where valid_regions is true
void get_all_adj_pairs(
    Graph const &g, std::vector<std::pair<int, int>> &adj_pairs_vec,
    arma::subview_col<arma::uword> const &vertex_region_ids,
    std::vector<bool> const valid_regions
){
    int V = g.size();

    // Set to put adjacent pairs in 
    std::set<std::pair<int, int>> adj_region_pairs;

    // iterate over all vertices in g
    for (int i = 0; i < V; i++) {
        // Find out which region this vertex corresponds to
        int region_num_i = vertex_region_ids(i);

        // check if its a valid region and if not continue
        if(!valid_regions.at(region_num_i)){
            continue;
        }

        std::vector<int> nbors = g.at(i);

        // now iterate over its neighbors
        for (int nbor : nbors) {
            // find which region neighbor corresponds to
            int region_num_j = vertex_region_ids(nbor);

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



//' Returns a unordered_map mapping valid pairs of adjacent regions to 
//' the length of the boundary between them
//'
//' Finds all pairs of adjacent regions where at least one region is true in 
//' the `check_adj_to_regions` parameter and returns a hash map mapping the pairs
//' of valid adjacent regions to the length of the boundary between them.
//'
//'
//' @title Get Valid Adjacent Regions to Boundary Lengths Map
//'
//' @param g A graph (adjacency list) passed by reference
//' @param plan A plan object
//' @param check_adj_to_regions A vector of length plan.num_regions of regions 
//' tracking whether or not we should check for edges adjacent to vertices in
//' that region. For example, `check_adj_to_regions[i] == true` means we will 
//' attempt to find edges adjacent to any vertex in region i
//' If all values are true then the function will find all pairs of adjacent regions
//'
//' @details No modifications to inputs made
//'
//' @return A hash map mapping (std::pair<int,int> to int) that maps a pair of region ids
//' in sorted order (ie smallest first) to the length of the boundary between them in `g`
//'
std::unordered_map<std::pair<int, int>, int, bounded_hash> get_pairs_and_boundary_len_map(
    Graph const &g, Plan const &plan,
    std::vector<bool> const &check_adj_to_regions
){
    // NOTE: In the case where its one district split you maybe don't even need
    // a hash map since you can just index by the not remainder region but nbd
    // for now
    
    // Initialize unordered_map with num_region * 2.5 buckets
    // Hueristic. Bc we know planar graph has at most 3|V| - 6 edges
    int init_bucket_size = std::ceil(2.5*plan.num_regions);

    // create the hash map
    std::unordered_map<std::pair<int, int>, int, bounded_hash> region_pair_map(
        init_bucket_size, bounded_hash(plan.num_regions-1)
        );

    int V = g.size();


    for (int i = 0; i < V; i++) {
        // Find out which region this vertex corresponds to
        int region_num_i = plan.region_ids(i);

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!check_adj_to_regions.at(region_num_i)){
            continue;
        }

        // get neighbors
        std::vector<int> nbors = g.at(i);

        // now iterate over its neighbors
        for (int nbor : nbors) {
            // find which region neighbor corresponds to
            int region_num_j = plan.region_ids(nbor);

            // if they are different regions then region i and j are adjacent 
            // as they share an edge across
            if (region_num_i != region_num_j) {
                // if region j is invalid then we don't need to worry about double counting so just add it
                if(!check_adj_to_regions.at(region_num_j)){
                    region_pair_map[
                        {std::min(region_num_j,region_num_i), std::max(region_num_j,region_num_i)}
                        ]++;
                }else{ // else if both valid then edge will appear twice so we only count if region i
                // smaller
                    if (region_num_i < region_num_j) {
                    region_pair_map[{region_num_i, region_num_j}]++;
                    }
                }
            }
        }
    }

    // adj_pairs_vec
    return region_pair_map;
}


//' Returns a vector of the triple (smaller region id, bigger region id, boundary len)
//' for all valid pairs of adjacent regions in the plan. (Either all adjacent regions if
//' doing generalized region splits or just adjacent to the remainder if only doing 
//' one district splits.)
//'
//'
//' @title Get All Valid Adjacent Regions and their Boundary Length
//'
//' @param g A graph (adjacency list) passed by reference
//' @param plan A plan object
//' @param split_district_only If true only gets regions adjacent to the remainder but if 
//' false then gets all adjacent regions in the plan
//'
//' @details No modifications to inputs made
//'
//' @return A vector of integer arrays of size 3 where the values are
//' (smaller region id, bigger region id, boundary len)
//'
std::vector<std::array<int, 3>> get_valid_adj_regions_and_boundary_lens_vec(
    Graph const &g, Plan const &plan,
    bool const split_district_only 
){

    // make vector for if we just want adjacent to remainder or all
    // adj pairs
    std::vector<bool> valid_regions;
    if(split_district_only && plan.num_regions < plan.ndists){
        // if splitting district only then only find adjacent to remainder
        // which is region 2 because its the bigger one
        valid_regions.resize(plan.num_regions, false);
        valid_regions.at(plan.remainder_region) = true;
    }else{
        valid_regions.resize(plan.num_regions, true);
    }

    // Get all valid adj pairs and the boundary length
    auto region_pair_map = get_pairs_and_boundary_len_map(
    g, plan, valid_regions);

    // Now put them into a vector 
    std::vector<std::array<int, 3>> adj_pairs_and_boundary_lens;
    adj_pairs_and_boundary_lens.reserve(region_pair_map.size());

    for (const auto& entry : region_pair_map) {
        const auto region_pair = entry.first; // get the region pair 
        const int boundary_len = entry.second; // get the boundary length

        

        // Make it so elements of array are 
        // - smaller region id
        // - bigger region id
        // - boundary length
        adj_pairs_and_boundary_lens.emplace_back(
            std::array<int, 3>{region_pair.first, region_pair.second, boundary_len});
    }

    return adj_pairs_and_boundary_lens;

}



//' Returns a vector of the triple (smaller region id, bigger region id, boundary len)
//' for all valid pairs of adjacent regions in the plan. (Either all adjacent regions if
//' doing generalized region splits or just adjacent to the remainder if only doing 
//' one district splits.)
//'
//'
//' @title Get All Valid Adjacent Regions and their Boundary Length
//'
//' @param g A graph (adjacency list) passed by reference
//' @param plan A plan object
//' @param split_district_only If true only gets regions adjacent to the remainder but if 
//' false then gets all adjacent regions in the plan
//'
//' @details No modifications to inputs made
//'
//' @return A vector of integer arrays of size 3 where the values are
//' (smaller region id, bigger region id, boundary len)
//'
std::vector<std::array<int, 3>> get_valid_adj_regions_and_boundary_lens(
    Graph const &g, Plan const &plan,
    bool const split_district_only,
    std::vector<bool> const &check_adj_to_regions,
    std::vector<std::vector<bool>> const valid_merge_pairs
){

    // Get all valid adj pairs and the boundary length
    auto region_pair_map = get_pairs_and_boundary_len_map(
    g, plan, check_adj_to_regions);

    // Now put them into a vector 
    std::vector<std::array<int, 3>> adj_pairs_and_boundary_lens;
    adj_pairs_and_boundary_lens.reserve(region_pair_map.size());

    for (const auto& entry : region_pair_map) {
        const auto region_pair = entry.first; // get the region pair 
        const int boundary_len = entry.second; // get the boundary length

    
        // Make it so elements of array are 
        // - smaller region id
        // - bigger region id
        // - boundary length
        adj_pairs_and_boundary_lens.emplace_back(
            std::array<int, 3>{region_pair.first, region_pair.second, boundary_len});
    }

    return adj_pairs_and_boundary_lens;

}



// Returns a sampler over a vector of adjacent pairs where the probability 
// of a pair is decided according to `selection_type`
//'
//' Current supported options are
//'     - uniform - Every pair has equal probability
//'     - district_pair - double district pairs have weight 1000, one district is 10,
//'         and two multidistricts have 1/(1+sum of their dvals)
//'
//' @title Get Sampler over Adj Regions List
//'
//' @param plan A plan object
//' @param adj_pairs_and_boundary_lens A vector where each pair is 
//' (adj region1, adj region2, boundary length between 2 regions)
//' @param selection_type A string controlling the function to use
//' in assigning the unnormalized weight to each pair
//'
//' @details No modifications to inputs made
//'
//' @return A sampler where index i has probability proportional to the weight 
//' given to that pair 
//'
std::discrete_distribution<>  get_adj_pair_sampler(
    Plan const &plan,
    std::vector<std::array<int, 3>> const &adj_pairs_and_boundary_lens,
    std::string const &selection_type
){
    // make a vector for the unnormalized weights 
    std::vector<double> unnormalized_sampling_weights(adj_pairs_and_boundary_lens.size());

    // if uniform then just make all the weights 1 over the size 
    if(selection_type == "uniform"){
        std::fill(
            unnormalized_sampling_weights.begin(), 
            unnormalized_sampling_weights.end(), 
            static_cast<double>(adj_pairs_and_boundary_lens.size())
            );
    }else if(selection_type == "district_pair"){
        // if district pair then give huge weight to pairs of districts 
        for (size_t i = 0; i < adj_pairs_and_boundary_lens.size(); i++)
        {
            int region1_id = adj_pairs_and_boundary_lens.at(i).at(0);
            int region2_id = adj_pairs_and_boundary_lens.at(i).at(1);
            int region1_dval = plan.region_sizes(region1_id); int region2_dval = plan.region_sizes(region2_id);

            // check if both are districts
            if(region1_dval == 1 && region2_dval == 1){
                unnormalized_sampling_weights.at(i) = 1000.0;
            }else if(region1_dval == 1 || region2_dval == 1){
                // else check if at least one is a district
                unnormalized_sampling_weights.at(i) = 10.0;
            }else{
                // if both multidistricts then do 1/1+(sum of two dvals)
                // This penalizes bigger pairs 
                unnormalized_sampling_weights.at(i) = 1/(1+ static_cast<double>(region1_dval+region2_dval));
            }
        }
        
    }else if(selection_type == "multidistrict_pair"){
        // if multidistrict pair then give huge weight to pairs of districts 
        for (size_t i = 0; i < adj_pairs_and_boundary_lens.size(); i++)
        {
            int region1_id = adj_pairs_and_boundary_lens.at(i).at(0);
            int region2_id = adj_pairs_and_boundary_lens.at(i).at(1);
            int region1_dval = plan.region_sizes(region1_id); int region2_dval = plan.region_sizes(region2_id);

            // check if both are districts
            if(region1_dval == 1 && region2_dval == 1){
                unnormalized_sampling_weights.at(i) = .5;
            }else if(region1_dval == 1 || region2_dval == 1){
                // else check if at least one is a district
                unnormalized_sampling_weights.at(i) = 1.0;
            }else{
                // if both multidistricts then do 1/1+(sum of two dvals)
                // This penalizes bigger pairs 
                unnormalized_sampling_weights.at(i) = 1000.0;
            }
        }
        
    }else{
        throw Rcpp::exception("No valid adjacent pair sampler provided!");
    }

    // now create the sampler 
    std::discrete_distribution<> index_sampler(
            unnormalized_sampling_weights.begin(),
            unnormalized_sampling_weights.end()
        );

    // Return 1 over the log of the size 
    return index_sampler;
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
                            const arma::subview_col<arma::uword> &vertex_region_ids,
                            int const region1_id,
                            int const region2_id
){
     int V = g.size(); // get number of vertices
     double count = 0; // count of number of edges across the two regions

     // Iterate over every vertex in the graph
     for (int i = 0; i < V; i++) {
         /* Check vertex if in first region, if not continue
          Since edges appear twice in adjacency list to avoid double
          counting we will only count those where first edge is in
          region 1
          */
         if (vertex_region_ids(i) != region1_id) continue;

         // Get vertice's neighbors
         std::vector<int> nbors = g[i];

         // Since edges show up twice to avoid double counting we will only count
         // edges where first one is region1_id and second is region2_id

         // Now check if neighbors are in second region
         for (int nbor : nbors) {
             if (vertex_region_ids(nbor) != region2_id)
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
    const int old_region_boundary_length, const int new_region_boundary_length, 
    const double current_pair_merge_prob, const double new_pair_merge_prob, 
    Plan &current_plan, Plan &new_plan
){
    // get the log boundary lengths
    double old_log_boundary_length = std::log(static_cast<double>(old_region_boundary_length));

    double new_log_boundary_length = std::log(static_cast<double>(new_region_boundary_length));

    // Get the merge probability which for now is uniform
    double log_old_merge_prob = std::log(current_pair_merge_prob);

    double log_new_merge_prob = std::log(new_pair_merge_prob);

    // double log_mh_ratio_numerator = old_log_boundary_length + log_old_merge_prob;
    // double log_mh_ratio_denominator = new_log_boundary_length + log_new_merge_prob;

    double log_mh_ratio_numerator = old_log_boundary_length + log_new_merge_prob;
    double log_mh_ratio_denominator = new_log_boundary_length + log_old_merge_prob;

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
        const std::vector<bool> &valid_region_sizes_to_split,
        const int region1_id, const int region2_id
){

    // count total dvals of all the multidistricts
    int total_multi_ds = 0;

    if(DEBUG_VERBOSE) Rprintf("Possible options: ");

    // Iterate over all regions
    for(int region_id = 0; region_id < plan.num_regions; region_id++) {
        int region_size = plan.region_sizes(region_id);
        // ignore if this is not a region size that could be split 
        if(!valid_region_sizes_to_split[region_size]) continue;

        // collect info if multidistrict and not the two we started with
        if(region_size > 1 && region_id != region1_id && region_id != region2_id){
            // Add that regions d value to the total
            if (DEBUG_VERBOSE) Rprintf(" %d ", region_size);
            total_multi_ds += region_size;
        }
    }

    // Now get the sum of dnk values of two regions aka the unioned old region
    int unioned_region_dnk = plan.region_sizes(region1_id) + plan.region_sizes(region2_id);

    if (DEBUG_VERBOSE) Rprintf(" (unioned) %d\n", unioned_region_dnk);

    // update the total number of multi district dvals with the value the union region would have been
    total_multi_ds += unioned_region_dnk;

    

    // so prob of picking is the sum of the two regions dnk over the dnk of all
    // multidistricts

    double log_prob = std::log(
        static_cast<double>(unioned_region_dnk)
    ) - std::log(
            static_cast<double>(total_multi_ds)
    );

    // Rprintf("For regions (%d,%d), size (%d,%d) - Prob is %d/%d, so log is %.5f\n",
    // region1_id, region2_id, 
    // (int) plan.region_sizes(region1_id), (int) plan.region_sizes(region2_id),
    // unioned_region_dnk, total_multi_ds,
    // log_prob);

    return log_prob;

}




// computes the backwards kernel that is uniform in 
// the number of ancestors 
double compute_uniform_graph_plan_adj_log_incremental_weight(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function, Plan const &plan, 
    const TreeSplitter &edge_splitter, bool compute_log_splitting_prob,
    bool is_final_plan
){


    // find the index of the two newly split regions 
    auto most_recent_split_regions = plan.get_most_recently_split_regions();
    int region1_id = most_recent_split_regions.first;
    int region2_id = most_recent_split_regions.second;


    // Find all pairs and boundary length 
    // IN FUTURE CAN REPLACE WITH count valid adj regions and log oundary
    // that is two passes through as opposed to this which requires a hash as well
    auto region_pairs_and_log_boundary_len = plan.get_valid_adj_regions_and_eff_log_boundary_lens(
        map_params, splitting_schedule,
        edge_splitter
    );


    const double log_num_valid_adj_regions = std::log(
        static_cast<double>(region_pairs_and_log_boundary_len.size())
    );

    double log_boundary_len; 
    // need to iterate over the map until we find the pair 
    for (auto const &pair_tuple: region_pairs_and_log_boundary_len)
    {
        if(std::get<0>(pair_tuple) == region1_id && std::get<1>(pair_tuple) == region2_id){
            log_boundary_len = std::get<2>(pair_tuple);
        }
    }

    double log_splitting_prob = 0;
    // for one district split the probability that region was chosen to be split is always 1
    if(compute_log_splitting_prob){
        // in generalized region split find probability you would have 
        // picked to split the union of the the two regions 
        log_splitting_prob = get_log_retroactive_splitting_prob(plan, 
        splitting_schedule.valid_region_sizes_to_split, region1_id, region2_id);
        if(DEBUG_VERBOSE) Rprintf("Computed split prob %f\n", std::exp(log_splitting_prob));
    }

    double log_constraint_ratio = 0;

    double log_region1_score, log_region2_score, log_merged_region_score;

    log_region1_score = log_region2_score = log_merged_region_score = 0;

    // compute if any constraints 
    if(scoring_function.any_constraints){
        // compute scoring functions
        log_region1_score = scoring_function.compute_log_region_score(
            plan, region1_id, is_final_plan
        );
        log_region2_score = scoring_function.compute_log_region_score(
            plan, region2_id, is_final_plan
        );
        log_merged_region_score = scoring_function.compute_log_merged_region_score(
            plan, region1_id, region2_id, is_final_plan
        );
    }



    // The weight is 
    //      - Numerator: e^J(new_region1+new_region2)
    //      - Denominator: multid_selection_prob * boundary_len * num_valid_adj_regions * e^J(old_region)

    const double log_numerator = log_region1_score + log_region2_score;
    const double log_denom = log_num_valid_adj_regions + log_boundary_len + log_merged_region_score + log_splitting_prob;

    int region1_size = plan.region_sizes(region1_id);
    int region2_size = plan.region_sizes(region2_id);
    if(DEBUG_VERBOSE){
    Rprintf("Doing (%d,%d) - sizes (%d, %d): len %f, split prob %f, ratio %f!\n", 
        region1_id, region2_id, region1_size, region2_size, std::exp(log_boundary_len), 
        std::exp(log_splitting_prob), std::exp(log_numerator - log_merged_region_score));

    Rprintf("Numerator - %f * %f\n", 
        std::exp(log_region1_score), std::exp(log_region2_score)
    );

    Rprintf("Denominator - %f * %f * %f * %f \n", 
        std::exp(log_num_valid_adj_regions), std::exp(log_boundary_len),
        std::exp(log_merged_region_score), std::exp(log_splitting_prob)
    );}

    // multiply the boundary length and selection probability by adding the logs
    // now exponentiate and add to the sum
    double incremental_weight = log_numerator - log_denom;

    // Check its not - infinity
    if(incremental_weight == -std::numeric_limits<double>::infinity()){
        throw Rcpp::exception("Error! weight is negative infinity for some reason \n");
    }

    return incremental_weight;
}


void get_all_plans_uniform_adj_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    const std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    std::vector<double> &unnormalized_sampling_weights,
    int verbosity
){
    int M = (int) plans_ptr_vec.size();

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {

        double log_incr_weight = compute_uniform_graph_plan_adj_log_incremental_weight(
            map_params, splitting_schedule,
            scoring_function, 
            *plans_ptr_vec.at(i), 
            *tree_splitters_ptr_vec.at(i),
            compute_log_splitting_prob,
            is_final_plans
        );

        log_incremental_weights(i) = log_incr_weight;
        unnormalized_sampling_weights[i] = std::exp(log_incr_weight);
    });

    // Wait for all the threads to finish
    pool.wait();


    return;
}




// NOT FULLY TESTED but taken from district_graph function which was tested
//' Creates the region level graph of a plan
//'
//' Given a plan object this returns a graph of the regions in the plan using
//' the region ids as indices
//'
//' @title Get Region-Level Graph
//'
//' @param g The graph of the entire map
//' @param plan A plan object
//'
//' @details No modifications to inputs made
//'
//' @return the log of the probability the specific value of `region_to_split` was chosen
//'
Graph get_region_graph(const Graph &g, const Plan &plan) {
    int V = g.size();
    Graph out;

    // make a matrix where entry i,j represents if region i and j are adjacent
    std::vector<std::vector<bool>> gr_bool(
            plan.num_regions, std::vector<bool>(plan.num_regions, false)
        );

    // iterate over all vertices in g
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        // Find out which region this vertex corresponds to
        int region_num_i = plan.region_ids[i];

        // now iterate over its neighbors
        for (int nbor : nbors) {
            // find which region neighbor corresponds to
            int region_num_j = plan.region_ids[nbor];
            // if they are different regions mark matrix true since region i
            // and region j are adjacent as they share an edge across
            if (region_num_i != region_num_j) {
                gr_bool.at(region_num_i).at(region_num_j) = true;
            }
        }
    }


    // Now build the region level graph
    for (int i = 0; i < plan.num_regions; i++) {
        // create vector of i's neighbors
        std::vector<int> tmp;
        for (int j = 0; j < plan.num_regions; j++) {
            // check if i and j are adjacent, if so add j
            if (gr_bool.at(i).at(j)) {
                tmp.push_back(j);
            }
        }
        out.push_back(tmp);
    }

    return out;
}




double get_log_retroactive_splitting_prob_for_joined_tree(
    MapParams const &map_params,
    Plan const &plan, Tree &ust, const TreeSplitter &edge_splitter,
    std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
    const int region1_root, const int region2_root,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try
){
    int region1_id = plan.region_ids(region1_root);
    int region2_id = plan.region_ids(region2_root);
    int total_merged_region_pop = plan.region_pops.at(region1_id)+plan.region_pops.at(region2_id); 
    int total_merged_region_size = plan.region_sizes(region1_id)+plan.region_sizes(region2_id);

    // Get all the valid edges in the joined tree 
    std::vector<EdgeCut> valid_edges = get_valid_edges_in_joined_tree(
        map_params, plan.forest_graph, ust,
        visited, pops_below_vertex,
        region1_id, region1_root,
        region2_id, region2_root,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_pop, total_merged_region_size
    );

    // find the index of the actual edge we cut 
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(
        region1_root, region2_root, region1_root, 
        plan.region_sizes(region2_id), plan.region_pops.at(region2_id),
        plan.region_sizes(region1_id), plan.region_pops.at(region1_id)
    );

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    int actual_cut_edge_index = std::distance(valid_edges.begin(), it);
    // REprintf("Actual Cut Edge at Index %d\n", actual_cut_edge_index);

    return edge_splitter.get_log_selection_prob(map_params, valid_edges, actual_cut_edge_index);
}





// eventually need to modify to allow presaved options
// OLD DOCUMENTATION FROM GRAPH THING NEED TO UPDATE
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
double compute_log_optimal_weights(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function, Plan const &plan, 
    const TreeSplitter &edge_splitter, bool compute_log_splitting_prob,
    bool is_final_plan
){
    // boolean for whether or not to compute the splitting probability of merged regions
    // dont need to do when

    double incremental_weight = 0.0;

    // get region pair to effective boundary length map
    auto region_pair_log_eff_boundary_map = plan.get_valid_adj_regions_and_eff_log_boundary_lens(
        map_params, splitting_schedule, edge_splitter
    );

    // iterate over the pairs 

    // Now iterate over adjacent region pairs and add splitting and pop temper
    for (const auto& pair_tuple: region_pair_log_eff_boundary_map){
        const int region1_id = std::get<0>(pair_tuple); // get the smaller region id  
        const int region2_id = std::get<1>(pair_tuple); // get the bigger region id  
        const double eff_log_boundary_len = std::get<2>(pair_tuple); // get the effective boundary length

        double log_splitting_prob = 0.0;
        // If not compute_log_splitting_prob then its just log(1) = 0
        if(compute_log_splitting_prob){
            // in generalized region split find probability you would have 
            // picked to split the union of the the two regions 
            log_splitting_prob = get_log_retroactive_splitting_prob(
                plan, splitting_schedule.valid_region_sizes_to_split, 
                region1_id, region2_id
            );
        }

        // add constraint ratio
        double log_constraint_ratio = 0;

        // compute if any constraints 
        if(scoring_function.any_constraints){
        // compute scoring functions
        const double log_region1_score = scoring_function.compute_log_region_score(
            plan, region1_id, is_final_plan
        );
        const double log_region2_score = scoring_function.compute_log_region_score(
            plan, region2_id, is_final_plan
        );
        const double log_merged_region_score = scoring_function.compute_log_merged_region_score(
            plan, region1_id, region2_id, is_final_plan
        );
        // log ratio is log score merged - (log region1 + log region2)
        log_constraint_ratio = log_merged_region_score - (log_region1_score + log_region2_score);
        }

        // int region1_size = plan.region_sizes(region1_id);
        // int region2_size = plan.region_sizes(region2_id);
        // Rprintf("Adding (%d,%d) - sizes (%d, %d): len %f, split prob %f, ratio %f!\n", 
        //     region1_id, region2_id, region1_size, region2_size, std::exp(eff_log_boundary_len), 
        //     std::exp(log_splitting_prob), std::exp(log_constraint_ratio));


        // Now the term is 
        //      - log_splitting_prob
        //      - eff_log_boundary_len
        //      -log_constraint_ratio
        incremental_weight += std::exp(log_splitting_prob + eff_log_boundary_len + log_constraint_ratio);

    }

    // Check its not infinity
    if(incremental_weight == -std::numeric_limits<double>::infinity()){
        Rcpp::stop("Error! weight is negative infinity for some reason \n");
    }

    // Rprintf("Total Boundary %d and log incre = %f\n",
    // total_boundary_len, -std::log(incremental_weight));

    // now return the log of the inverse of the sum
    return -std::log(incremental_weight);

}





//' NEED TO UPDATE THIS IS OLD DOCUMENTATION FOR GRAPH STUFF
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
//' @param plans_ptr_vec A vector of plans to compute the log unnormalized weights
//' of
//' @param split_district_only whether or not to compute the weights under 
//' the district only split scheme or not. If `split_district_only` is true
//' then uses optimal weights from one-district split scheme.
//' @param log_incremental_weights A vector of the log incremental weights
//' computed for the plans. The value of `log_incremental_weights[i]` is
//' the log incremental weight for `plans_ptr_vec[i]`
//' @param unnormalized_sampling_weights A vector of the unnormalized sampling
//' weights to be used with sampling the `plans_ptr_vec` in the next iteration of the
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
void compute_all_plans_log_optimal_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    const std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    std::vector<double> &unnormalized_sampling_weights,
    int verbosity
){
    const int nsims = static_cast<int>(plans_ptr_vec.size());
    const int check_int = 50; // check for interrupts every _ iterations

    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&] (int i) {
        // REprintf("I=%d\n", i);
        double log_incr_weight = compute_log_optimal_weights(
            map_params, splitting_schedule,
            scoring_function, 
            *plans_ptr_vec.at(i), 
            *tree_splitters_ptr_vec.at(i),
            compute_log_splitting_prob,
            is_final_plans
        );

        log_incremental_weights(i) = log_incr_weight;
        unnormalized_sampling_weights[i] = std::exp(log_incr_weight);

        if (verbosity >= 3) {
            ++bar;
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();


    return;
}