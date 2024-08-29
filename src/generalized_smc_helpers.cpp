#include "generalized_smc_helpers.h"


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
//' @return sum of weights squared over sum of squared weights (sum(wgt)^2 / sum(wgt^2))
//'
double compute_n_eff(const std::vector<double> &log_wgt) {
    double sum_wgt = 0.0;
    double sum_wgt_squared = 0.0;

    for (const double& log_w : log_wgt) {
        double wgt = std::exp(log_w);
        sum_wgt += wgt;
        sum_wgt_squared += wgt * wgt;
    }


    return std::exp(
        (2 * std::log(sum_wgt)) - std::log(sum_wgt_squared)
    );
}



//' Returns the log of the Count of the number of edges across two regions in
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
//' @param plan A plan object
//' @param region1_label The label of the first region
//' @param region2_label The label of the second region
//'
//' @details Modifications
//'    - If two new valid regions are split then the plan object is updated accordingly
//'    - If two new valid regions are split then the new_region_labels is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'
//' @return True if two valid regions were split off false otherwise
//'
 double region_log_boundary(const Graph &g, const Plan &plan,
                            std::string const&region1_label,
                            std::string const&region2_label
 ) {
     int V = g.size();

     double count = 0; // count of number of edges across the two regions
     // Iterate over every vertex
     for (int i = 0; i < V; i++) {
         /* Check vertex if in first region, if not continue
          Since edges appear twice in adjcancy list to avoid double
          counting we will only count those where first edge is in
          region 1
          */
         if (plan.region_labels.at(i) != region1_label) continue;


         // Get vertices neighbors
         std::vector<int> nbors = g[i];

         // Since edges show up twice to avoid double counting we will only count
         // edges where first one is region1_label and second is region2_label



         // Now check if neighbors are in second regoin
         for (int nbor : nbors) {
             if (plan.region_labels.at(nbor) != region2_label)
                 continue;
             // otherwise, boundary with root -> ... -> i -> nbor
             count += 1.0;
         }
     }

     return std::log(count);
 }




//' Selects a multidistrict with probability proportional to its d_nk value and
//' returns the log probability of the selected region
//'
//' Given a plan object with at least one multidistrict this function randomly
//' selects a multidistrict with probability proporitional to its d_nk value
//' (relative to all multidistricts) and returns the log of the probability that
//' region was chosen.
//'
//'
//' @title Choose multidistrict to split
//'
//' @param plan A plan object
//' @param region_to_split a string that will be updated by reference with the
//' name of the region selected to split
//'
//' @details Modifications
//'    - sets `region_to_split` to be the region that was selected
//'
//' @return the log of the probability the specific value of `region_to_split` was chosen
//'
double choose_multidistrict_to_split(
        Plan const&plan, std::string &region_to_split){

    if(plan.num_multidistricts < 1){
        Rprintf("ERROR: Trying to find multidistrict to split when there are none!\n");
    }

    // count total
    int total_multi_ds = 0;

    // make vectors with cumulative d value and region label for later
    std::vector<int> multi_d_vals;
    std::vector<std::string> region_labels;

    // Iterate over all regions
    for (const auto &[key, value]: plan.r_to_d_map ) {

        // collect info if multidistrict
        if(value > 1){
            // Add that regions d value to the total
            total_multi_ds += value;
            // add the count and label to vector
            multi_d_vals.push_back(value);
            region_labels.push_back(key);
        }
    }


    // Now pick an index proportational to d_nk value
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(multi_d_vals.begin(), multi_d_vals.end());

    int idx = d(gen);

    region_to_split = region_labels[idx];
    double log_prob = std::log(
        static_cast<double>(multi_d_vals[idx])
        ) - std::log(
                static_cast<double>(total_multi_ds)
        );

    return log_prob;
}





/*
 * Make the region adjacency graph for `plan` from the overall precinct graph `g`
 * where the regions are represented by their integer id representation form
 */
// NOT TESTED but taken from district_graph function which was tested
Graph get_region_graph(const Graph &g, const Plan &plan) {
    int V = g.size();
    Graph out;

    std::vector<std::vector<bool>> gr_bool(
            plan.num_regions, std::vector<bool>(plan.num_regions, false)
        );


    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        // Find out which region this vertex corresponds to
        int region_num_i = plan.region_num_ids[i];

        // now iterate over its neighbors
        for (int nbor : nbors) {
            // find which region neighbor corresponds to
            int region_num_j = plan.region_num_ids[nbor];
            // if they are different regions mark matrix true
            if (region_num_i != region_num_j) {
                gr_bool.at(region_num_i).at(region_num_j) = true;
            }
        }

    }


    for (int i = 0; i < plan.num_regions; i++) {
        std::vector<int> tmp;
        for (int j = 0; j < plan.num_regions; j++) {
            if (gr_bool.at(i).at(j)) {
                tmp.push_back(j);
            }
        }
        out.push_back(tmp);
    }

    return out;
}

// Given a plan and two regions get the probability their log 'union' would have
// been selected to be split
double get_log_retroactive_splitting_prob(
        const Plan &plan,
        const std::string region1_label, const std::string region2_label
){
    // We just want

    // count total
    int total_multi_ds = 0;


    // Iterate over all regions
    for (const auto &[key, value]: plan.r_to_d_map ) {
        // collect info if multidistrict and not the two we started with
        if(value > 1 && key != region1_label && key != region2_label){
            // Add that regions d value to the total
            total_multi_ds += value;
        }
    }

    // Now get the sum of dnk values of two regions
    int unioned_region_dnk = plan.r_to_d_map.at(region1_label) + plan.r_to_d_map.at(region2_label);
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


// Compute the incremental log weight
double compute_log_incremental_weight(const Graph &g, const Plan &plan){
    Graph rg = get_region_graph(g, plan);


    // check
    if((int)rg.size() != plan.num_regions){
        Rprintf("SOMETHING WENT WRONG IN COMPPUTE LOG INCREMENT WEIGHJT\n");
    }

    // First get all adjacent region pairs

    // Set to store unique pairs of adjacent vertices
    std::set<std::pair<int, int>> adj_region_pairs;
    // Iterate through the adjacency list
    for (int u = 0; u < rg.size(); u++){
        // iterate over neighbors
        for(auto v: rg[u]){
            if (u < v) {
                adj_region_pairs.emplace(u, v);
            } else {
                adj_region_pairs.emplace(v, u);
            }
        }
    }



    double incremental_weight = 0.0;


    // Now iterate over all pairs
    for (const auto& edge : adj_region_pairs) {
        // convert to region string label representation
        std::string region1_label = plan.num_id_to_str_label_map.at(edge.first);
        std::string region2_label = plan.num_id_to_str_label_map.at(edge.second);
        // get log of boundary
        double log_boundary =  region_log_boundary(g, plan, region1_label, region2_label);
        // double log_boundary = 0;
        // get prob of selecting union of adj regions
        double log_splitting_prob = get_log_retroactive_splitting_prob(plan, region1_label, region2_label);
        // NO e^J TERM YET
        // add the logs, exponentiate and add the sum
        incremental_weight += std::exp(log_boundary + log_splitting_prob);
    }


    // Check its not infinity
    if(incremental_weight == -std::numeric_limits<double>::infinity()){
        Rprintf("Error! weight is negative infinity for some reason \n");
    }


    // now return log of the weight
    return -std::log(incremental_weight);

}
