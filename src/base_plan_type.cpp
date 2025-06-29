/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/01
 * Purpose: Implements base Plan (ie graph partition) class 
 ********************************************************/

#include "base_plan_type.h"

bool DEBUG_BASE_PLANS_VERBOSE = false;

// checks the inputted plan has the number of regions it claims it does
// checks the sizes and that the labels make sense.
// if makes sense then it counts the number of districts and multidistricts
void Plan::check_inputted_region_sizes(int ndists, bool split_district_only) const{

    // check sum of first num_region elements is ndists and it matches expected
    // number of districts 
    int num_districts_implied_by_sizes_mat = 0;
    int total_size_implied_by_sizes_mat = 0;
    for (size_t i = 0; i < num_regions; i++)
    {
        // make sure each regions dval is non-zero
        if(region_sizes[i] <= 0){
            throw Rcpp::exception("Region size input for region is 0 or less!");
        }

        total_size_implied_by_sizes_mat += region_sizes[i];
        if(region_sizes[i] == 1) num_districts_implied_by_sizes_mat++;

        // add check that if split district only then only last one has size > 1
        if (split_district_only && i != num_regions-1)
        {
            if(region_sizes[i] != 1) throw Rcpp::exception(
                "For partial plan the remainder does not have region id=num_regions-1!"
                ); 
        }
    }


    if(total_size_implied_by_sizes_mat != ndists){
        throw Rcpp::exception("Sum of sizez in region sizes mat does equal ndists!");
    } 

    // check the other entries are zero
    for (int i = num_regions; i < ndists; i++){
        // check that the values are 0,...,num_regions - 1
        if(region_sizes[i] != 0){
            REprintf("Expected Region %d to be zero it was actually %d!\n",
                i, (int) region_sizes[i]);
            throw Rcpp::exception(
                "Plan didn't have the correct expected region size matrix!\n"
                );
        }
    }

    return;
}


void Plan::check_inputted_region_ids(int ndists) const{
    // Use std::unordered_set to store unique elements
    std::unordered_set<int> unique_ids;unique_ids.reserve(num_regions);
    // get unique labels
    for (size_t i = 0; i < region_ids.size(); ++i) {
        unique_ids.insert(region_ids[i]); // Insert each element of the subview column
    }
    int actual_num_regions = static_cast<int>(unique_ids.size());

    // make sure the size is the number of regions 
    if(actual_num_regions != num_regions){
        REprintf("Expected %d regions in plan but there were actually %d!\n",
            num_regions, actual_num_regions);
        throw Rcpp::exception(
            "Plan didn't have the expected number of regions!\n"
            );
    }

    // now check the labels are as expected
    std::vector<int> sorted_labels(unique_ids.begin(), unique_ids.end());
    std::sort(sorted_labels.begin(), sorted_labels.end());
    for (int i = 0; i < num_regions; i++)
    {
        // check that the values are 0,...,num_regions - 1
        if(sorted_labels.at(i) != i){
            REprintf("Expected Region ID %d but it was actually %d!\n",
                i, sorted_labels.at(i));
            throw Rcpp::exception(
                "Plan didn't have the expected region ids!\n"
                );
        }
    }
    
    return;
}



// Constructs existing parital plan
Plan::Plan(int const num_regions,
    const arma::uvec &pop,
    PlanVector &this_plan_region_ids, 
    RegionSizes &this_plan_region_sizes,
    IntPlanAttribute &this_plan_region_pops,
    IntPlanAttribute &this_plan_order_added
):
    region_ids(this_plan_region_ids),
    region_sizes(this_plan_region_sizes),
    region_pops(this_plan_region_pops),
    region_added_order(this_plan_order_added),
    num_regions(num_regions),
    region_order_max(num_regions+2)
{

    // Create other region-level information 
    // fill first num_regions entries with 1,...,num_regions 
    std::iota(
        std::begin(region_added_order), 
        std::begin(region_added_order) + num_regions, 
        1); 
    
    // // compute the population for each of the regions 
    for (size_t v = 0; v < region_ids.size(); v++)
    {
        region_pops[region_ids[v]] += pop(v);
    }
};

// assumes that the inputted attributes are all zero 
Plan::Plan(int const total_seats,
    int const total_pop,
    PlanVector &this_plan_region_ids, 
    RegionSizes &this_plan_region_sizes,
    IntPlanAttribute &this_plan_region_pops,
    IntPlanAttribute &this_plan_order_added
):
    region_ids(this_plan_region_ids),
    region_sizes(this_plan_region_sizes),
    region_pops(this_plan_region_pops),
    region_added_order(this_plan_order_added),
    num_regions(1),
    region_order_max(total_seats + 1)    
{
    region_pops[0] = total_pop;
    region_sizes[0] = total_seats;
    region_added_order[0] = 1;
};



// returns the region ids of the two most recently split regions
// DO NOT CALL THIS WHEN THERE IS ONLY 1 REGION!
std::pair<int, int> Plan::get_most_recently_split_regions() const{
    int largest_index = -1, second_largest_index = -1;
    int largest_value = -1 * num_regions, second_largest_value = -1* num_regions;

    // Iterate through the vector
    for (std::size_t i = 0; i < num_regions; ++i) {
        if (region_added_order[i] > largest_value) {
            // Update second-largest before updating largest
            second_largest_value = largest_value;
            second_largest_index = largest_index;

            largest_value = region_added_order[i];
            largest_index = i;
        } else if (region_added_order[i] > second_largest_value) {
            // Update second-largest only
            second_largest_value = region_added_order[i];
            second_largest_index = i;
        }
    }

    // always make region1 the region with the smaller id
    int region1_id = std::min(largest_index,second_largest_index);
    int region2_id = std::max(largest_index,second_largest_index);

    return std::make_pair(region1_id, region2_id);
}



// Takes a plan and reorders the regions according to the order the regions
// were split
// IT IS VERY IMPORTANT THAT THE TWO PLANS NOT POINT TO THE SAME arma::umat or everything breaks
// This breaks the dummy plan and does not copy it to be the original
void Plan::reorder_plan_by_oldest_split(
    Plan &dummy_plan) {


    // Recall that the region_added_order attribute is a vector that stores info
    // on the relative split order. So if entry i is greater than entry j that 
    // means region i was split more recently than region j

    // For example if you had a plan with 6 regions and a region_added_order vector set to
    // {4,1,7,3,9,2} then that means 
    //  - Region 0 was added 3rd
    //  - Region 1 was added 1st
    //  - Region 2 was added 5th
    //  - Region 3 was added 4rd
    //  - Region 4 was added 6th
    //  - Region 5 was added 2nd 

    // want to return [3, 0, 4, 2, 5, 1]
    // std::vector<int> rel_order = {4,1,7,3,9,2};

    // get the relative order vector of the n regions
    std::vector<int> rel_order(
        this->region_added_order.begin(), 
        this->region_added_order.begin() + this->num_regions
        );


    // Create a vector of indices [0, 1, 2, ..., n-1]
    std::vector<int> indices(rel_order.size());
    std::iota(indices.begin(), indices.end(), 0);

    /*
    Sort the indices based on the values in rel_order so that means 
    indices[i] < indices[j] iff rel_order[i] < rel_order[j] so that means
    indices[0] is the old (ie not updated) label of the oldest split region, 
    In general indices[i] == k  means that the region with the old label of
    k will have a new id of i

    For our example with rel_order = {4,1,7,3,9,2} then that means the indices
    would be sorted to be {1, 5, 3, 0, 2, 4} so we know old region 1 is now zero,
    old region 5 is now 1, etc. Alternatively interpret as indices[i] == k means
    new region i was old region k
    */
    std::sort(indices.begin(), indices.end(), [&rel_order](int a, int b) {
        return rel_order[a] < rel_order[b];
    });


    /*
    Create a vector mapping old region id to the new ordered region id. So 
    `old_region_id_to_new_vec[i]` is the new region id of the region with old
    id i. In other words it means region i should be relabeled as old_region_id_to_new_vec[i]
    */
   
    std::vector<int> old_region_id_to_new_vec(rel_order.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        old_region_id_to_new_vec[indices[i]] = i;
    }

    dummy_plan.region_ids.copy(region_ids);
    // First we relabel all the region vertex ids
    for (size_t i = 0; i < this->region_ids.size(); i++)
    {
        // Send region id i to old_region_id_to_new_vec[i]
        this->region_ids[i] = old_region_id_to_new_vec.at(dummy_plan.region_ids[i]);
    }


    
    // Make some dummy plan attributes a shallow copy of the plan
    dummy_plan.region_sizes.copy(region_sizes);
    dummy_plan.region_pops.copy(region_pops);

    // Now we reorder the region dvals and population 
    for (size_t i = 0; i < this->num_regions; i++)
    {
        // Recall indices[i] == k means the old region with id k now has id i
        // so we want to set the value at the new region id `i` to the value it
        // had in the old one which is `indices[i]`
        int old_region_id = indices[i];
        int new_region_id = i;

        this->region_sizes[new_region_id] = dummy_plan.region_sizes[old_region_id];
        this->region_pops[new_region_id] = dummy_plan.region_pops[old_region_id];
        // Since the regions are in order of added reset the order added to just be 1,..., n
        this->region_added_order[i] = i+1;
    }

    // reset the max region counter 
    this->region_order_max = this->num_regions + 6;

}


// Returns the number of districts and multidistricts in the plan
std::pair<int, int> Plan::get_num_district_and_multidistricts() const{
    int num_districts = 0;
    int num_multidistricts = 0;
    for (size_t region_id = 0; region_id < num_regions; region_id++)
    {
        if(region_sizes[region_id] == 1) num_districts++;
        if(region_sizes[region_id] > 1) num_multidistricts++;
    }

    return std::make_pair(num_districts, num_multidistricts);
    
}

// Prints our object using Rcout. Should be used in Rcpp call
void Plan::Rprint(bool verbose) const{
    auto region_counts = get_num_district_and_multidistricts();
    int num_districts = region_counts.first;
    int num_multidistricts = region_counts.second;
    Rcpp::Rcout << "Plan with " << num_regions << " regions, " << num_districts
                      << " districts, " << num_multidistricts << " multidistricts and "
                      << std::accumulate(region_sizes.begin(), region_sizes.end(), 0) << " total seats and "
                      << region_ids.size() << " Vertices.\n";


    Rcpp::Rcout << "Region Level Values:[";
    for(int region_id = 0; region_id < num_regions; region_id++){
        Rcpp::Rcout << "(Region " << region_id <<
            ", Size=" << static_cast<int>(region_sizes[region_id]) << ", pop= " <<
            region_pops[region_id] <<"), ";
    }
    Rcpp::Rcout << "]\n";

    if(verbose){
        REprintf("Plan IDs: c(");
        for (size_t i = 0; i < region_ids.size(); ++i) {
            REprintf("%u", region_ids[i]);
            if (i + 1 != region_ids.size())
                REprintf(", ");
        }
        REprintf(")\n Plan Sizes: c(");
        for (size_t i = 0; i < num_regions; ++i) {
            REprintf("%u", region_sizes[i]);
            if (i + 1 != num_regions)
                REprintf(", ");
        }
        REprintf(")\n");
    }

    
}


// this shallow copies two plans of the same size
void Plan::shallow_copy(Plan const &plan_to_copy){
    // copy the scalars 
    num_regions = plan_to_copy.num_regions;
    region_order_max = plan_to_copy.region_order_max;
    // copy plan vector
    region_ids.copy(plan_to_copy.region_ids);
    // copy the region sizes vector
    region_sizes.copy(plan_to_copy.region_sizes);
    // copy population
    region_pops.copy(plan_to_copy.region_pops);
    // copy order added tracker
    region_added_order.copy(plan_to_copy.region_added_order);
    // if forest graph bigger than 1 copy that 
    if(plan_to_copy.forest_graph.size() > 0){
        for (auto i = 0; i < forest_graph.size(); ++i) {
            forest_graph[i].assign(
                plan_to_copy.forest_graph[i].begin(), 
                plan_to_copy.forest_graph[i].end()
            );
        }
    }

    // if linking edges exist then copy that
    if(plan_to_copy.linking_edges.size() > 0){
        linking_edges = plan_to_copy.linking_edges;
    }
    return;
}



// Compute the log number of spanning trees on a region 
double Plan::compute_log_region_spanning_trees(MapParams const &map_params,
    int const region_id) const{
    double log_st = 0;
    // comput tau for each county intersect region
    for (int county_num = 1; county_num <= map_params.num_counties; county_num++) {
        log_st += compute_log_region_and_county_spanning_tree(
            map_params.g, map_params.counties, county_num,
            region_ids, region_id
        );
    }
    // Add county level multigraph tau
    log_st += compute_log_county_level_spanning_tree(
        map_params.g, map_params.counties, map_params.num_counties,
        region_ids,
        region_id
    );

    return log_st;
}



double Plan::compute_log_plan_spanning_trees(MapParams const &map_params) const{
    double log_st = 0;
    // compute tau from each region
    for (size_t region_id = 0; region_id < num_regions; region_id++)
    {
        log_st += compute_log_region_spanning_trees(
            map_params, region_id
        );
    }
    return log_st;
}


// Compute the log number of spanning trees on a merged region 
double Plan::compute_log_merged_region_spanning_trees(MapParams const &map_params,
    int const region1_id, int const region2_id) const{
    double log_st = 0;
    // comput tau for each county intersect region
    for (int county_num = 1; county_num <= map_params.num_counties; county_num++) {
        log_st += compute_log_region_and_county_spanning_tree(
            map_params.g, map_params.counties, county_num,
            region_ids, region1_id, region2_id
        );
    }
    // Add county level multigraph tau
    log_st += compute_log_county_level_spanning_tree(
        map_params.g, map_params.counties, map_params.num_counties,
        region_ids, region1_id, region2_id
    );

    return log_st;
}


double Plan::compute_log_linking_edge_count(
    PlanMultigraph &plan_multigraph
) const{
    return compute_log_region_multigraph_spanning_tree(
        plan_multigraph.pair_map.get_multigraph_counts(num_regions)
    );
};

// counts global number of county splits 
// Definition accoridng to McCartarn & Imai
int Plan::count_county_splits(MapParams const &map_params, std::vector<bool> &visited) const{
    // if 1 county then its just num_regions - 1
    if(map_params.num_counties == 1) return num_regions - 1;
    // 
    int num_splits = 0;
    int num_connected_components = 0; // counts the total number of connected components in all county intersect region
    // We iterate over all vertices. For each vertex we explore all
    // connected neighbors with the same district and county 

    // first reset visited
    std::fill(visited.begin(), visited.end(), false);

    for (int u = 0; u < map_params.V; u++)
    {
        // skip if we've already visited 
        if(visited[u]) continue;
        // else we've encountered a new component so increase the count
        ++num_connected_components;
        // Now we traverse all neighbors with the same county and region 
        int current_county = map_params.counties(u);
        int current_region = region_ids[u];
        std::queue<int> vertex_queue;
        vertex_queue.push(u);

        while(!vertex_queue.empty()){
            // get from queue
            int v = vertex_queue.front(); vertex_queue.pop();
            int v_region = region_ids[v];
            int v_county = map_params.counties(v);

            if(current_county != v_county || current_region != v_region){
                REprintf("BIG TIME ERROR IN COUNT SPLITS!!\n");
                throw Rcpp::exception("WOAHHHHH\n");
            }
            // mark this as visited 
            visited[v] = true;
            // go through children 
            for(auto const child_vertex: map_params.g[v]){
                // add the children if same region and county and not visited before 
                if(map_params.counties(child_vertex) == current_county && 
                   region_ids[child_vertex] == current_region &&
                   !visited[child_vertex]){
                    // mark as visited to avoid being added later  
                    visited[child_vertex] = true;
                    vertex_queue.push(child_vertex);
                }                
            }
        }

    }

    // REprintf("%d components and %d regions = %d \n", num_connected_components, map_params.num_counties,
    //     num_connected_components - map_params.num_counties);

    return num_connected_components - map_params.num_counties;

}


int Plan::count_merged_county_splits(MapParams const &map_params, std::vector<bool> &visited,
    int const region1_id, int const region2_id) const{
        // if 1 county then its just num_regions - 2
        if(map_params.num_counties == 1) return num_regions - 2;
        // 
        int num_splits = 0;
        int num_connected_components = 0; // counts the total number of connected components in all county intersect region
        // We iterate over all vertices. For each vertex we explore all
        // connected neighbors with the same district and county 
    
        // first reset visited
        std::fill(visited.begin(), visited.end(), false);
    
        for (int u = 0; u < map_params.V; u++)
        {
            // skip if we've already visited 
            if(visited[u]) continue;
            // else we've encountered a new component so increase the count
            ++num_connected_components;
            // Now we traverse all neighbors with the same county and region 
            int current_county = map_params.counties(u);
            int current_region = region_ids[u];
            // If its region2 we pretend its region 1
            if(current_region == region2_id) current_region = region1_id;
            std::queue<int> vertex_queue;
            vertex_queue.push(u);
    
            while(!vertex_queue.empty()){
                // get from queue
                int v = vertex_queue.front(); vertex_queue.pop();
                int v_region = region_ids[v];
                if(v_region == region2_id) v_region = region1_id;
                int v_county = map_params.counties(v);
    
                if(current_county != v_county || current_region != v_region){
                    REprintf("BIG TIME ERROR IN COUNT SPLITS!!\n");
                    throw Rcpp::exception("Error in count merge county splits\n");
                }
                // mark this as visited 
                visited[v] = true;
                // go through children 
                for(auto const child_vertex: map_params.county_restricted_graph[v]){
                    int child_region = region_ids[child_vertex];
                    // if its in the merged region we make it 1 by default 
                    if(child_region == region2_id) child_region = region1_id;
                    // add the children if same region and county and not visited before 
                    if(map_params.counties[child_vertex] == current_county && 
                       child_region == current_region &&
                       !visited[child_vertex]){
                        // mark as visited to avoid being added later  
                        visited[child_vertex] = true;
                        vertex_queue.push(child_vertex);
                    }                
                }
            }
    
        }
    
        // REprintf("%d components and %d regions = %d \n", num_connected_components, map_params.num_counties,
        //     num_connected_components - map_params.num_counties);
    
        return num_connected_components - map_params.num_counties;
    
    }

//' Selects a valid multidistrict to split uniformly at random 
//'
//' Given a plan object with at least one multidistrict this function randomly
//' selects a valid multidistrict to split with uniform probability (ie it is
//' one over the number of valid multidistricts to split in the plan.) and 
//' returns the log of the probability that region was chosen.
//'
//'
//' @param region_to_split an integer that will be updated by reference with the
//' id number of the region selected to split
//' @param valid_region_sizes_to_split A 1-indexed vector mapping region sizes to
//' whether or not they can be split. So `valid_region_sizes_to_split[r] == true`
//' means that multidistricts of size r can be split.  
//'
//' @details `region_id_to_split` is set to the id of the region selected
//'
//' @return the region id that was chosen to be split 
//'
int Plan::choose_multidistrict_to_split(
        std::vector<bool> const &valid_region_sizes_to_split, 
        RNGState &rng_state,
        double const selection_alpha
    ) const{

    // make vectors with cumulative d value and region label for later
    std::vector<int> valid_region_ids, associated_region_sizes;

    for(int region_id = 0 ; region_id < num_regions; region_id++) {
        auto region_size = region_sizes[region_id];
        // if valid then add id to vector 
        if(valid_region_sizes_to_split[region_size]){
            // add the count and label to vector
            valid_region_ids.push_back(region_id);
            associated_region_sizes.push_back(region_size);
        }
    }
    auto num_candidates = valid_region_ids.size();

    // If one just return that 
    if(num_candidates == 1) return valid_region_ids[0];

    arma::vec region_wgts(valid_region_ids.size());

    for (size_t i = 0; i < valid_region_ids.size(); i++)
    {
        region_wgts(i) = std::pow(associated_region_sizes[i], selection_alpha);
    }
    int idx = rng_state.r_int_unnormalized_wgt(region_wgts); 
    int region_id_to_split = valid_region_ids.at(idx);

    return region_id_to_split;
}




//' Attempts to draw a spanning tree on a region using Wilson's algorithm
//'
//' Attempts to draw a spanning tree on a specific region of a plan using 
//' Wilson's algorithm. The function will try `attempts_to_make` times to 
//' draw a tree and if all those fail then it will return false. If its 
//' successful true will be returned and the tree will be copied into the spanning forest.
//'
//' @title Attempt to Find a Valid Spanning Tree Edge to Cut into Two New Regions 
//'
//' @param map_params Map parameters (adj graph, population, etc.)
//' @param ust_sampler Uniform tree sampler object 
//' @param region_to_draw_tree_on The id of the region to draw the tree on
//' @param rng_state thread safe random number generation object
//' @param attempts_to_make How many attempts at calling wilson's algorithm to make 
//' before giving up.
//'
//' @details Modifications
//'    - `ust_sampler` is modified in place to store the new tree is sampling is 
//'         successful
//'
//' @return `true` if tree was successfully draw, `false` otherwise.
//'
//' @keyword internal
//' @noRd
std::pair<bool, int> Plan::draw_tree_on_region(
    const MapParams &map_params, const int region_to_draw_tree_on,
    Tree &ust, std::vector<bool> &visited, std::vector<bool> &ignore, int &root,
    RNGState &rng_state, int const attempts_to_make) {

    int num_attempts = 0;
    bool tree_drawn = false;
    for (size_t attempt_num = 0; attempt_num < attempts_to_make; attempt_num++)
    {
        ++num_attempts;
        // attempt to draw a tree 
        // Mark it as ignore if its not in the region to split
        for (int i = 0; i < map_params.V; i++){
            ignore[i] = region_ids[i] != region_to_draw_tree_on;
        }

        // Get a uniform spanning tree drawn on that region
        clear_tree(ust);
        // Get a tree
        int result = sample_sub_ust(map_params.g, ust, 
            map_params.V, root, visited, ignore, 
            map_params.pop, map_params.lower, map_params.upper, 
            map_params.counties, map_params.cg,
            rng_state);

        tree_drawn = result == 0;

        // if successful return 
        if(tree_drawn) break;
    }

    // return immediately if false 
    if(!tree_drawn) return std::make_pair(false, num_attempts);

    // update root region id
    // and its forest vertices
    int n_desc = ust[root].size();
    // clear this vertices neighbors in the graph and reserve size for children
    forest_graph[root].clear(); forest_graph[root].reserve(n_desc);

    // make a queue of vertex, parent 
    std::queue<std::pair<int,int>> vertex_queue;
    // add roots children to queue 
    for(auto const &child_vertex: ust[root]){
        vertex_queue.push({child_vertex, root});
        forest_graph[root].push_back(child_vertex);
    }
    
    // update all the children
    while(!vertex_queue.empty()){
        // get and remove head of queue 
        auto queue_pair = vertex_queue.front();
        int vertex = queue_pair.first;
        int parent_vertex = queue_pair.second;
        vertex_queue.pop();

        // clear this vertices neighbors in the graph and reserve size for children and parent 
        forest_graph[vertex].clear(); 
        forest_graph[vertex].reserve(ust[vertex].size()+1);
        // add the edge from vertex to parent 
        forest_graph[vertex].push_back(parent_vertex);
        
        for(auto const &child_vertex: ust[vertex]){
            // add children to queue
            vertex_queue.push({child_vertex, vertex});
            // add this edge from vertex to its children 
            forest_graph[vertex].push_back(child_vertex);
        }
    }

    return std::make_pair(true, num_attempts);

}





// Updates the region level info
//' Creates new regions and updates the `Plan` object using a cut tree
//'
//' Takes a cut spanning tree `ust` and variables on the two new regions
//' induced by the cuts and creates space/updates the information on those
//' two new regions in the `plan` object. This function increases the number
//' of regions aspect by 1 and updates the region level information and all
//' other variables changed by adding a new region. 
//'
//' It also sets `plan.remainder_region` equal to `new_region2_id` if 
//' split_district_only is true. 
//'
//'
//' @title Create and update new plan regions from cut tree
//'
//' @param ust A cut (ie has two partition pieces) directed spanning tree
//' passed by reference
//' @param plan A plan object
//' @param split_district_only Whether or not this was split according to a 
//' one district split scheme (as in does the remainder need to be updated)
//' @param old_split_region_id The id of the region that was split into the two
//' new ones. Region1 will be set to this id
//' @param new_region_id The id that region2 will be set 
//' @param new_region1_tree_root The vertex of the root of one piece of the cut
//' tree. This always corresponds to the region with the smaller dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region1_dval The dval associated with the new region 1
//' @param new_region1_pop The population associated with the new region 1
//' @param new_region2_tree_root The vertex of the root of other piece of the cut
//' tree. This always corresponds to the region with the bigger dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region2_dval The dval associated with the new region 2
//' @param new_region2_pop The population associated with the new region 2
//' @param new_region1_id The id the new region 1 was assigned in the plan
//' @param new_region2_id The id the new region 2 was assigned in the plan
//'
//' @details Modifications
//'    - `plan` is updated in place with the two new regions
//'    - `new_region1_id` is set to the id new region1 was assigned
//'    which is just the `old_split_region_id`
//'    - `new_region2_id` is set to the id new region2 was assigned 
//'    which is just `plan.num_regions-1`
//'
void Plan::update_region_info_from_cut(
        EdgeCut cut_edge, 
        const int split_region1_id, const int split_region2_id,
        bool const add_region
){
    // Get information on the two new regions cut
    int split_region1_tree_root, split_region2_tree_root;
    int split_region1_size, split_region2_size;
    int split_region1_pop, split_region2_pop;

    cut_edge.get_split_regions_info(
        split_region1_tree_root, split_region1_size, split_region1_pop,
        split_region2_tree_root, split_region2_size, split_region2_pop
    );

    // update plan with new regions if needed
    if(add_region) num_regions++; // increase region count by 1

    // make the first new region get max plus one
    region_order_max++;
    int new_region1_order_added_num = region_order_max;
    // Now make the second new region max plus one again
    region_order_max++;
    int new_region2_order_added_num = region_order_max;


    // Now update the region level information
    // updates the new region 1
    region_sizes[split_region1_id] = split_region1_size;
    region_added_order[split_region1_id] = new_region1_order_added_num;
    region_pops[split_region1_id] = split_region1_pop;

    // updates the new region 2
    // New region 2's id is the highest id number so push back
    region_sizes[split_region2_id] = split_region2_size;
    region_added_order[split_region2_id] = new_region2_order_added_num;
    region_pops[split_region2_id] = split_region2_pop;


}



void Plan::update_from_successful_split(
    TreeSplitter const &tree_splitter,
    USTSampler &ust_sampler, EdgeCut const &cut_edge,
    int const new_region1_id, int const new_region2_id,
    bool const add_region
){
    // now update the region level information from the edge cut
    update_region_info_from_cut(
        cut_edge,
        new_region1_id, new_region2_id,
        add_region
    );
    
    // Now update the vertex level information
    update_vertex_and_plan_specific_info_from_cut(
        tree_splitter,
        ust_sampler, cut_edge, 
        new_region1_id, new_region2_id,
        add_region
    ); 
     
}

// Attempts Gets valid pairs to merge for MCMC
// Returns false if multigraph was not successfully built 
std::pair<bool, std::vector<std::pair<RegionID,RegionID>>> Plan::attempt_to_get_valid_mergesplit_pairs(
        PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule
) const{
    bool const result = plan_multigraph.build_plan_multigraph(*this);
    // return false if not successful
    if(!result) return std::make_pair(false, std::vector<std::pair<RegionID,RegionID>>{}); 

    // else remove all the invalid size and mergesplit pairs 
    plan_multigraph.remove_invalid_size_pairs(*this, splitting_schedule);
    plan_multigraph.remove_invalid_mergesplit_pairs(*this);

    return std::make_pair(true, plan_multigraph.pair_map.hashed_pairs);

}

std::vector<std::pair<RegionID,RegionID>> Plan::get_valid_smc_merge_regions(
        PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule
) const{
    // build the multigraph 
    plan_multigraph.build_plan_multigraph(*this);

    std::vector<std::pair<RegionID,RegionID>> output_pairs;
    output_pairs.reserve(plan_multigraph.pair_map.num_hier_smc_merge_valid_pairs);

    // go through and only add pairs with valid size pairs 
    for(auto const a_pair: plan_multigraph.pair_map.hashed_pairs){
        // check both region sizes are valid split sizes and they are ok to merge
        bool invalid_sizing = (
            !splitting_schedule.valid_split_region_sizes[region_sizes[a_pair.first]] ||
            !splitting_schedule.valid_split_region_sizes[region_sizes[a_pair.second]] ||
            !splitting_schedule.valid_merge_pair_sizes[region_sizes[a_pair.first]][region_sizes[a_pair.second]]
        );

        if(invalid_sizing) continue;

        // now check hierarchical merge is ok 
        auto result =  plan_multigraph.pair_map.get_value(a_pair.first, a_pair.second);
        // if yes then add to output
        if(result.second.merge_is_hier_valid){
            // if yes then add to output 
            output_pairs.emplace_back(a_pair);
        }
    }
    
    return output_pairs;
}


// Prints relevant info - for debugging
void RegionPairHash::Rprint() const{
    REprintf("Pair Map has %d Elements:\n", num_hashed_pairs);
    for(auto const &a_pair: hashed_pairs){
        auto val = get_value(a_pair.first, a_pair.second);
        REprintf("    Regions: (%u, %u) | %s Hierarchically Adjacent | %d within county edges, %d across county edges\n",
            a_pair.first, a_pair.second, (val.second.admin_adjacent ? "YES" : "NOT"),
            val.second.within_county_edges, val.second.across_county_edges 
        );
    }
};

// Returns multigraph counts for a pair map that 
// has been built on a plan
RegionMultigraphCount RegionPairHash::get_multigraph_counts(int const num_regions)
const{
    RegionMultigraphCount region_multigraph(num_regions);
    // iterate over all pairs and add the proper counts 
    for(auto const key_val_pair: get_all_values()){
        // skip if not hierarchically valid merge 
        if(!key_val_pair.second.merge_is_hier_valid) continue;
        int boundary_len;
        if(key_val_pair.second.admin_adjacent){
            // if administratively adjacent only count within same county
            boundary_len = key_val_pair.second.within_county_edges;
        }else{
            // else only count across county 
            boundary_len = key_val_pair.second.across_county_edges;
        }

        // we increase the count of edges
        region_multigraph[key_val_pair.first.first][key_val_pair.first.second] = boundary_len;
        region_multigraph[key_val_pair.first.second][key_val_pair.first.first] = boundary_len;
    }

    return region_multigraph;
}


// Need to take care because some previously inelgible hier merge pairs now 
// become ok to merge 
RegionMultigraphCount RegionPairHash::get_merged_multigraph_counts(
            int const num_regions, std::vector<RegionID> &merge_index_reshuffle,
            RegionID const region1_id, RegionID const region2_id 
) const{

    RegionMultigraphCount merged_region_multigraph(num_regions-1);

    int merged_reindex = num_regions-2;
    for (int current_reindex = 0, i = 0; i < num_regions; i++){
        if(i == region1_id || i == region2_id){
            merge_index_reshuffle[i] = merged_reindex;
        }else{
            merge_index_reshuffle[i] = current_reindex;
            ++current_reindex;
        }
    }

    return(merged_region_multigraph);


    // for(auto const key_val_pair: get_all_values()){
    //     // get the regions in the pair 
    //     auto reindexed_region1_id = merge_index_reshuffle[key_val_pair.first.first];
    //     auto reindexed_region2_id = merge_index_reshuffle[key_val_pair.first.second];
    //     // Skip if the same region since that means its the merged pair itself 
    //     if(reindexed_region1_id == reindexed_region2_id) continue;
    //     // swap if out of order 
    //     if(reindexed_region2_id < reindexed_region1_id) std::swap(reindexed_region1_id, reindexed_region2_id);
    //     //
    //     // skip if not hierarchically valid merge 
    //     if(!key_val_pair.second.merge_is_hier_valid) continue;
    //     int boundary_len;
    //     if(key_val_pair.second.admin_adjacent){
    //         // if administratively adjacent only count within same county
    //         boundary_len = key_val_pair.second.within_county_edges;
    //     }else{
    //         // else only count across county 
    //         boundary_len = key_val_pair.second.across_county_edges;
    //     }

    //     // we increase the count of edges
    //     merged_region_multigraph[key_val_pair.first.first][key_val_pair.first.second] = boundary_len;
    //     merged_region_multigraph[key_val_pair.first.second][key_val_pair.first.first] = boundary_len;
    // }
}

void swap_pair_maps(RegionPairHash &a, RegionPairHash &b) {
    std::swap(a.num_hashed_pairs, b.num_hashed_pairs);
    std::swap(a.values, b.values);
    std::swap(a.hashed, b.hashed);
    std::swap(a.hashed_pairs, b.hashed_pairs);
};



PlanMultigraph::PlanMultigraph(MapParams const &map_params):
    map_params(map_params),
    counties_on(map_params.num_counties > 1),
    vertices_visited(map_params.V),
    county_component(counties_on ? map_params.ndists : 0, 0), 
    component_split_counts(counties_on ? map_params.ndists : 0, 0),
    component_region_counts(counties_on ? map_params.ndists : 0, 0),
    region_overlap_counties(counties_on ? map_params.ndists : 0),
    num_county_region_components(0),
    pair_map(map_params.ndists){

};


void PlanMultigraph::Rprint() const{
    pair_map.Rprint();
}


void PlanMultigraph::Rprint_detailed(Plan const &plan){
    REprintf("Pair Map has %d Elements:\n", pair_map.num_hashed_pairs);
    for(auto const &a_pair: pair_map.hashed_pairs){
        auto val = pair_map.get_value(a_pair.first, a_pair.second);
        REprintf("    Regions: (%u, %u) | Sizes (%u, %u) | %s Hierarchically Adjacent | %d within county edges, %d across county edges\n",
            a_pair.first, a_pair.second, 
            plan.region_sizes[a_pair.first], plan.region_sizes[a_pair.second],
            (val.second.admin_adjacent ? "YES" : "NOT"),
            val.second.within_county_edges, val.second.across_county_edges 
        );
    }
}




// indexing (i,j) in a matrix 
// assuming the length of a row is row_length
// i is the row index
// j is the column index so j < row_length
inline size_t get_county_component_lookup_index(
    int const county_id, RegionID const region_id, int const num_regions
){
    return county_id * num_regions + region_id;
}

// checks if a plan is hierarchically connected and if it is then it 
// returns the number of county region components
std::pair<bool, int> PlanMultigraph::is_hierarchically_connected(
    Plan const &plan, std::vector<bool> component_lookup
){
    // if no counties then its just one giant component 
    if(!counties_on) return std::make_pair(true, 1);
    // reset the vertices visited vector 
    std::fill(vertices_visited.begin(), vertices_visited.end(), false);
    // component_lookup should be length num_countes * num_regions
    // fill this vector with all false
    std::fill(component_lookup.begin(), component_lookup.end(), false);
    // reset the number of county region components
    num_county_region_components = 0;

    // start the counter at zero
    int component_counter = 0;

    // loop through the graph 
    for (auto v = 0; v < map_params.V; v++)
    {
        // skip if we've already visited this vertex 
        if(vertices_visited[v]) continue;
        
        // get current county 
        int current_county = map_params.counties[v]-1; 
        auto current_region = plan.region_ids[v];
        // for lookup we pretend table is num_counties x num_regions
        auto current_component_lookup_index = get_county_component_lookup_index(current_county, current_region, plan.num_regions);

        // check if we've already declared this component 
        // if default value then we haven't encountered this component yet 
        if(!component_lookup[current_component_lookup_index]){
            // increase number of connected components
            ++num_county_region_components;
        }else{
            // else we have seen this before in which case its not hierarchically connected so immediately return false
            return std::make_pair(false, -1);
        }

        // Now we traverse all neighbors in the same component 
        std::queue<int> vertex_queue;
        vertex_queue.push(v);

        while(!vertex_queue.empty()){
            int u = vertex_queue.front(); vertex_queue.pop();
            int u_region = plan.region_ids[u];
            int u_county = map_params.counties[u];

            // mark this as visited 
            vertices_visited[u] = true;
            // now visit each of the vertices children 
            for(auto const child_vertex: map_params.g[u]){
                // ignore if we already visited 
                if(vertices_visited[child_vertex]) continue;
                // Now only add if region and county are the same 
                if(u_region = plan.region_ids[child_vertex] && 
                   u_county == map_params.counties[child_vertex]){
                    // if same then same component mark as visited to avoid being added later  
                    vertices_visited[child_vertex] = true;
                    vertex_queue.push(child_vertex);
                }
            }
        }
    }
    
    // if we've passed through the whole graph then now return true 
    return std::make_pair(true, num_county_region_components);

}


void PlanMultigraph::build_plan_non_hierarchical_multigraph(
    Plan const &plan
){
    // reset the hash map 
    pair_map.reset();

    for (int v = 0; v < map_params.V; v++)
    {
        auto v_region = plan.region_ids[v];

        for (auto const u: map_params.g[v])
        {
            auto u_region = plan.region_ids[u];
            if(v_region < u_region){
                // if true then count 
                  pair_map.count_graph_edge(
                        v_region, u_region, true, 0
                    );
            }
        }
    }

    pair_map.num_hier_smc_merge_valid_pairs = pair_map.num_hashed_pairs;
    
    return;
}


bool PlanMultigraph::build_plan_hierarchical_multigraph(
    Plan const &plan
){
    // reset the hash map 
    pair_map.reset();
    // reset the vertices visited vector 
    std::fill(vertices_visited.begin(), vertices_visited.end(), false);
    // reset the region component assignments and countes 
    std::fill(county_component.begin(), county_component.end(), 0);
    std::fill(component_split_counts.begin(), component_split_counts.end(), 0);
    std::fill(component_region_counts.begin(), component_region_counts.end(), 0);
    

    // reset the number of county region components
    num_county_region_components = 0;
    // counter for labelling components of administratively connected
    // quotient graph
    int county_connected_component_counter = 0;

    // Now we walk through the graph 
    for (int w = 0; w < map_params.V; w++)
    {
        // skip if already visited 
        if(vertices_visited[w]) continue;

        // else we've encountered a new component 
        int num_current_county_region_components = 1;
        int num_current_counties = 1;

        std::queue<int> other_counties_vertices;
        std::queue<int> current_county_diff_region_vertices;
        std::queue<int> current_county_region_vertices;
        // add this vertex and mark as visited 
        current_county_region_vertices.push(w);
        // assign this region to the component
        county_component[plan.region_ids[w]] = county_connected_component_counter;

        auto current_county = map_params.counties[w];
        auto current_region = plan.region_ids[w];

        while(
            !current_county_region_vertices.empty() ||
            !current_county_diff_region_vertices.empty() ||
            !other_counties_vertices.empty()
        ){
            
            int v;

            // first see if anything in this region and county 
            if (!current_county_region_vertices.empty()) {
                v = current_county_region_vertices.front(); 
                current_county_region_vertices.pop();
            } else if (!current_county_diff_region_vertices.empty()) {
                // else we've got a vertex in a different region in the same county
                v = current_county_diff_region_vertices.front(); 
                current_county_diff_region_vertices.pop();
                // if we haven't visited this vertex yet that means we're reaching 
                // this component for the first time 
                if(!vertices_visited[v]){
                    // then change the current region and increase the count 
                    current_region = plan.region_ids[v];
                    num_current_county_region_components++;
                    // increase global count 
                    num_county_region_components++;
                    // mark the region as being in current component
                    county_component[plan.region_ids[v]] = county_connected_component_counter;
                    // REprintf("Assigned Region %u to component %d\n",
                    //     plan.region_ids[v]+1, 
                    // county_connected_component_counter);
                }
            } else if (!other_counties_vertices.empty()) {
                // else we've got a vertex in a different county
                v = other_counties_vertices.front(); 
                other_counties_vertices.pop();
                // if we haven't visited this vertex yet that means we're reaching 
                // this county for the first time 
                if(!vertices_visited[v]){
                    // then change the current region and increase the count 
                    current_county = map_params.counties[v];
                    num_current_county_region_components++;
                    num_current_counties++;
                    // increase global count 
                    num_county_region_components++;
                    // mark the region as being in current component
                    county_component[plan.region_ids[v]] = county_connected_component_counter;
                    // REprintf("Assigned Region %u to component %d\n",
                    //     plan.region_ids[v]+1, 
                    // county_connected_component_counter);
                }
            }

            // skip if we already visited 
            if (vertices_visited[v]) continue;
            // mark as visited now 
            vertices_visited[v] = true;


            auto v_region = plan.region_ids[v];
            auto v_county = map_params.counties[v];

            for (int u : map_params.g[v]) {
                // get region and county 
                auto u_region = plan.region_ids[u];
                auto u_county = map_params.counties[u];
                bool same_county = u_county == v_county;

                // if v_region < u_region then count the edge
                if(v_region < u_region){
                    pair_map.count_graph_edge(
                        v_region, u_region, same_county, u_county - 1
                    );
                }
                // if already visited don't add to the queue again
                if(vertices_visited[u]) continue;
                // else add to queue
                if(u_county == v_county){
                    if(u_region == v_region){
                        // if same region and county add to 
                        // queue of same region and county vertices
                        current_county_region_vertices.push(u);
                    }else{
                        // else if same county different region add to that 
                        // queue 
                        current_county_diff_region_vertices.push(u);
                    }
                }else if(u_region == v_region){
                    // if differerent counties then only add if same region
                    // we do not explore edges across both regions and counties 
                    other_counties_vertices.push(u);
                }
            }
        }

        // now we've visited every vertex in this county adj component so get the count 
        component_split_counts[county_connected_component_counter] = num_current_county_region_components - num_current_counties;

        if(DEBUG_BASE_PLANS_VERBOSE){
            REprintf("Component %d has %d components, %d counties so %d splits!\n",
                county_connected_component_counter, 
                num_current_county_region_components,
                num_current_counties,
                num_current_county_region_components - num_current_counties);
        }


        // now increment the counter by 1
        county_connected_component_counter++;

        // if too many global splits then auto reject
        if(num_county_region_components - map_params.num_counties >= plan.num_regions){
            return false;
        }
    }

    // now for each component add up the number of regions in the component
    for (size_t region_id = 0; region_id < plan.num_regions; region_id++)
    {
        ++component_region_counts[county_component[region_id]];
    }

    // now check for each component the number of splits is the number of regions
    // minus 1 
    for (size_t component_id = 0; component_id < county_connected_component_counter; component_id++)
    {
        if(DEBUG_BASE_PLANS_VERBOSE){
            REprintf("Component %d has %d splits and %d regions!\n",
                component_id, component_split_counts[component_id],
                component_region_counts[component_id]);
        }

        // error check, this shouldn't be possible 
        if(component_split_counts[component_id] < component_region_counts[component_id] - 1){
            REprintf("CODE ERROR!!!! Somehow county adj component has %d regions but %d splits!\n",
            component_region_counts[component_id], component_split_counts[component_id]);
        }else if(component_split_counts[component_id] >= component_region_counts[component_id]){
            return false;
        }
    }

    pair_map.num_hier_smc_merge_valid_pairs = pair_map.num_hashed_pairs;

    // Now mark which pairs are hierarchically ok to merge 
    for(auto const a_pair: pair_map.hashed_pairs){
        // if different components then its ok. Since merge ok by default then continue
        if(county_component[a_pair.first] != county_component[a_pair.second]) continue;

        // else if they're the same component check if admin adjacent
        auto hash_index = pair_map.pair_hash(a_pair.first, a_pair.second);
        // if admin adjacent do nothing 
        if(pair_map.values[hash_index].admin_adjacent) continue;
        // else reset values mark merge as false 
        pair_map.values[hash_index].merge_is_hier_valid = false;
        --pair_map.num_hier_smc_merge_valid_pairs;
    }

    return true;    
    
}

bool PlanMultigraph::build_plan_multigraph(
    Plan const &plan
){
    if(counties_on){
        return build_plan_hierarchical_multigraph(plan);
    }else{
        build_plan_non_hierarchical_multigraph(plan);
        return true;
    }
}


void PlanMultigraph::remove_invalid_size_pairs(
    Plan const &plan, SplittingSchedule const &splitting_schedule
){
    // remove invalid sizes 
    pair_map.hashed_pairs.erase(
        std::remove_if(pair_map.hashed_pairs.begin(), pair_map.hashed_pairs.end(), 
        [&](std::pair<RegionID, RegionID> a_pair) { 
            bool invalid_sizing = (
                !splitting_schedule.valid_split_region_sizes[plan.region_sizes[a_pair.first]] ||
                !splitting_schedule.valid_split_region_sizes[plan.region_sizes[a_pair.second]] ||
                !splitting_schedule.valid_merge_pair_sizes[plan.region_sizes[a_pair.first]][plan.region_sizes[a_pair.second]]
            );
            // if invalid then reset data in pair map 
            if(invalid_sizing){
                auto hash_index = pair_map.pair_hash(a_pair.first, a_pair.second);
                // else reset values and decrease hash table count 
                --pair_map.num_hashed_pairs;
                pair_map.hashed[hash_index] = false;
                pair_map.values[hash_index] = PairHashData();
            }

            // remove if any size is invalid split size or if merged size can't be split
            return invalid_sizing;
        }
    ), pair_map.hashed_pairs.end());
}


void PlanMultigraph::remove_invalid_hierarchical_merge_pairs(
    Plan const &plan
){
    if (!counties_on) return;

    // Now remove the pairs that are in the same component but not administratively adjacent 
    // remove invalid sizes 
    pair_map.hashed_pairs.erase(
        std::remove_if(pair_map.hashed_pairs.begin(), pair_map.hashed_pairs.end(), 
        [&](std::pair<RegionID, RegionID> a_pair) { 
            auto hash_index = pair_map.pair_hash(a_pair.first, a_pair.second);
            if(pair_map.values[hash_index].merge_is_hier_valid) return false;
            
            // else reset values and decrease hash table count 
            --pair_map.num_hashed_pairs;
            pair_map.hashed[hash_index] = false;
            pair_map.values[hash_index] = PairHashData();
            return true;

        }
    ), pair_map.hashed_pairs.end());

    return;

}



void PlanMultigraph::remove_invalid_mergesplit_pairs(
    Plan const &plan
){
    if (!counties_on) return;
    // reset the sets
    for (size_t i = 0; i < region_overlap_counties.size(); i++)
    {
        region_overlap_counties[i].clear();
    }
    // iterate through the pairs and note the counties they overlap in
    for(auto const &key_val_pair: pair_map.get_all_values()){
        if(key_val_pair.second.admin_adjacent){
            // if they overlap then add that county to both of the regions
            region_overlap_counties[key_val_pair.first.first].insert(key_val_pair.second.shared_county);
            region_overlap_counties[key_val_pair.first.second].insert(key_val_pair.second.shared_county);
        }
    }

    // Now we go through and delete and adjacent regions where
    // - they ARE NOT admin adjacent and they are both contain the same county 

    pair_map.hashed_pairs.erase(
        std::remove_if(pair_map.hashed_pairs.begin(), pair_map.hashed_pairs.end(), 
        [&](std::pair<RegionID, RegionID> a_pair) { 

            // if different components then its ok 
            if(county_component[a_pair.first] != county_component[a_pair.second]) return false;
            // else if they're the same component check if admin adjacent
            auto hash_index = pair_map.pair_hash(a_pair.first, a_pair.second);
            // if admin adjacent do nothing because thats ok
            if(pair_map.values[hash_index].admin_adjacent) return false;

            bool invalid_merge = false;
            // now check if they overlap in any of the two counties 
            if(region_overlap_counties[a_pair.first].size() < region_overlap_counties[a_pair.second].size()){
                // check if any element of the smaller set is in the larger set 
                for(auto const a_county: region_overlap_counties[a_pair.first]){
                    if(region_overlap_counties[a_pair.second].find(a_county) != region_overlap_counties[a_pair.second].end()){
                        invalid_merge = true;
                        break;
                    }
                }
            }else{
                for(auto const a_county: region_overlap_counties[a_pair.second]){
                    if(region_overlap_counties[a_pair.first].find(a_county) != region_overlap_counties[a_pair.first].end()){
                        invalid_merge = true;
                        break;
                    }
                }
            }
            if(invalid_merge){
                // if invalid merge deteced then remove
                --pair_map.num_hashed_pairs;
                pair_map.hashed[hash_index] = false;
                pair_map.values[hash_index] = PairHashData();
            }
            return invalid_merge;

        }
    ), pair_map.hashed_pairs.end());

    return;
    
}



bool PlanMultigraph::is_hierarchically_valid(
    Plan const &plan, std::vector<bool> component_lookup
){
    // if no counties then always valid 
    if(!counties_on) return true;
    // first check if hierarchically connected
    auto result = is_hierarchically_connected(plan, component_lookup);
    // if not immediately return false
    if(!result.first) return false;
    // if number of splits is greater than number of regions minus 1 reject
    if(result.second - map_params.num_counties >= plan.num_regions) return false;

    // Now we know its hierarchically connected and has at most num_regions-1
    // splits so we just need to check for cycles in the administratively 
    // adjacent quotient graph 
    return build_plan_hierarchical_multigraph(plan);
}


void swap_plan_multigraphs(PlanMultigraph &a, PlanMultigraph &b) {
    std::swap(a.county_component, b.county_component);
    std::swap(a.component_split_counts, b.component_split_counts);
    std::swap(a.component_region_counts, b.component_region_counts);
    std::swap(a.region_overlap_counties, b.region_overlap_counties);
    std::swap(a.num_county_region_components, b.num_county_region_components);
    swap_pair_maps(a.pair_map, b.pair_map);
}