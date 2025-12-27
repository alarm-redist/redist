/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/01
 * Purpose: Implements base Plan (ie graph partition) class 
 ********************************************************/

#include "base_plan_type.h"

bool constexpr DEBUG_BASE_PLANS_VERBOSE = false;
bool constexpr DEBUG_LOG_LINK_EDGE_VERBOSE = false;


bool Plan::check_region_pop_valid(MapParams const &map_params, int const region_id) const{
    auto region_pop = region_pops[region_id];
    auto region_size = region_sizes[region_id];
    auto region_pop_lb = map_params.lower * region_size;
    auto region_pop_ub = map_params.upper * region_size;

    return region_pop_lb <= region_pop && region_pop <= region_pop_ub;
}

std::pair<bool, std::vector<int>> Plan::all_region_pops_valid(MapParams const &map_params) const{

    bool all_valid = true;
    std::vector<int> bad_regions;

    for (int region_id = 0; region_id < num_regions; region_id++)
    {
        if(!check_region_pop_valid(map_params, region_id)){
            all_valid = false;
            bad_regions.push_back(region_id);
        }        
    }

    return std::make_pair(all_valid, bad_regions);
}


std::pair<bool, std::vector<int>> Plan::all_regions_connected(
        Graph const &g, CircularQueue<int> &vertex_queue,
        std::vector<bool> &vertices_visited,
        std::vector<bool> &regions_visited
){
    std::fill(regions_visited.begin(), regions_visited.end(), false);
    std::fill(vertices_visited.begin(), vertices_visited.end(), false);
    vertex_queue.clear();
    int const V = g.size();

    std::set<int> disconnected_regions;

    for (int v = 0; v < V; v++)
    {
        // skip if we've already visited this
        if(vertices_visited[v]){
            continue;
        }
        auto v_region = region_ids[v];
        // If we've already visited this region but not a vertex
        // then the region is disconnected
        if(regions_visited[v_region] && !vertices_visited[v]){
            disconnected_regions.insert(v_region);
        }
        // mark this vertex as visited
        vertices_visited[v] = true;
        // mark this region as visited
        regions_visited[v_region] = true;
        vertex_queue.push(v);
        // now visit all its neighbors
        while (!vertex_queue.empty())
        {
            auto u = vertex_queue.pop();
            auto u_region = region_ids[u];
            // mark this as visited 
            vertices_visited[u] = true;

            // add unvisited neighbors
            for(auto const &u_nbor : g[u]){
                if(vertices_visited[u_nbor]) continue;
                auto u_nbor_region = region_ids[u_nbor];
                if(u_nbor_region == v_region){
                    vertex_queue.push(u_nbor);
                    vertices_visited[u_nbor] = true;
                }
            }
        }
    }

    std::vector<int> disconnected_region_output(
        disconnected_regions.begin(),
        disconnected_regions.end()
    );
    bool any_disconnected = disconnected_region_output.size() > 0;

    return std::make_pair(any_disconnected, disconnected_region_output);
    
}

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
    // gather the counties to iterate over 
    std::set<int> region_counties;
    // // if one county just that 
    if(map_params.num_counties == 1){
        // if there's only 1 county 
        region_counties.insert(map_params.counties[0]);
    }else{
        // just add the counties in that region 
        for (size_t v = 0; v < map_params.V; v++)
        {
            if(region_ids[v] == region_id){
                region_counties.insert(
                    map_params.counties[v]
                );
            }
        }
    }

    // comput tau for each county intersect region
    for (auto const county_num: region_counties) {
        log_st += compute_log_region_and_county_spanning_tree_eigen_tri(
            map_params.g, map_params.counties, county_num,
            region_ids, region_id, region_id
        );
        // REprintf("Called for %d!\n", county_num);
        // log_st += compute_log_region_and_county_spanning_tree(
        //     map_params.g, map_params.counties, county_num,
        //     region_ids, region_id, region_id
        // );

        // REprintf("%f vs %f, diff %.20f and %d\n", 
        //     new_tau, old_tau, 
        //     std::fabs(new_tau - old_tau),
        //     old_tau == new_tau);
    }
    // Add county level multigraph tau
    log_st += compute_log_county_level_spanning_tree(
        map_params.g, map_params.counties, map_params.num_counties,
        region_ids,
        region_id, region_id
    );

    return log_st;
}

// Compute the log number of spanning trees on a merged region 
double Plan::compute_log_merged_region_spanning_trees(MapParams const &map_params,
    int const region1_id, int const region2_id) const{
    double log_st = 0;
    // gather the counties to iterate over 
    std::set<int> region_counties;
    // // if one county just that 
    if(map_params.num_counties == 1){
        // if there's only 1 county 
        region_counties.insert(map_params.counties[0]);
    }else{
        // just add the counties in that region 
        for (size_t v = 0; v < map_params.V; v++)
        {
            if(region_ids[v] == region1_id || region_ids[v] == region2_id){
                region_counties.insert(
                    map_params.counties[v]
                );
            }
        }
    }

    // comput tau for each county intersect region
    for (auto const county_num: region_counties) {
        log_st += compute_log_region_and_county_spanning_tree_eigen_tri(
            map_params.g, map_params.counties, county_num,
            region_ids, region1_id, region2_id
        );
        // REprintf("Called for %d!\n", county_num);
        // log_st += compute_log_region_and_county_spanning_tree(
        //     map_params.g, map_params.counties, county_num,
        //     region_ids, region1_id, region2_id
        // );
    }
    // Add county level multigraph tau
    log_st += compute_log_county_level_spanning_tree(
        map_params.g, map_params.counties, map_params.num_counties,
        region_ids, region1_id, region2_id
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

    Tree county_tree = init_tree(map_params.num_counties);
    TreePopStack county_stack(map_params.num_counties);
    arma::uvec county_pop(map_params.num_counties, arma::fill::zeros);
    std::vector<std::vector<int>> county_members(map_params.num_counties, std::vector<int>{});
    std::vector<bool> c_visited(map_params.num_counties, true);
    std::vector<int> cty_pop_below(map_params.num_counties, 0);
    std::vector<std::array<int, 3>> county_path;
    std::vector<int> path;


    int num_attempts = 0;
    bool tree_drawn = false;
    auto the_region_size = region_sizes[region_to_draw_tree_on];
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
        int result = sample_sub_ust(
            map_params, ust, root, 
            map_params.lower * the_region_size, map_params.upper * the_region_size,
            visited, ignore, county_tree, county_stack, county_pop, county_members, 
            c_visited, cty_pop_below, county_path, path,
            rng_state
        );            

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
        PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule,
        ScoringFunction const &scoring_function, bool const is_final_split
) const{
    bool const result = plan_multigraph.build_plan_multigraph(region_ids, num_regions);
    // return false if not successful
    if(!result) return std::make_pair(false, std::vector<std::pair<RegionID,RegionID>>{}); 
    // plan_multigraph.Rprint_detailed(*this);
    // else remove all the invalid size and mergesplit pairs and invalid hard constraint
    plan_multigraph.remove_invalid_size_pairs(*this, splitting_schedule);
    plan_multigraph.remove_invalid_hard_constraint_pairs(*this, scoring_function, is_final_split);
    plan_multigraph.remove_invalid_mergesplit_pairs(*this);
    // plan_multigraph.Rprint_detailed(*this);

    return std::make_pair(true, plan_multigraph.pair_map.hashed_pairs);

}

std::vector<std::pair<RegionID,RegionID>> Plan::get_valid_smc_merge_regions(
        PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule,
        ScoringFunction const &scoring_function, bool const is_final_split
) const{
    // build the multigraph 
    plan_multigraph.build_plan_multigraph(region_ids, num_regions);
    

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
        auto result = plan_multigraph.pair_map.get_value(a_pair.first, a_pair.second);
        if(!result.second.merge_is_hier_valid) continue;

        // now check hard constraints are satisfied
        bool const hard_constr_result = scoring_function.merged_plan_ok(
            *this, a_pair.first, a_pair.second, is_final_split
        );

        if(!hard_constr_result) continue;

        // we've passed all checks so now we can add 
        output_pairs.emplace_back(a_pair);
    }
    
    return output_pairs;
}


// Prints relevant info - for debugging
void RegionPairHash::Rprint(std::vector<int> const &county_component) const{
    REprintf("Pair Map has %d Elements:\n", num_hashed_pairs);
    for(auto const &a_pair: hashed_pairs){
        auto val = get_value(a_pair.first, a_pair.second);
        REprintf("    Regions: (%u, %u) | Components (%d, %d) %s | %s Admin Adjacent | %d within county edges, %d across county edges\n",
            a_pair.first, a_pair.second,
            county_component[a_pair.first], county_component[a_pair.second],
            (val.second.same_admin_component ? "SAME" : "NOT SAME"),
             (val.second.admin_adjacent ? "YES" : "NOT"),
            val.second.within_county_edges, val.second.across_county_edges 
        );
    }
};




double PlanMultigraph::compute_non_hierarchical_log_multigraph_tau(
    int const num_regions, ScoringFunction const &scoring_function
){
    // if two regions then its just the log(boundary length) between
    // Since any minor of laplacian is just the degree of vertex 0 or 1
    if(num_regions == 2){
        // Get the value 
        auto val = pair_map.get_value(0, 1).second;
        if(val.admin_adjacent){
            return std::log(static_cast<double>(val.within_county_edges));
        }else{
            return std::log(static_cast<double>(val.across_county_edges));
        }
    }

    // else go through and build the laplacian for regions 0 through num_regions-1
    arma::mat laplacian_minor(num_regions-1, num_regions-1, arma::fill::zeros);

    // Now we iterate through the pairs 
    for (auto const a_pair: pair_map.hashed_pairs){
        // Get the value 
        auto val = pair_map.get_value(a_pair.first, a_pair.second).second;
        // Since we're assuming non hierarchical there only one 
        // component so all edges are within county edges


        int const edges = val.within_county_edges;;
        // Now for both pairs we need to 
        // 1. Update the degree of the vertex
        // 2. Increase the count of edges between the regions 

        // But if either region has index num_regions-1
        // we don't count its degree or care about edge count

        // Check if neither pair is id num_regions-1
        if(a_pair.first != num_regions-1 && a_pair.second != num_regions-1){
            // increase the degree of both vertices 
            laplacian_minor(a_pair.first, a_pair.first) += edges;
            laplacian_minor(a_pair.second, a_pair.second) += edges;
            // subtract edges between them 
            laplacian_minor(a_pair.first, a_pair.second) -= edges;
            laplacian_minor(a_pair.second, a_pair.first) -= edges;
        }else if(a_pair.first != num_regions-1){
            // increase the degree 
            laplacian_minor(a_pair.first, a_pair.first) += edges;
            // because other region is merged one ignore 
        }else if(a_pair.second != num_regions-1){
            // increase the degree 
            laplacian_minor(a_pair.second, a_pair.second) += edges;
            // because other region is merged one ignore 
        }

    }

    return arma::log_det_sympd(laplacian_minor);
};


double PlanMultigraph::compute_non_hierarchical_merged_log_multigraph_tau(
    int const num_regions, 
    RegionID const region1_id, RegionID const region2_id,
    ScoringFunction const &scoring_function
){
    // If 2 regions then merged is 1 region and there's exactly 1 way to draw
    // linking edges (since its empty set) so log(1) = 0
    if(num_regions == 2) return 0.0;

    /* 
    we make merged reindex which essentially just reindexes all the regions to 
    the two merged regions are the new largest value and everything else gets
    potentially shifted down. For example imagine 5 regions
    and we merge 1 and 3 so the reindex vector is [0,3,1,3,2]
     */
    int const merged_reindex = num_regions-2;
    for (int current_reindex = 0, i = 0; i < num_regions; i++){
        if(i == region1_id || i == region2_id){
            region_reindex_vec[i] = merged_reindex;
        }else{
            region_reindex_vec[i] = current_reindex;
            ++current_reindex;
        }
    }

    /*
    Recall the graph laplacian for a multigraph with V vertices is a VxV matric
    where
        - v_ii is the degree of the vertex (ie the total number of edges where 
        v_ii is a vertex)
        - v_ij is -m_ij where m_ij is the number of edges between vertex i and j
    
    Kirchenoff's theorem says the determinant of any minor of laplcacian is number
    of spanning trees. A minor is laplacian with any row and column deleted. 

    So for the merged plan it has num_regions-1 vertices and since we only compute 
    determinant of the minor we actually have a num_regions-2 x num_regions-2 matrix
    where we delete the row and column corresponding to the merged region. 
     */
    arma::mat merged_laplacian_minor(num_regions-2, num_regions-2, arma::fill::zeros);



    // Now we iterate through the pairs 
    for (auto const a_pair: pair_map.hashed_pairs){
        // Get the value 
        auto val = pair_map.get_value(a_pair.first, a_pair.second).second;

        auto reshuffled_pair1 = region_reindex_vec[a_pair.first];
        auto reshuffled_pair2 = region_reindex_vec[a_pair.second];

        // If its (region1, region2) then skip because merged 
        if(reshuffled_pair1 == reshuffled_pair2) continue;

        // since non hierarchcial every merge pair is valid and all edges are
        // within county edges

        int edges = val.within_county_edges;
        // Now for both pairs we need to 
        // 1. Update the degree of the vertex
        // 2. Increase the count of edges between the regions 

        // But if either region is part of the merged region then 
        // we don't count its degree or care about edge count

        // Check if neither pair is merged region
        if(reshuffled_pair1 != merged_reindex && reshuffled_pair2 != merged_reindex){
            // increase the degree of both vertices 
            merged_laplacian_minor(reshuffled_pair1, reshuffled_pair1) += edges;
            merged_laplacian_minor(reshuffled_pair2, reshuffled_pair2) += edges;
            // subtract edges between them 
            merged_laplacian_minor(reshuffled_pair1, reshuffled_pair2) -= edges;
            merged_laplacian_minor(reshuffled_pair2, reshuffled_pair1) -= edges;
        }else if(reshuffled_pair1 != merged_reindex){
            // increase the degree 
            merged_laplacian_minor(reshuffled_pair1, reshuffled_pair1) += edges;
            // because other region is merged one ignore 
        }else if(reshuffled_pair2 != merged_reindex){
            // increase the degree 
            merged_laplacian_minor(reshuffled_pair2, reshuffled_pair2) += edges;
            // because other region is merged one ignore 
        }

    }

    return arma::log_det_sympd(merged_laplacian_minor); 
}



double PlanMultigraph::compute_log_multigraph_tau(
    int const num_regions, 
    ScoringFunction const &scoring_function
){
    if(counties_on){
        return compute_hierarchical_log_multigraph_tau(
            num_regions, scoring_function
        );
    }else{
        return compute_non_hierarchical_log_multigraph_tau(
            num_regions, scoring_function
        );
    }
}


double PlanMultigraph::compute_merged_log_multigraph_tau(
            int const num_regions, 
            RegionID const region1_id, RegionID const region2_id,
            ScoringFunction const &scoring_function
){
    if(counties_on){
        return compute_hierarchical_merged_log_multigraph_tau(
            num_regions, 
            region1_id, region2_id,
            scoring_function
        );
    }else{
        return compute_non_hierarchical_merged_log_multigraph_tau(
            num_regions, 
            region1_id, region2_id,
            scoring_function
        );
    }
}

double PlanMultigraph::compute_hierarchical_log_multigraph_tau(
    int const num_regions, 
    ScoringFunction const &scoring_function
) {
    // if two regions then its just the log(boundary length) between
    // Since any minor of laplacian is just the degree of vertex 0 or 1
    if(num_regions == 2){
        // Get the value 
        auto val = pair_map.get_value(0, 1).second;
        if(val.admin_adjacent){
            return std::log(static_cast<double>(val.within_county_edges));
        }else{
            return std::log(static_cast<double>(val.across_county_edges));
        }
    }

    // We need to compute this hierarchically. First we compute within each 
    // connected component then we compute across all of them 
    // First we need to sort all the pairs

    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    Rprint();
    }

    auto all_pairs = pair_map.get_all_values(true);

    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    REprintf("%d Components | Pre-Sorted Pairs: \n", num_county_connected_components);
    for(auto const &val: all_pairs){
        REprintf("Regions (%u, %u) | Components (%d, %d) | Shared Status %s\n",
            val.first.first, val.first.second, 
            county_component[val.first.first], county_component[val.first.second],
            (val.second.same_admin_component ? "SHARED" : "NOT SHARED"));
    }
    }

    // Sort the pairs to 
    //  1. All pairs that are in different components are at the end
    //  2. Among pairs in the same component sort them by their component  
    std::sort(all_pairs.begin(), all_pairs.end(),
        [&](const auto &a, const auto &b) {
            const auto &[region_1a, region_2a] = a.first;
            const PairHashData &data_a = a.second;

            const auto &[region_1b, region_2b] = b.first;
            const PairHashData &data_b = b.second;

            // Rule 1: same_admin_component == true goes to the front
            if (data_a.same_admin_component != data_b.same_admin_component) {
                // recall true =1, false = 0 so this says if data_a and data_b 
                // different then a is greater iff its true 
                return data_a.same_admin_component > data_b.same_admin_component;
            }

            // Rule 2: for matching entries, sort by county_component[region1]
            if (data_a.same_admin_component) {
                return county_component[region_1a] < county_component[region_1b];
            }

            return false; // otherwise, keep existing relative order
        });


    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    REprintf("%d Components | NOW SORTED Pairs: \n", num_county_connected_components);
    for(auto const &val: all_pairs){
        REprintf("Regions (%u, %u) | Components (%d, %d) | Shared Status %s\n",
            val.first.first, val.first.second, 
            county_component[val.first.first], county_component[val.first.second],
            (val.second.same_admin_component ? "SHARED" : "NOT SHARED"));
    }
    }


    double log_tau = 0.0;

    // current index through all_pairs
    int curr_index = 0;

    // Now we go through each component and compute log spanning trees 
    for (int component_id = 0; component_id < num_county_connected_components; component_id++)
    {

        // stop loop if we've reached pairs across components
        if(curr_index >= all_pairs.size() || !all_pairs[curr_index].second.same_admin_component) break;

        int num_component_regions = component_region_counts[component_id];
        
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
            REprintf("Component %d which has %d regions!\n", component_id, num_component_regions);
        }

        // If 1 component then there only 1 tree so log(1) = 0
        if(num_component_regions <= 1){
            // do nothing
            // REprintf("1 Region in Component! Log tau unchanged at %f\n", log_tau);
            // do nothing 
        }else if(num_component_regions == 2){
            // If its two then its just the degree of either vertex which is the boundary length
            // Since its the same component its just within component edges
            log_tau += std::log(static_cast<double>(all_pairs[curr_index].second.within_county_edges));
            // REprintf("2 Regions in Component! Log tau now at %f\n", log_tau);
            // now increment index
            ++curr_index;
        }else{
            if(DEBUG_LOG_LINK_EDGE_VERBOSE){
            REprintf("3 or more regions in component!\n");
            }
            // else more than two regions so we need to reindex the regions 
            // This ensures each region in the component is indexed from 
            // 0, ... num regions in the component 
            for (int current_reindex = 0, a_region_id = 0; a_region_id < num_regions; a_region_id++){
                // if in the component then reindex that region
                if(county_component[a_region_id] == component_id){
                    region_reindex_vec[a_region_id] = current_reindex;
                    ++current_reindex;
                    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
                    REprintf("Mapping %d to %d!\n", a_region_id, region_reindex_vec[a_region_id]);
                    }
                }

            }

            // Now make the graph laplacian matrix 
            arma::mat laplacian_minor(num_component_regions-1, num_component_regions-1, arma::fill::zeros);

            // now we iterate through all pairs where both are in this component
            while(
                curr_index < all_pairs.size() &&
                county_component[all_pairs[curr_index].first.first] == component_id && 
                county_component[all_pairs[curr_index].first.second] == component_id 
            ){

                const auto &[pair_region1, pair_region2] = all_pairs[curr_index].first;
                const PairHashData &pair_val = all_pairs[curr_index].second;

                // Because we filter out hier merge invalid we don't need to check

                int edges = 0;

                // else valid merge so get edges we count
                if(pair_val.admin_adjacent){
                    edges = all_pairs[curr_index].second.within_county_edges;
                }else{
                    REprintf("Region Pair (%d, %d) - ", pair_region1, pair_region2);
                    REprintf("Error! in BASE PLAN TYPE NON ADMIN ADJACENT HIER MERGE VALID!\n");
                    throw Rcpp::exception("Bro!]n");
                }

                if(DEBUG_LOG_LINK_EDGE_VERBOSE){
                    REprintf("Region Pair (%d, %d) | %d | %d edges!\n ", 
                        pair_region1, pair_region2, 
                        pair_val.admin_adjacent,
                        edges);
                }

                auto reindexed_pair_region1 = region_reindex_vec[pair_region1];
                auto reindexed_pair_region2 = region_reindex_vec[pair_region2];

                // Now for both pairs we need to 
                // 1. Update the degree of the vertex
                // 2. Increase the count of edges between the regions 

                // But if either region has index num_regions-1
                // we don't count its degree or care about edge count

                // Check if neither pair is id num_regions-1
                if(reindexed_pair_region1 != num_component_regions-1 && 
                   reindexed_pair_region2 != num_component_regions-1){
                    // increase the degree of both vertices 
                    laplacian_minor(reindexed_pair_region1, reindexed_pair_region1) += edges;
                    laplacian_minor(reindexed_pair_region2, reindexed_pair_region2) += edges;
                    // subtract edges between them 
                    laplacian_minor(reindexed_pair_region1, reindexed_pair_region2) -= edges;
                    laplacian_minor(reindexed_pair_region2, reindexed_pair_region1) -= edges;
                }else if(reindexed_pair_region1 != num_component_regions-1){
                    // increase the degree 
                    laplacian_minor(reindexed_pair_region1, reindexed_pair_region1) += edges;
                    // because other region is merged one ignore 
                }else if(reindexed_pair_region2 != num_component_regions-1){
                    // increase the degree 
                    laplacian_minor(reindexed_pair_region2, reindexed_pair_region2) += edges;
                    // because other region is merged one ignore 
                }
                // increase current index by 1
                ++curr_index;
            }
            if(DEBUG_LOG_LINK_EDGE_VERBOSE){
            REprintf("Printing Laplacian Minor!\n");
            laplacian_minor.print();
            }
            // Now add log det
            log_tau += arma::log_det_sympd(laplacian_minor);
        }
    }

    if(curr_index >= all_pairs.size() && num_county_connected_components > 1){
        REprintf("ERROR!!\n %u pairs, curr index %d but %d num admin components\n",
            all_pairs.size(), curr_index, num_county_connected_components);
            Rprint();
    REprintf("%d Components | NOW SORTED Pairs: \n", num_county_connected_components);
    for(auto const &val: all_pairs){
        REprintf("Regions (%u, %u) | Components (%d, %d) | Shared Status %s\n",
            val.first.first, val.first.second, 
            county_component[val.first.first], county_component[val.first.second],
            (val.second.same_admin_component ? "SHARED" : "NOT SHARED"));
    }
    for (int i = 0; i < num_county_connected_components; i++)
    {
        REprintf("Component %d - %d Regions\n", i, component_region_counts[i]);
    }
        throw Rcpp::exception("! Hier\n");
    }

    // If only one connected component then we stop 
    if(num_county_connected_components == 1){
        return log_tau;
    }else if(num_county_connected_components == 2){
        int across_component_edges = 0;
        // If two then we just sum the edges across and take log 
        for (int i = curr_index; i < all_pairs.size(); i++)
        {
            across_component_edges += all_pairs[i].second.across_county_edges;
        }
        // add
        log_tau += std::log(static_cast<double>(across_component_edges));
        return log_tau;
    }

    // Now we compute spanning trees across components 
    arma::mat component_laplacian_minor(
        num_county_connected_components-1, 
        num_county_connected_components-1, 
        arma::fill::zeros);

    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    REprintf("Current index now %d with %u pairs\n", curr_index, all_pairs.size());
    }
    
    // Iterate over the remaining pairs that are across components 
    for (int i = curr_index; i < all_pairs.size(); i++)
    {
        const auto &[pair_region1, pair_region2] = all_pairs[i].first;
        const PairHashData &pair_val = all_pairs[i].second;

        int component1_id = county_component[pair_region1];
        int component2_id = county_component[pair_region2];

        // 
        if(pair_val.same_admin_component){
            REprintf("BIG ERROR SAME COMPONENT WHEN SHOULD BE DIFF!\n");
        }

        int edges = pair_val.across_county_edges;
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("Components (%d,%d) - %d edges across!\n",
            component1_id,component2_id, edges);
        }
        // Now for both pairs we need to 
        // 1. Update the degree of the vertex
        // 2. Increase the count of edges between the regions 

        // But if either region has index num_regions-1
        // we don't count its degree or care about edge count

        // Check if neither pair is id num_regions-1
        if(component1_id != num_county_connected_components-1 && 
            component2_id != num_county_connected_components-1){
            // increase the degree of both vertices 
            component_laplacian_minor(component1_id, component1_id) += edges;
            component_laplacian_minor(component2_id, component2_id) += edges;
            // subtract edges between them 
            component_laplacian_minor(component1_id, component2_id) -= edges;
            component_laplacian_minor(component2_id, component1_id) -= edges;
        }else if(component1_id != num_county_connected_components-1){
            // increase the degree 
            component_laplacian_minor(component1_id, component1_id) += edges;
            // because other region is merged one ignore 
        }else if(component2_id != num_county_connected_components-1){
            // increase the degree 
            component_laplacian_minor(component2_id, component2_id) += edges;
            // because other region is merged one ignore 
        }
    }
    
    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("Printing Component Laplacian Minor\n");
    component_laplacian_minor.print();
    }

    // Now add log det
    log_tau += arma::log_det_sympd(component_laplacian_minor);

    return log_tau;
}


// Need to take care because some previously inelgible hier merge pairs now 
// become ok to merge 
double PlanMultigraph::compute_hierarchical_merged_log_multigraph_tau(
            int const num_regions, 
            RegionID const region1_id, RegionID const region2_id,
            ScoringFunction const &scoring_function
){
    // if two regions then their merge is just 1 region which has no linking 
    // edges. Size of empty set is 0 so return log(1) = 0
    if(num_regions <= 2){
        return 0.0;
    }

    // Determine if the merge is joining two admin connected components
    auto region1_component = county_component[region1_id];
    auto region2_component = county_component[region2_id];
    bool const two_components_merged = region1_component != region2_component;
    int merged_num_admin_connected_components = num_county_connected_components - two_components_merged;


    // We need to compute this hierarchically. First we compute within each 
    // connected component then we compute across all of them 
    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("We are merging %u and %u! in unmereged components( %d, %d) | There are %d Unmerged Components and %d merged components | \n", 
            region1_id, region2_id,
            region1_component, region2_component,
            num_county_connected_components,
            merged_num_admin_connected_components
        );
        Rprint();
    }

    // First get all pairs ignoring invalid merges 
    auto all_pairs = pair_map.get_all_values(true);

    // We'll use the reindex vector to first reindex component ids 
    // Since there can't be more components than regions 
    // this will just map component two to component 1
    for (int a_region_id = 0; a_region_id < num_regions; a_region_id++){
        auto const this_region_component = county_component[a_region_id];
        // all components but region 2's get mapped to themselves
        if(this_region_component == region2_component){
            county_component_reindex[a_region_id] = region1_component;
        }else{
            county_component_reindex[a_region_id] = this_region_component;
        }
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("Region %d maps to component %d now! (was %d) \n", 
            a_region_id, 
            county_component_reindex[a_region_id], 
            county_component[a_region_id]
            );
        }
    }

    if(two_components_merged){
        // If we're joining two components then we need to filter out pairs
        // that are not admin adjacent but were previously in different components
        all_pairs.erase(
            std::remove_if(all_pairs.begin(), all_pairs.end(), 
            [&](const std::pair<std::pair<RegionID, RegionID>, PairHashData> &an_el) { 
                const auto &[region_1a, region_2a] = an_el.first;
                const PairHashData &data_a = an_el.second;
            
                // We remove if both are same component and not admin adjacent
                bool const invalid_merge_now = (
                    county_component_reindex[region_1a] == county_component_reindex[region_2a] &&
                    !data_a.admin_adjacent
                );

                // remove if any size is invalid split size or if merged size can't be split
                return invalid_merge_now;
            }
        ), all_pairs.end());
    }else{        
        // We're merging two regions in the same admin connected component so just 
        // remove that pair
        all_pairs.erase(
            std::remove_if(all_pairs.begin(), all_pairs.end(), 
            [&](auto an_el) { 
                const auto &[region_1a, region_2a] = an_el.first;
                return region_1a == region1_id && region_2a == region2_id;
            }
        ), all_pairs.end());
    }

    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    REprintf("%d Components | Pre-Sorted Pairs: \n", merged_num_admin_connected_components);
    for(auto const &val: all_pairs){
        REprintf("Regions (%u, %u) | Components (%d, %d) | Shared Status %s | Edges %d\n",
            val.first.first, val.first.second, 
            county_component_reindex[val.first.first], county_component_reindex[val.first.second],
            (val.second.same_admin_component ? "SHARED" : "NOT SHARED"),
            (val.second.same_admin_component ? val.second.within_county_edges : val.second.across_county_edges));
    }
    }


    // Sort the pairs to 
    //  1. All pairs that are in different components are at the end
    //  2. Among pairs in the same component sort them by their component 
        // need to account for one less component in sorting 
    std::sort(all_pairs.begin(), all_pairs.end(),
        [&](const auto &a, const auto &b) {
            const auto &[region_1a, region_2a] = a.first;

            bool const data_a_same_component = (
                county_component_reindex[region_1a] == county_component_reindex[region_2a]
            );

            const auto &[region_1b, region_2b] = b.first;


            bool const data_b_same_component = (
                county_component_reindex[region_1b] == county_component_reindex[region_2b]
            );

            // Rule 1: same_admin_component == true goes to the front
            if (data_a_same_component != data_b_same_component) {
                // recall true =1, false = 0 so this says if data_a and data_b 
                // different then a is greater iff its true 
                return data_a_same_component > data_b_same_component;
            }

            // Rule 2: for matching entries, sort by reindexed county_component[region1]
            if (data_a_same_component) {
                return county_component_reindex[region_1a] < county_component_reindex[region_1b];
            }

            return false; // otherwise, keep existing relative order
        });

    
    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    REprintf("%d Components | NOW SORTED Pairs: \n", merged_num_admin_connected_components);
    for(auto const &val: all_pairs){
        REprintf("Regions (%u, %u) | Components (%d, %d) | Shared Status %s | Edges %d\n",
            val.first.first, val.first.second, 
            county_component_reindex[val.first.first], county_component_reindex[val.first.second],
            (val.second.same_admin_component ? "SHARED" : "NOT SHARED"),
            (val.second.same_admin_component ? val.second.within_county_edges : val.second.across_county_edges));
    }
    }


    double log_tau = 0.0;

    // current index through all_pairs
    int curr_index = 0;

    // Now we go through each component and compute log spanning trees 
    for (int component_id = 0; component_id < num_county_connected_components; component_id++)
    {
        // Number of components in the region is same if not component 1
        // If component 1 then we add up the two regions 
        int num_component_regions;
        if(two_components_merged){
            if(component_id == region2_component){
                // if two components merged then we ignore the second one 
                continue;
            }else if(component_id == region1_component){
                // number of components is the sum of the regions 
                // minus 1 because merged 
                num_component_regions = component_region_counts[component_id] + component_region_counts[region2_component] - 1;
            }else{
                num_component_regions = component_region_counts[component_id];
            }
        }else{
            num_component_regions = component_region_counts[component_id];
            // If merged component subtract a region
            if(component_id == region1_component){
                --num_component_regions;
            }
        }
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
            REprintf("Component %d which has %d regions!\n", component_id, num_component_regions);
        }

        // stop loop if we've reached pairs across components
        if(curr_index >= all_pairs.size() || 
            (county_component_reindex[all_pairs[curr_index].first.first]  != 
             county_component_reindex[all_pairs[curr_index].first.second])
            ) break;
        


        // If 1 component then there only 1 tree so log(1) = 0
        if(num_component_regions <= 1){
            // do nothing
            if(DEBUG_LOG_LINK_EDGE_VERBOSE){
            REprintf("1 Region in Component! Log tau unchanged at %f\n", log_tau);
            }
            // do nothing 
        }else if(num_component_regions == 2){
            // If its two then its just the degree of either vertex which is the boundary length
            // Since its the same component its just within component edges

            // If this is the merged component then there are potentially more than 2 pairs
            if(component_id == region1_component){
                int merged_boundary = 0;
                while(
                    curr_index < all_pairs.size() &&
                    county_component_reindex[all_pairs[curr_index].first.first] == component_id && 
                    county_component_reindex[all_pairs[curr_index].first.second] == component_id 
                ){
                    merged_boundary += all_pairs[curr_index].second.within_county_edges;
                    ++curr_index;
                }
                log_tau += std::log(static_cast<double>(merged_boundary));
                if(DEBUG_LOG_LINK_EDGE_VERBOSE){
                REprintf("2 Regions in Merged Component! | Boundary length %d |Log tau now at %f\n", 
                    merged_boundary,
                    log_tau);
                }
            }else{
                // else if not the merged component there is actually only one pair 
                log_tau += std::log(static_cast<double>(all_pairs[curr_index].second.within_county_edges));
                if(DEBUG_LOG_LINK_EDGE_VERBOSE){
                REprintf("2 Regions in Non-merged Component! | Boundary length %d |Log tau now at %f\n", 
                    all_pairs[curr_index].second.within_county_edges,
                    log_tau);
                }
                // now increment index
                ++curr_index;
            }
        }else{
            if(DEBUG_LOG_LINK_EDGE_VERBOSE){
            REprintf("3 or more regions in component!\n");
            }
            if(component_id == region1_component){
                // We will reindex the merged region to the last index in the component 
                for (int current_reindex = 0, a_region_id = 0; a_region_id < num_regions; a_region_id++){
                    if(a_region_id == region1_id || a_region_id == region2_id){
                        // the regions we're merging get reindex to the last region in the component 
                        region_reindex_vec[a_region_id] = num_component_regions - 1;
                    }else if(county_component_reindex[a_region_id] == component_id){
                        region_reindex_vec[a_region_id] = current_reindex;
                        ++current_reindex;
                    }
                    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
                    REprintf("Mapping %d to %d!\n", a_region_id, region_reindex_vec[a_region_id]);
                    }
                } 
            }else{
                // else more than two regions so we need to reindex the regions 
                // This ensures each region in the component is indexed from 
                // 0, ... num regions in the component 
                for (int current_reindex = 0, a_region_id = 0; a_region_id < num_regions; a_region_id++){
                    // if in the component then reindex that region
                    if(county_component_reindex[a_region_id] == component_id){
                        region_reindex_vec[a_region_id] = current_reindex;
                        ++current_reindex;
                        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
                        REprintf("Mapping %d to %d!\n", a_region_id, region_reindex_vec[a_region_id]);
                        }
                    }
                } 
            }

            // Now make the graph laplacian matrix 
            arma::mat laplacian_minor(num_component_regions-1, num_component_regions-1, arma::fill::zeros);

            // now we iterate through all pairs where both are in this component
            while(
                curr_index < all_pairs.size() &&
                county_component_reindex[all_pairs[curr_index].first.first] == component_id && 
                county_component_reindex[all_pairs[curr_index].first.second] == component_id 
            ){

                const auto &[pair_region1, pair_region2] = all_pairs[curr_index].first;
                const PairHashData &pair_val = all_pairs[curr_index].second;

                // Because we filter out hier merge invalid we don't need to check

                int edges = 0;
                // else valid merge so get edges we count
                if(pair_val.admin_adjacent){
                    edges = all_pairs[curr_index].second.within_county_edges;
                }else{
                    REprintf("Region Pair (%d, %d) - ", pair_region1, pair_region2);
                    REprintf("Error! in BASE PLAN TYPE NON ADMIN ADJACENT HIER MERGE VALID!\n");
                    throw Rcpp::exception("Bro!]n");
                }

                if(DEBUG_LOG_LINK_EDGE_VERBOSE){
                    REprintf("Region Pair (%d, %d) | %d | %d edges!\n ", 
                        pair_region1, pair_region2, 
                        pair_val.admin_adjacent,
                        edges);
                }

                auto reindexed_pair_region1 = region_reindex_vec[pair_region1];
                auto reindexed_pair_region2 = region_reindex_vec[pair_region2];
                // Now for both pairs we need to 
                // 1. Update the degree of the vertex
                // 2. Increase the count of edges between the regions 

                // But if either region has index num_regions-1
                // we don't count its degree or care about edge count

                // Check if neither pair is id num_regions-1
                if(reindexed_pair_region1 != num_component_regions-1 && 
                   reindexed_pair_region2 != num_component_regions-1){
                    // increase the degree of both vertices 
                    laplacian_minor(reindexed_pair_region1, reindexed_pair_region1) += edges;
                    laplacian_minor(reindexed_pair_region2, reindexed_pair_region2) += edges;
                    // subtract edges between them 
                    laplacian_minor(reindexed_pair_region1, reindexed_pair_region2) -= edges;
                    laplacian_minor(reindexed_pair_region2, reindexed_pair_region1) -= edges;
                }else if(reindexed_pair_region1 != num_component_regions-1){
                    // increase the degree 
                    laplacian_minor(reindexed_pair_region1, reindexed_pair_region1) += edges;
                    // because other region is merged one ignore 
                }else if(reindexed_pair_region2 != num_component_regions-1){
                    // increase the degree 
                    laplacian_minor(reindexed_pair_region2, reindexed_pair_region2) += edges;
                    // because other region is merged one ignore 
                }
                // increase current index by 1
                ++curr_index;
            }
            if(DEBUG_LOG_LINK_EDGE_VERBOSE){
            REprintf("Printing Laplacian Minor!\n");
            laplacian_minor.print();
            }
            // Now add log det
            log_tau += arma::log_det_sympd(laplacian_minor);
        }
    }

    if(curr_index >= all_pairs.size() && merged_num_admin_connected_components > 1){
        REprintf("ERROR!!\n %u pairs, curr index %d but %d num merged admin components\n",
            all_pairs.size(), curr_index, merged_num_admin_connected_components);
            Rprint();
    REprintf("%d Components | NOW SORTED Pairs: \n", merged_num_admin_connected_components);
    for(auto const &val: all_pairs){
        REprintf("Regions (%u, %u) | Components (%d, %d) | Shared Status %s\n",
            val.first.first, val.first.second, 
            county_component[val.first.first], county_component[val.first.second],
            (val.second.same_admin_component ? "SHARED" : "NOT SHARED"));
    }
    for (int i = 0; i < merged_num_admin_connected_components; i++)
    {
        REprintf("Component %d - %d Regions\n", i, component_region_counts[i]);
    }
        throw Rcpp::exception("! Hier\n");
    }


    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    REprintf("Remaining Across Merged Component Pairs: \n", all_pairs.size() - curr_index - 1);
    for (size_t i = curr_index; i < all_pairs.size(); i++)
    {
        auto val = all_pairs[i];
        REprintf("Regions (%u, %u) | Components (%d, %d) | Shared Status %s | Edges %d\n",
            val.first.first, val.first.second, 
            county_component_reindex[val.first.first], county_component_reindex[val.first.second],
            (val.second.same_admin_component ? "SHARED" : "NOT SHARED"),
            (val.second.same_admin_component ? val.second.within_county_edges : val.second.across_county_edges));
    }
    
    }

    // If only one connected component then we stop 
    if(merged_num_admin_connected_components == 1){
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("1 Component in entire Merged Plan! Log tau unchanged at %f\n", log_tau);
        }
        return log_tau;
    }else if(merged_num_admin_connected_components == 2){
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("2 Components in entire Merged Plan!\n");
        }
        int across_component_edges = 0;
        // If two then we just sum the edges across and take log 
        for (int i = curr_index; i < all_pairs.size(); i++)
        {
            across_component_edges += all_pairs[i].second.across_county_edges;
        }
        // add
        log_tau += std::log(static_cast<double>(across_component_edges));
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("Adding log of %d, Log tau now at %f\n", across_component_edges, log_tau);
        }
        return log_tau;
    }

    // Now we compute spanning trees across components 
    arma::mat component_laplacian_minor(
        merged_num_admin_connected_components-1, 
        merged_num_admin_connected_components-1, 
        arma::fill::zeros);

    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
    REprintf("Current index now %d with %d pairs\n", curr_index, all_pairs.size());
    }

    // Now we don't care about specific component id so reindex things 
    // to make the merged regions 
    /* 
    we make merged reindex which essentially just reindexes all the regions to 
    the two merged regions are the new largest value and everything else gets
    potentially shifted down. For example imagine 5 regions
    and we merge 1 and 3 so the reindex vector is [0,3,1,3,2]
     */
    int const merged_component = merged_num_admin_connected_components-1;
    for (int current_reindex = 0, a_component_id = 0; a_component_id < num_county_connected_components; a_component_id++){
        // all components but region 2's get mapped to themselves
        if(a_component_id == region2_component || a_component_id == region1_component){
            county_component_reindex[a_component_id] = merged_component;
        }else{
            county_component_reindex[a_component_id] = current_reindex;
            ++current_reindex;
        }
        // REprintf("Mapping %d to %d!\n", i, merge_index_reshuffle[i]);
    }

    
    // Iterate over the remaining pairs that are across components 
    for (int i = curr_index; i < all_pairs.size(); i++)
    {
        const auto &[pair_region1, pair_region2] = all_pairs[i].first;
        const PairHashData &pair_val = all_pairs[i].second;

        int component1_id = county_component_reindex[county_component[pair_region1]];
        int component2_id = county_component_reindex[county_component[pair_region2]];

        // 
        if(pair_val.same_admin_component){
            REprintf("BIG ERROR SAME COMPONENT WHEN SHOULD BE DIFF!\n");
        }

        int edges = pair_val.across_county_edges;
        if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("Components (%d,%d) - %d edges across!\n",
            component1_id, component2_id, edges);
        }
        // Now for both pairs we need to 
        // 1. Update the degree of the vertex
        // 2. Increase the count of edges between the regions 

        // But if either region has index num_regions-1
        // we don't count its degree or care about edge count

        // Check if neither pair is id num_regions-1
        if(component1_id != merged_num_admin_connected_components-1 && 
            component2_id != merged_num_admin_connected_components-1){
            // increase the degree of both vertices 
            component_laplacian_minor(component1_id, component1_id) += edges;
            component_laplacian_minor(component2_id, component2_id) += edges;
            // subtract edges between them 
            component_laplacian_minor(component1_id, component2_id) -= edges;
            component_laplacian_minor(component2_id, component1_id) -= edges;
        }else if(component1_id != merged_num_admin_connected_components-1){
            // increase the degree 
            component_laplacian_minor(component1_id, component1_id) += edges;
            // because other region is merged one ignore 
        }else if(component2_id != merged_num_admin_connected_components-1){
            // increase the degree 
            component_laplacian_minor(component2_id, component2_id) += edges;
            // because other region is merged one ignore 
        }
    }
    
    if(DEBUG_LOG_LINK_EDGE_VERBOSE){
        REprintf("Printing Component Laplacian Minor\n");
    component_laplacian_minor.print();
    }

    // Now add log det
    log_tau += arma::log_det_sympd(component_laplacian_minor);

    return log_tau;

}

void swap_pair_maps(RegionPairHash &a, RegionPairHash &b) {
    std::swap(a.num_hashed_pairs, b.num_hashed_pairs);
    std::swap(a.values, b.values);
    std::swap(a.hashed, b.hashed);
    std::swap(a.hashed_pairs, b.hashed_pairs);
};



PlanMultigraph::PlanMultigraph(MapParams const &map_params, bool const need_to_compute_multigraph_taus):
    map_params(map_params),
    counties_on(map_params.num_counties > 1),
    vertices_visited(map_params.V),
    county_component(map_params.ndists, 0), 
    component_split_counts(counties_on ? map_params.ndists : 0, 0),
    component_region_counts(counties_on ? map_params.ndists : 0, 0),
    region_overlap_counties(counties_on ? map_params.ndists : 0),
    num_county_region_components(0),
    num_county_connected_components(0),
    other_counties_vertices(map_params.num_edges+1),
    current_county_diff_region_vertices(map_params.num_edges+1),
    current_county_region_vertices(map_params.num_edges+1),
    pair_map(map_params.ndists),
    need_to_compute_multigraph_taus(need_to_compute_multigraph_taus),
    county_component_reindex(need_to_compute_multigraph_taus ? map_params.ndists : 0, 0),
    region_reindex_vec(need_to_compute_multigraph_taus ? map_params.ndists : 0, 0),
    WAIT_laplacian_minor(0,0),
    WAIT_merged_laplacian_minor(0,0){

};

void PlanMultigraph::prep_for_calculations(int const num_regions){
    // if no multigraph tau needed or two regions just return
    if(!need_to_compute_multigraph_taus || num_regions <= 2){
        return;
    }else{
        // else resize the minor matrices 
        WAIT_laplacian_minor = arma::mat(num_regions-1, num_regions-1, arma::fill::none);
        WAIT_merged_laplacian_minor = arma::mat(num_regions-2, num_regions-2, arma::fill::none);
        return;
    }
}


void PlanMultigraph::Rprint() const{
    pair_map.Rprint(county_component);
}


void PlanMultigraph::Rprint_detailed(Plan const &plan){
    REprintf("Pair Map has %d Elements:\n", pair_map.num_hashed_pairs);
    return;
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
            // mark as true 
            component_lookup[current_component_lookup_index] = true;
        }else{
            // else we have seen this before in which case its not hierarchically connected so immediately return false
            return std::make_pair(false, -1);
        }

        // Now we traverse all neighbors in the same component 
        current_county_region_vertices.clear();
        current_county_region_vertices.push(v);

        while(!current_county_region_vertices.empty()){
            int u = current_county_region_vertices.pop();
            int u_region = plan.region_ids[u];
            int u_county = map_params.counties[u];

            // mark this as visited 
            vertices_visited[u] = true;
            // now visit each of the vertices children 
            for(auto const child_vertex: map_params.g[u]){
                // ignore if we already visited 
                if(vertices_visited[child_vertex]) continue;
                // Now only add if region and county are the same 
                if(u_region == plan.region_ids[child_vertex] && 
                   u_county == map_params.counties[child_vertex]){
                    // if same then same component mark as visited to avoid being added later  
                    vertices_visited[child_vertex] = true;
                    current_county_region_vertices.push(child_vertex);
                }
            }
        }
    }
    
    // if we've passed through the whole graph then now return true 
    return std::make_pair(true, num_county_region_components);

}


void PlanMultigraph::build_plan_non_hierarchical_multigraph(
    PlanVector const &region_ids
){
    // reset the hash map 
    pair_map.reset();

    for (int v = 0; v < map_params.V; v++)
    {
        auto v_region = region_ids[v];

        for (auto const u: map_params.g[v])
        {
            auto u_region = region_ids[u];
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
    PlanVector const &region_ids, int const num_regions
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
    // quotient graph. reset that
    num_county_connected_components = 0;

    // clear the queues 
    other_counties_vertices.clear();
    current_county_diff_region_vertices.clear();
    current_county_region_vertices.clear();

    // Now we walk through the graph 
    for (int w = 0; w < map_params.V; w++)
    {
        // skip if already visited 
        if(vertices_visited[w]) continue;

        // else we've encountered a new component 
        int num_current_county_region_components = 1;
        int num_current_counties = 1;


        // add this vertex and mark as visited 
        current_county_region_vertices.push(w);
        // assign this region to the component
        county_component[region_ids[w]] = num_county_connected_components;

        auto current_county = map_params.counties[w];
        auto current_region = region_ids[w];

        while(
            !current_county_region_vertices.empty() ||
            !current_county_diff_region_vertices.empty() ||
            !other_counties_vertices.empty()
        ){
            
            int v;

            // first see if anything in this region and county 
            if (!current_county_region_vertices.empty()) {
                v = current_county_region_vertices.pop();
            } else if (!current_county_diff_region_vertices.empty()) {
                // else we've got a vertex in a different region in the same county
                v = current_county_diff_region_vertices.pop();
                // if we haven't visited this vertex yet that means we're reaching 
                // this component for the first time 
                if(!vertices_visited[v]){
                    // then change the current region and increase the count 
                    current_region = region_ids[v];
                    num_current_county_region_components++;
                    // increase global count 
                    num_county_region_components++;
                    // mark the region as being in current component
                    county_component[region_ids[v]] = num_county_connected_components;
                    // REprintf("Assigned Region %u to component %d\n",
                    //     plan.region_ids[v]+1, 
                    // num_county_connected_components);
                }
            } else if (!other_counties_vertices.empty()) {
                // else we've got a vertex in a different county
                v = other_counties_vertices.pop();
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
                    county_component[region_ids[v]] = num_county_connected_components;
                    // REprintf("Assigned Region %u to component %d\n",
                    //     plan.region_ids[v]+1, 
                    // num_county_connected_components);
                }
            }

            // skip if we already visited 
            if (vertices_visited[v]) continue;
            // mark as visited now 
            vertices_visited[v] = true;


            auto v_region = region_ids[v];
            auto v_county = map_params.counties[v];

            for (int u : map_params.g[v]) {
                // get region and county 
                auto u_region = region_ids[u];
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
        component_split_counts[num_county_connected_components] = num_current_county_region_components - num_current_counties;

        if(DEBUG_BASE_PLANS_VERBOSE){
            REprintf("Component %d has %d components, %d counties so %d splits!\n",
                num_county_connected_components, 
                num_current_county_region_components,
                num_current_counties,
                num_current_county_region_components - num_current_counties);
        }


        // now increment the counter by 1
        num_county_connected_components++;

        // if too many global splits then auto reject
        if(num_county_region_components - map_params.num_counties >= num_regions){
            return false;
        }
    }

    // now for each component add up the number of regions in the component
    for (size_t region_id = 0; region_id < num_regions; region_id++)
    {
        ++component_region_counts[county_component[region_id]];
    }

    // now check for each component the number of splits is the number of regions
    // minus 1 
    for (size_t component_id = 0; component_id < num_county_connected_components; component_id++)
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
        // update to say same componet 
        pair_map.values[hash_index].same_admin_component = true;
        // if admin adjacent do nothing 
        if(pair_map.values[hash_index].admin_adjacent) continue;
        // else reset values mark merge as false 
        pair_map.values[hash_index].merge_is_hier_valid = false;
        --pair_map.num_hier_smc_merge_valid_pairs;
    }


    return true;    
    
}

bool PlanMultigraph::build_plan_multigraph(
    PlanVector const &region_ids, int const num_regions
){
    if(counties_on){
        return build_plan_hierarchical_multigraph(region_ids, num_regions);
    }else{
        build_plan_non_hierarchical_multigraph(region_ids);
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


void PlanMultigraph::remove_invalid_hard_constraint_pairs(
    Plan const &plan, ScoringFunction const &scoring_function, bool const is_final_split
){
    // do nothing if no hard constraints 
    if (!scoring_function.any_hard_constraints) return;


    pair_map.hashed_pairs.erase(
        std::remove_if(pair_map.hashed_pairs.begin(), pair_map.hashed_pairs.end(), 
        [&](std::pair<RegionID, RegionID> a_pair) { 
            // check if merging the region invalidates hard constraints
            bool failed_constraint = !scoring_function.merged_plan_ok(
                plan, a_pair.first, a_pair.second, is_final_split
            );

            // if invalid then reset data in pair map 
            if(failed_constraint){
                auto hash_index = pair_map.pair_hash(a_pair.first, a_pair.second);
                // else reset values and decrease hash table count 
                --pair_map.num_hashed_pairs;
                pair_map.hashed[hash_index] = false;
                pair_map.values[hash_index] = PairHashData();
            }

            // remove if doesn't satisfy hard constraint 
            return failed_constraint;
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
    return build_plan_hierarchical_multigraph(plan.region_ids, plan.num_regions);
}


void swap_plan_multigraphs(PlanMultigraph &a, PlanMultigraph &b) {
    std::swap(a.county_component, b.county_component);
    std::swap(a.component_split_counts, b.component_split_counts);
    std::swap(a.component_region_counts, b.component_region_counts);
    std::swap(a.region_overlap_counties, b.region_overlap_counties);
    std::swap(a.num_county_region_components, b.num_county_region_components);
    std::swap(a.num_county_connected_components, b.num_county_connected_components);
    swap_pair_maps(a.pair_map, b.pair_map);
}





constexpr bool MERGED_TREE_SPLITTING_VERBOSE = false; // Compile-time constant

std::vector<EdgeCut> TreeSplitter::get_all_valid_pop_edge_cuts_in_directed_tree(
    const MapParams &map_params, 
    Tree const &ust, const int root, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    int const region_population, int const region_size,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try
) const {

    // reset pops_below_vertex and valid edges thing
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    std::fill(no_valid_edges_vertices.begin(), no_valid_edges_vertices.end(), false);
    std::vector<EdgeCut> valid_edges = get_all_valid_edges_in_directed_tree(
        ust, root, map_params.pop, stack,
        pops_below_vertex, no_valid_edges_vertices,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        region_population, region_size,
        map_params.lower, map_params.upper, map_params.target);


    return valid_edges;
}


std::pair<bool, EdgeCut> TreeSplitter::attempt_to_find_edge_to_cut(
    const MapParams &map_params, ScoringFunction const &scoring_function, RNGState &rng_state,
    Plan const &plan, int const split_region1, int const split_region2,
    Tree const &ust, const int root, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    int const region_population, int const region_size,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    bool save_selection_prob
) {
    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, 
        ust, root, stack,
        pops_below_vertex, no_valid_edges_vertices,
        region_population, region_size,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try
    );

    int num_valid_edges  = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if(num_valid_edges == 0){
        return std::make_pair(false, EdgeCut());
    }else{ // else have derived class choose according to its rule
        return select_edge_to_cut(scoring_function, ust, rng_state, valid_edges, save_selection_prob);
    }
}

// returns edge cut and log probability it was chosen
std::pair<bool, EdgeCut> TreeSplitter::select_edge_to_cut(
        ScoringFunction const &scoring_function, Tree const &ust,
        RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
        bool save_selection_prob
    ) const {
    auto num_valid_edges = valid_edges.size();

    // if no valid edges reject immediately 
    if(num_valid_edges == 1){
        // if only 1 just return that
        // selection prob is just 1 so don't touch
        // if(save_selection_prob){
        //     Rprintf("Save true: %d valid, only 1 edge, log prob is %f \n", 
        //         num_valid_edges, valid_edges[0].log_prob);
        // }
        return std::make_pair(true, valid_edges[0]);
    }

    // get the weights 
    arma::vec unnormalized_wgts(num_valid_edges);

    for (size_t i = 0; i < num_valid_edges; i++)
    {
        unnormalized_wgts(i) = compute_unnormalized_edge_cut_weight(
            valid_edges[i]
        );
    }
    
    
    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    double log_selection_prob = 0.0;
    if(save_selection_prob){
        selected_edge_cut.log_prob = std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
        // Rprintf("Save, %d valid, log prob is %f and %f\n", num_valid_edges, selected_edge_cut.log_prob, 
        //     std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts)));
    }

    return std::make_pair(true, selected_edge_cut);
}


// Takes a vector of valid edge cuts and returns the log probability 
    // the one an index idx would have been chosen 
double TreeSplitter::get_log_selection_prob(
        std::vector<EdgeCut> &valid_edges,
        int idx
) const{
    auto num_valid_edges = valid_edges.size();
    // get the weights 
    double weight_sum = 0.0;
    // get idx weight
    double idx_weight = compute_unnormalized_edge_cut_weight(valid_edges[idx]);

    // get sum of weights 
    for (size_t i = 0; i < num_valid_edges; i++)
    {
        weight_sum += compute_unnormalized_edge_cut_weight(
            valid_edges[i]
        );
    }

    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return std::log(idx_weight) - std::log(weight_sum);
}


double TreeSplitter::get_log_retroactive_splitting_prob_for_joined_tree(
    MapParams const &map_params, ScoringFunction const &scoring_function,
    VertexGraph const &forest_graph, TreePopStack &stack,
    std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
    const int region1_root, const int region2_root,
    Plan const &plan,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try
){
    const int region1_population = plan.region_pops[plan.region_ids[region1_root]];
    const int region2_population = plan.region_pops[plan.region_ids[region2_root]];
    
    const int region1_size = plan.region_sizes[plan.region_ids[region1_root]];
    const int region2_size = plan.region_sizes[plan.region_ids[region2_root]];
    int total_merged_region_size = region1_size+region2_size;


    // Get all the valid edges in the joined tree 
    std::vector<EdgeCut> valid_edges = get_valid_edges_in_joined_tree(
        map_params, forest_graph, stack,
        pops_below_vertex, visited,
        region1_root, region1_population,
        region2_root, region2_population,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_size
    );


    // find the index of the actual edge we cut 
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(
        region1_root, region2_root, region1_root, 
        region2_size, region2_population,
        region1_size, region1_population
    );

    if(MERGED_TREE_SPLITTING_VERBOSE){
    REprintf("Finding Merge prob for (%d, %d) - %u valid edges!\n", 
        region1_root, region2_root, valid_edges.size());
    }

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    int actual_cut_edge_index = std::distance(valid_edges.begin(), it);
    if(MERGED_TREE_SPLITTING_VERBOSE){
    REprintf("Actual Cut Edge at Index %d and so prob is %f \n", 
        actual_cut_edge_index,
        get_log_selection_prob(valid_edges, actual_cut_edge_index));
    }

    return get_log_selection_prob(valid_edges, actual_cut_edge_index);
}


void NaiveTopKSplitter::update_single_int_param(int int_param){
    if(int_param <= 0) throw Rcpp::exception("Splitting k must be at least 1!\n");
    k_param = int_param;
}

std::pair<bool, EdgeCut> NaiveTopKSplitter::select_edge_to_cut(
    ScoringFunction const &scoring_function, Tree const &ust,
    RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
    bool save_selection_prob
) const{

    int num_valid_edges  = static_cast<int>(valid_edges.size());
    // if(num_valid_edges > k_param){
    //     REprintf("k was %d but found %d valid edges\n", k_param, num_valid_edges);
    //     // throw Rcpp::exception("K not big enough!\n");
    // }

    int idx = rng_state.r_int(k_param);
    // if we selected k greater than number of edges failure
    if(idx >= num_valid_edges){
        return std::make_pair(false, EdgeCut()); 
    }else{
        // we always store selection probability since its so cheap to compute
        EdgeCut selected_edge_cut = valid_edges[idx];
        selected_edge_cut.log_prob = - std::log(k_param);
        return std::make_pair(true, selected_edge_cut);
    }

}



std::pair<bool, EdgeCut> UniformValidSplitter::select_edge_to_cut(
    ScoringFunction const &scoring_function, Tree const &ust,
    RNGState &rng_state,std::vector<EdgeCut> &valid_edges,
    bool save_selection_prob
) const{
    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if only 1 edge just return that
    if(num_valid_edges == 1) return std::make_pair(true, valid_edges[0]);

    // pick one unif at random 
    int idx = rng_state.r_int(num_valid_edges);
    // we always store selection probability since its so cheap to compute
    EdgeCut selected_edge_cut = valid_edges[idx];
    selected_edge_cut.log_prob = - std::log(num_valid_edges);

    return std::make_pair(true, selected_edge_cut);
}



double ExpoWeightedSplitter::compute_unnormalized_edge_cut_weight(
    EdgeCut const &edge_cut
) const{
    std::array<double, 2> devs = edge_cut.compute_abs_pop_deviances(target);
    double bigger_dev = std::max(devs.at(0), devs.at(1));
    return std::exp(-alpha*bigger_dev);
}


double ExpoWeightedSmallerDevSplitter::compute_unnormalized_edge_cut_weight(
    EdgeCut const &edge_cut
) const{
    std::array<double, 2> devs = edge_cut.compute_abs_pop_deviances(target);
    double smaller_dev = std::min(devs.at(0), devs.at(1));
    return std::exp(-alpha*smaller_dev);
}


double PopTemperSplitter::compute_unnormalized_edge_cut_weight(
    EdgeCut const &edge_cut
) const{
    double region1_pop_temper = compute_log_pop_temper(target, pop_temper, ndists,
        edge_cut.cut_above_pop, edge_cut.cut_above_region_size
    );
    double region2_pop_temper = compute_log_pop_temper(target, pop_temper, ndists,
        edge_cut.cut_below_pop, edge_cut.cut_below_region_size
    );
    // larger population deviation means bigger pop temper but we want smaller 
    // so we add then do exp(- sum) 
    return std::exp(-(region1_pop_temper +region2_pop_temper));

}


std::pair<bool, EdgeCut> ExperimentalSplitter::select_edge_to_cut(
        ScoringFunction const &scoring_function, Tree const &ust,
        RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
        bool save_selection_prob
    ) const {
    auto num_valid_edges = valid_edges.size();

    // if no valid edges reject immediately 
    if(num_valid_edges == 1){
        // if only 1 just return that
        return std::make_pair(true, valid_edges[0]);
    }

    // get the weights 
    arma::vec unnormalized_wgts = compute_almost_best_weights_on_smaller_dev_edges(
        valid_edges, epsilon, target);
    
    
    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    double log_selection_prob = 0.0;
    if(save_selection_prob){
        selected_edge_cut.log_prob = std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
    }

    return std::make_pair(true, selected_edge_cut);
}


double ExperimentalSplitter::get_log_selection_prob(
    std::vector<EdgeCut> &valid_edges,
    int idx
    ) const{
    // get the weights 
    arma::vec unnormalized_wgts = compute_almost_best_weights_on_smaller_dev_edges(
        valid_edges, epsilon, target);
    
    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return log(unnormalized_wgts(idx)) - log(arma::sum(unnormalized_wgts));
}


// std::pair<bool, EdgeCut> ConstraintSplitter::select_edge_to_cut(
//     ScoringFunction const &scoring_function, Tree const &ust,
//     RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
//     bool save_selection_prob
// ){
//     int num_valid_edges = static_cast<int>(valid_edges.size());
//     // if only 1 edge just return that
//     if(num_valid_edges == 1) return std::make_pair(true, valid_edges[0]);


//     // get the weights 
//     arma::vec unnormalized_wgts = compute_soft_constraint_edge_cut_weights(
//         valid_edges, scoring_function, ust,
//         region_ids, region_sizes, region_pops,
//         int const split_region_id1, int const split_region_id2
//     )
    
    
//     // select with prob proportional to the weights
//     int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
//     EdgeCut selected_edge_cut = valid_edges.at(idx);
//     // compute selection probability if needed
//     double log_selection_prob = 0.0;
//     if(save_selection_prob){
//         selected_edge_cut.log_prob = std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
//     }

//     return std::make_pair(true, selected_edge_cut);
// }


std::pair<bool, EdgeCut> ConstraintSplitter::attempt_to_find_edge_to_cut(
        const MapParams &map_params, ScoringFunction const &scoring_function, RNGState &rng_state,
        Plan const &plan, int const split_region1, int const split_region2,
        Tree const &ust, const int root, TreePopStack &stack,
        std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
        int const region_population, int const region_size,
        const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try,
        bool save_selection_prob
){
    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, 
        ust, root, stack,
        pops_below_vertex, no_valid_edges_vertices,
        region_population, region_size,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try
    );

    int num_valid_edges  = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if(num_valid_edges == 0){
        return std::make_pair(false, EdgeCut());
    }else if(num_valid_edges == 1){
        return std::make_pair(true, valid_edges[0]);
    }

    // copy over the current plan information 
    region_ids.copy(plan.region_ids);
    // copy the region sizes vector
    region_sizes.copy(plan.region_sizes);
    // copy population
    region_pops.copy(plan.region_pops);


    // get the weights 
    arma::vec unnormalized_wgts = compute_soft_constraint_edge_cut_weights(
        valid_edges, scoring_function, ust, plan.num_regions + 1,
        region_ids, region_sizes, region_pops,
        split_region1, split_region2, vertex_queue
    );
    
    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    double log_selection_prob = 0.0;
    if(save_selection_prob){
        selected_edge_cut.log_prob = std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
        // REprintf("Selection prob %f\n", selected_edge_cut.log_prob);
    }

    return std::make_pair(true, selected_edge_cut);
}



// assumes two trees in spanning forest have been joined
void assign_region_ids_from_joined_undirected_tree(
    VertexGraph const &forest_graph, PlanVector &region_ids,
    int const cut_vertex_root, int const cut_vertex_root_region_id,
    int const cut_vertex_parent, int const cut_parent_region_id,
    CircularQueue<std::pair<int,int>> &vertex_queue
){
    // clear the queue
    vertex_queue.clear();
    
    // Since tree is undirected we first start from cut_vertex_root
    // and just make sure to skip the parent 

    // update root and add its children to queue 
    region_ids[cut_vertex_root] = cut_vertex_root_region_id;
    for(auto const &child_vertex: forest_graph[cut_vertex_root]){
        // ignore if its the parent aka the cut edge
        if(child_vertex == cut_vertex_parent) continue;

        vertex_queue.push({child_vertex, cut_vertex_root});
    }

    // update all the children
    while(!vertex_queue.empty()){
        // get and remove head of queue 
        auto [vertex, vtx_parent] = vertex_queue.pop();
        // update region ids
        region_ids[vertex] = cut_vertex_root_region_id;
        // add children 
        for(auto const &child_vertex: forest_graph[vertex]){
            // if its the parent then skip it 
            if(child_vertex == vtx_parent) continue;
            vertex_queue.push({child_vertex, vertex});
        }
    }

    // now we update starting at the cut root vertex
    // update root and add its children to queue 
    region_ids[cut_vertex_parent] = cut_parent_region_id;
    for(auto const &child_vertex: forest_graph[cut_vertex_parent]){
        // ignore if its the child aka the cut edge
        if(child_vertex == cut_vertex_root) continue;

        vertex_queue.push({child_vertex, cut_vertex_parent});
    }

    // update all the children
    while(!vertex_queue.empty()){
        // get and remove head of queue 
        auto [vertex, parent_vtx] = vertex_queue.pop();
        // update region ids
        region_ids[vertex] = cut_parent_region_id;
        // add children 
        for(auto const &child_vertex: forest_graph[vertex]){
            // if its the parent then skip it 
            if(child_vertex == parent_vtx) continue;
            vertex_queue.push({child_vertex, vertex});
        }
    }

    return;
}


double ConstraintSplitter::get_log_retroactive_splitting_prob_for_joined_tree(
        MapParams const &map_params, ScoringFunction const &scoring_function,
        VertexGraph const &forest_graph, TreePopStack &stack,
        std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
        const int region1_root, const int region2_root,
        Plan const &plan,
        const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try
){
    const int region1_population = plan.region_pops[plan.region_ids[region1_root]];
    const int region2_population = plan.region_pops[plan.region_ids[region2_root]];
    
    const int region1_size = plan.region_sizes[plan.region_ids[region1_root]];
    const int region2_size = plan.region_sizes[plan.region_ids[region2_root]];
    int total_merged_region_size = region1_size+region2_size;

    auto const region1_id = plan.region_ids[region1_root];
    auto const region2_id = plan.region_ids[region2_root];


    // Get all the valid edges in the joined tree 
    std::vector<EdgeCut> valid_edges = get_valid_edges_in_joined_tree(
        map_params, forest_graph, stack,
        pops_below_vertex, visited,
        region1_root, region1_population,
        region2_root, region2_population,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_size
    );

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if only 1 valid edge then its log(1) = 0
    if(num_valid_edges == 1){
        return 0.0;
    }


    // find the index of the actual edge we cut 
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(
        region1_root, region2_root, region1_root, 
        region2_size, region2_population,
        region1_size, region1_population
    );

    if(MERGED_TREE_SPLITTING_VERBOSE){
    REprintf("Finding Merge prob for (%d, %d) - %u valid edges!\n", 
        region1_root, region2_root, valid_edges.size());
    }

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    // copy over the current plan information 
    region_ids.copy(plan.region_ids);
    // copy the region sizes vector
    region_sizes.copy(plan.region_sizes);
    // copy population
    region_pops.copy(plan.region_pops);

    int actual_cut_edge_index = std::distance(valid_edges.begin(), it);
    // copy the forest over
    dummy_forest = forest_graph;
    // add the actual removed edge back 
    dummy_forest[region1_root].push_back(region2_root);
    dummy_forest[region2_root].push_back(region1_root);


    std::vector<long double> unnormed_wgts;
    unnormed_wgts.reserve(valid_edges.size());

    for (size_t i = 0; i < valid_edges.size(); i++)
    {
        // update split info 
        region_sizes[region1_id] = valid_edges[i].cut_above_region_size;
        region_sizes[region2_id] = valid_edges[i].cut_below_region_size;

        region_pops[region1_id] = valid_edges[i].cut_above_pop;
        region_pops[region2_id] = valid_edges[i].cut_below_pop;
        
        // update the region ids 
        assign_region_ids_from_joined_undirected_tree(
            forest_graph, region_ids,
            valid_edges[i].cut_vertex, region1_id,
            valid_edges[i].cut_vertex_parent, region2_id,
            vertex_queue
        );

        // get the soft score 
        double const score = scoring_function.compute_full_split_plan_soft_score(
            plan.num_regions, region_ids, region_sizes, region_pops,
            region1_id, region2_id
        );

        unnormed_wgts.push_back(std::exp(-score));

        // REprintf("Soft score %f, unnormed weight %.30f vs  %.30f \n", score, unnormed_wgts[i], std::exp(-score));

    }

    auto sum = std::accumulate(unnormed_wgts.begin(), unnormed_wgts.end(), 0.0);
    auto log_sum = std::log(sum);
    // REprintf("Sum %.30f\n", sum);


    // compute selection probability if needed
    double log_selection_prob = std::log(unnormed_wgts[actual_cut_edge_index]) - log_sum;
    // REprintf("Actual Cut Edge at Index %d and so prob is %f \n", 
    //     actual_cut_edge_index, log_selection_prob);

    if(MERGED_TREE_SPLITTING_VERBOSE){
    REprintf("Actual Cut Edge at Index %d and so prob is %f \n", 
        actual_cut_edge_index, log_selection_prob);
    }

    return log_selection_prob;

    // get the weights
    arma::vec unnormalized_wgts(unnormed_wgts.size());

    for (size_t i = 0; i < unnormed_wgts.size(); i++)
    {
        unnormalized_wgts[i] = std::exp(
            std::log(unnormed_wgts[i]) - log_sum
        );
        REprintf("%.30f\n", unnormalized_wgts[i]);
    }

    unnormalized_wgts = arma::cumsum(unnormalized_wgts);
    
    REprintf("Weights are:\n");
    for (auto const v: unnormalized_wgts)
    {
        REprintf("%.20f, ", v);
    }
     REprintf("\n");
    


}