/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/01
 * Purpose: Implements base Plan (ie graph partition) class 
 ********************************************************/

#include "base_plan_type.h"

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

// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED

// TODO: Need to add option to pass region population

// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED
Plan::Plan(
    arma::subview_col<arma::uword> region_ids_col, 
    arma::subview_col<arma::uword> region_sizes_col, 
    int ndists, int num_regions, const arma::uvec &pop, bool split_district_only): 
    num_regions(num_regions),
    region_ids(region_ids_col.begin(), region_ids_col.end()),
    region_sizes(region_sizes_col.begin(), region_sizes_col.end())
{
    // check num_regions and num_districts inputs make sense
    if (ndists < 2) throw Rcpp::exception("Tried to create a plan with ndists < 2 regions!");
    if (region_sizes.size() != ndists) throw Rcpp::exception("The region dvals column passed in is not size ndists!");

    // set number of multidistricts, and V
    int const V = region_ids.size();
    region_order_max = ndists+1;


    // now check these
    if (num_regions > ndists) throw Rcpp::exception("Tried to create a plan object with more regions than ndists!");
    if (num_regions == 0) throw Rcpp::exception("Tried to create a plan with 0 regions");

    // Create other region-level information 
    region_added_order = std::vector<int>(ndists, -1);
    // fill first num_regions entries with 1,...,num_regions 
    std::iota(
        std::begin(region_added_order), 
        std::begin(region_added_order) + num_regions, 
        1); 
    
    region_pops = std::vector<int>(ndists, 0);
    // compute the population for each of the regions 
    for (size_t v = 0; v < V; v++)
    {
        region_pops.at(region_ids[v]) += pop(v);
    }
}


// Plan::Plan(const Plan& other)
//     : region_ids(other.region_ids), region_sizes(other.region_sizes) // Share the same reference
// {
//     // Copy simple members
//     ndists = other.ndists;
//     V = other.V;
//     num_regions = other.num_regions;
//     num_districts = other.num_districts;
//     num_multidistricts = other.num_multidistricts;
//     map_pop = other.map_pop;
//     remainder_region = other.remainder_region;

//     // Deep copy std::vector members
//     region_pops = other.region_pops;
//     region_added_order = other.region_added_order;
//     region_order_max = other.region_order_max;
// }


// Plan::Plan(const Plan& other){
//     return *this;
// }



// returns the region ids of the two most recently split regions
// DO NOT CALL THIS WHEN THERE IS ONLY 1 REGION!
std::pair<int, int> Plan::get_most_recently_split_regions() const{
    int largest_index = -1, second_largest_index = -1;
    int largest_value = -1 * num_regions, second_largest_value = -1* num_regions;

    // Iterate through the vector
    for (std::size_t i = 0; i < num_regions; ++i) {
        if (region_added_order.at(i) > largest_value) {
            // Update second-largest before updating largest
            second_largest_value = largest_value;
            second_largest_index = largest_index;

            largest_value = region_added_order.at(i);
            largest_index = i;
        } else if (region_added_order.at(i) > second_largest_value) {
            // Update second-largest only
            second_largest_value = region_added_order.at(i);
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

    dummy_plan.region_ids = region_ids;
    // First we relabel all the region vertex ids
    for (size_t i = 0; i < this->region_ids.size(); i++)
    {
        // Send region id i to old_region_id_to_new_vec[i]
        this->region_ids[i] = old_region_id_to_new_vec.at(dummy_plan.region_ids[i]);
    }


    
    // Make some dummy plan attributes a shallow copy of the plan
    dummy_plan.region_sizes = region_sizes;
    dummy_plan.region_pops = region_pops;

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
        this->region_added_order.at(i) = i+1;
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
void Plan::Rprint() const{
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
}


/*  
 * @title Count the Number of Valid Adjacent Region Pairs
 * 
 * Counts and returns the number of valid adjacent region pairs in a graph. 
 * We define a valid of pair of adjacent regions to be two regions that 
 *      1. Share at least one edge between them in `g`
 *      2. The regions sizes are valid to merge with eligibility determined
 *         by the `valid_merge_pairs` matrix
 * 
 * @param g A graph (adjacency list) passed by reference
 * @param plan A plan object
 * @param check_adj_to_regions A vector tracking whether or not we should 
 * check for edges adjacent to vertices in a region of a particular size. For
 * example, `check_adj_to_regions[i] == true` means we will attempt to find 
 * edges adjacent to any vertex in a region of size i. This vector is 1-indexed
 * meaning we don't need to subtract region size by 1 when accessing.
 * @param valid_merge_pairs A 2D `ndists+1` by `ndists+1` boolean matrix that
 * uses region size to check whether or not two regions are considered a valid
 * merge that can be counted in the map. For example `valid_merge_pairs[i][j]`
 * being true means that any regions where the sizes are (i,j) are considered
 * valid to merge. 2D matrix is 1 indexed (don't need to subtract region size)
 * @param existing_pair_map A hash map mapping pairs of regions to a double. 
 * This is an optional parameter and if its empty then its ignored. If its not
 * empty then pairs already in the hash won't be counted added to the output.
 * 
 * @details No modifications to inputs made
 * @return The number of valid adjacent region pairs
 */
int Plan::count_valid_adj_regions(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule
) const{
    int num_valid_adj_pairs = get_valid_adj_regions(
        map_params, splitting_schedule
    ).size();

    return num_valid_adj_pairs;
}


/*  
 * @title Count the Number of Valid Adjacent Region Pairs
 * 
 * Counts and returns the number of valid adjacent region pairs in a graph. 
 * We define a valid of pair of adjacent regions to be two regions that 
 *      1. Share at least one edge between them in `g`
 *      2. The regions sizes are valid to merge with eligibility determined
 *         by the `valid_merge_pairs` matrix
 * 
 * @param g A graph (adjacency list) passed by reference
 * @param plan A plan object
 * @param check_adj_to_regions A vector tracking whether or not we should 
 * check for edges adjacent to vertices in a region of a particular size. For
 * example, `check_adj_to_regions[i] == true` means we will attempt to find 
 * edges adjacent to any vertex in a region of size i. This vector is 1-indexed
 * meaning we don't need to subtract region size by 1 when accessing.
 * @param valid_merge_pairs A 2D `ndists+1` by `ndists+1` boolean matrix that
 * uses region size to check whether or not two regions are considered a valid
 * merge that can be counted in the map. For example `valid_merge_pairs[i][j]`
 * being true means that any regions where the sizes are (i,j) are considered
 * valid to merge. 2D matrix is 1 indexed (don't need to subtract region size)
 * @param existing_pair_map A hash map mapping pairs of regions to a double. 
 * This is an optional parameter and if its empty then its ignored. If its not
 * empty then pairs already in the hash won't be counted added to the output.
 * 
 * @details No modifications to inputs made
 * @return The number of valid adjacent region pairs
 */
std::vector<std::pair<int,int>> Plan::get_valid_adj_regions(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    bool const check_split_constraint
) const{
    // 2D matrix for tracking if regions are adjacent 
    // To save space we will only ever access with [smaller id][larger id]
    // So it will be `num_regions-1` rows and row i will have num_regions-1
    // elements and the second element will need to be 1 indexed so 
    // subtract off 1 when accessing
    std::vector<std::vector<bool>> is_valid_adj_pair(num_regions-1);

    std::vector<std::pair<int,int>> valid_adj_pairs;

    for (size_t i = 0; i < num_regions-1; i++)
    {
        // for smaller value i there are num_regions-(i+1) larger values 
        is_valid_adj_pair[i] = std::vector<bool>(num_regions-(i+1), false);
    }


    // TEMP
    std::vector<bool> visited(map_params.num_counties > 1 ? map_params.V : 0);

    int const V = map_params.V;

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = region_ids[v];
        auto v_region_size = region_sizes[v_region_num];

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!splitting_schedule.check_adj_to_regions.at(v_region_size)){
            continue;
        }

        // get neighbors
        std::vector<int> nbors = map_params.g[v];

        // now iterate over its neighbors
        for (int v_nbor : nbors) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = region_ids[v_nbor];

            // ignore if they are in the same region
            if(v_region_num == v_nbor_region_num) continue;
            // ignore if pair can't be merged
            auto v_nbor_region_size = region_sizes[v_nbor_region_num];
            if(!splitting_schedule.valid_merge_pair_sizes[v_region_size][v_nbor_region_size]) continue;

            // else they are different so regions are adj
            // check if we need to be aware of double counting. IE if both regions
            // are ones where check adjacent is true then its possible to double count edges
            bool double_counting_not_possible = !splitting_schedule.check_adj_to_regions[v_nbor_region_size];
            // Rprintf("Neighbor size is %d and Double counting is %d\n", 
            //     (int) v_nbor_region_size, (int) double_counting_not_possible);

            // Now add this edge if either we can't double count it 
            // or `v_region_num` is the smaller of the pair 
            if(double_counting_not_possible || v_region_num < v_nbor_region_num){
                // Rprintf("Counting (%d,%d)\n", v, v_nbor);
                int smaller_region_id  = std::min(v_nbor_region_num,v_region_num);
                int bigger_region_id  = std::max(v_nbor_region_num,v_region_num);

                // Check if we've counted this before 
                if(!is_valid_adj_pair[smaller_region_id][bigger_region_id-1]){
                    // If not increase the count by one and mark this as true 
                    is_valid_adj_pair[smaller_region_id][bigger_region_id-1] = true;
                    // if we care about merging or more than one county need to check which to remove
                    if(!check_split_constraint || map_params.num_counties <= 1){
                        // add the pair 
                        valid_adj_pairs.emplace_back(smaller_region_id, bigger_region_id);
                    }else{
                        int merged_plan_county_splits = count_merged_county_splits(map_params, visited,
                            smaller_region_id, bigger_region_id);
                        // only count if merging has num_regions - 2 splits or less
                        if(merged_plan_county_splits <= num_regions - 2){
                            valid_adj_pairs.emplace_back(smaller_region_id, bigger_region_id);
                        }
                    }
                }
            }
        }
    }

     

    return valid_adj_pairs;
}


// gets adjacent pairs ignoring county restrictions
std::pair<int, std::vector<std::pair<int,int>>> Plan::get_or_count_valid_adj_regions_ignore_counties(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    bool const count_only
) const{
    // 2D matrix for tracking if regions are adjacent 
    // To save space we will only ever access with [smaller id][larger id]
    // So it will be `num_regions-1` rows and row i will have num_regions-1
    // elements and the second element will need to be 1 indexed so 
    // subtract off 1 when accessing
    std::vector<std::vector<bool>> is_valid_adj_pair(num_regions-1);

    std::vector<std::pair<int,int>> valid_adj_pairs;
    int num_valid_adj_regions = 0;

    for (size_t i = 0; i < num_regions-1; i++)
    {
        // for smaller value i there are num_regions-(i+1) larger values 
        is_valid_adj_pair[i] = std::vector<bool>(num_regions-(i+1), false);
    }

    int const V = map_params.V;

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = region_ids[v];
        auto v_region_size = region_sizes[v_region_num];

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!splitting_schedule.check_adj_to_regions.at(v_region_size)){
            continue;
        }

        // now iterate over its neighbors
        for (int v_nbor : map_params.g[v]) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = region_ids[v_nbor];

            // ignore if they are in the same region
            if(v_region_num == v_nbor_region_num) continue;
            // ignore if pair can't be merged
            auto v_nbor_region_size = region_sizes[v_nbor_region_num];
            if(!splitting_schedule.valid_merge_pair_sizes[v_region_size][v_nbor_region_size]) continue;

            // else they are different so regions are adj
            // check if we need to be aware of double counting. IE if both regions
            // are ones where check adjacent is true then its possible to double count edges
            bool double_counting_not_possible = !splitting_schedule.check_adj_to_regions[v_nbor_region_size];
            // Rprintf("Neighbor size is %d and Double counting is %d\n", 
            //     (int) v_nbor_region_size, (int) double_counting_not_possible);

            // Now add this edge if either we can't double count it 
            // or `v_region_num` is the smaller of the pair 
            if(double_counting_not_possible || v_region_num < v_nbor_region_num){
                // Rprintf("Counting (%d,%d)\n", v, v_nbor);
                int smaller_region_id  = std::min(v_nbor_region_num,v_region_num);
                int bigger_region_id  = std::max(v_nbor_region_num,v_region_num);

                // Check if we've counted this before 
                if(!is_valid_adj_pair[smaller_region_id][bigger_region_id-1]){
                    // If not increase the count by one and mark this as true 
                    is_valid_adj_pair[smaller_region_id][bigger_region_id-1] = true;
                    ++num_valid_adj_regions;
                    // if we want the specific pairs then we also add that
                    if(!count_only){
                        // add the pair 
                        valid_adj_pairs.emplace_back(smaller_region_id, bigger_region_id);
                    }
                }
            }
        }
    }

    return std::make_pair(num_valid_adj_regions, valid_adj_pairs);
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
    MapParams const &map_params
) const{
    return compute_log_region_multigraph_spanning_tree(
        build_region_multigraph(map_params.g, region_ids, num_regions)
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
//' Takes a spanning tree `ust` drawn on a specific region and attempts to cut
//' it to produce two new regions using the generalized splitting procedure
//' outlined <PAPER HERE>. This function is based on `cut_districts` in `smc.cpp`
//' however the crucial difference is even if a cut is successful it does not
//' update the plan. Instead it just returns the information on the two new
//' regions if successful and the cut tree.
//'
//' It will only attempt to create regions where the size is between
//' min_potential_d and max_potential_d (inclusive). So the one district
//' split case is `min_potential_d=max_potential_d=1`.
//'
//' By convention the first new region (`new_region1`) will always be the region
//' with the smaller d-value (although they can be equal).
//'
//' @title Attempt to Find a Valid Spanning Tree Edge to Cut into Two New Regions 
//'
//' @param g A graph (adjacency list) passed by reference
//' @param ust A directed tree object (this will be cleared in the function so
//' whatever it was before doesn't matter)
//' @param counties Vector of county labels of each vertex in `g`
//' @param cg multigraph object 
//' @param region_to_draw_tree_on The id of the region to draw the tree on
//' @param visited A vector of length `V` that will be used to draw the tree.
//' This will be cleared in the function so whatever was in it before doesn't 
//' matter
//' @param ignore A vector of length `V` that will be used to draw the tree.
//' This will be cleared in the function so whatever was in it before doesn't 
//' matter
//' @param pop A vector of the population associated with each vertex in `g`
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param root Will be set to the root of the tree drawn
//'
//' @details Modifications
//'    - `ust` is modified in place to be the new tree
//'    - `root` is set to be the root of the tree drawn
//'
//' @return `true` if tree was successfully draw, `false` otherwise.
//'
//' @keyword internal
//' @noRd
bool Plan::draw_tree_on_region(const MapParams &map_params, const int region_to_draw_tree_on,
        Tree &ust, std::vector<bool> &visited, std::vector<bool> &ignore, int &root,
        RNGState &rng_state) {

    int V = map_params.g.size();

    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = region_ids[i] != region_to_draw_tree_on;
    }

    // Get a uniform spanning tree drawn on that region
    clear_tree(ust);
    // Get a tree
    int result = sample_sub_ust(map_params.g, ust, 
        V, root, visited, ignore, 
        map_params.pop, map_params.lower, map_params.upper, 
        map_params.counties, map_params.cg,
        rng_state);
    // result == 0 means it was successful
    return(result == 0);
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
    region_added_order.at(split_region1_id) = new_region1_order_added_num;
    region_pops.at(split_region1_id) = split_region1_pop;

    // updates the new region 2
    // New region 2's id is the highest id number so push back
    region_sizes[split_region2_id] = split_region2_size;
    region_added_order.at(split_region2_id) = new_region2_order_added_num;
    region_pops.at(split_region2_id) = split_region2_pop;


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

