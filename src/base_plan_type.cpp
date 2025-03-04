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
void Plan::check_inputted_region_sizes(bool split_district_only){

    // check sum of first num_region elements is ndists and it matches expected
    // number of districts 
    int num_districts_implied_by_sizes_mat = 0;
    int total_size_implied_by_sizes_mat = 0;
    for (size_t i = 0; i < num_regions; i++)
    {
        // make sure each regions dval is non-zero
        if(region_sizes(i) <= 0){
            throw Rcpp::exception("Region size input for region is 0 or less!");
        }

        total_size_implied_by_sizes_mat += region_sizes(i);
        if(region_sizes(i) == 1) num_districts_implied_by_sizes_mat++;

        // add check that if split district only then only last one has size > 1
        if (split_district_only && i != num_regions-1)
        {
            if(region_sizes(i) != 1) throw Rcpp::exception(
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
        if(region_sizes(i) != 0){
            REprintf("Expected Region %d to be zero it was actually %d!\n",
                i, (int) region_sizes(i));
            throw Rcpp::exception(
                "Plan didn't have the correct expected region size matrix!\n"
                );
        }
    }

    num_districts = num_districts_implied_by_sizes_mat;
    num_multidistricts = num_regions - num_districts;
    return;
}


void Plan::check_inputted_region_ids(){
    // Use std::unordered_set to store unique elements
    std::unordered_set<int> unique_ids;
    // get unique labels
    for (size_t i = 0; i < region_ids.n_elem; ++i) {
        unique_ids.insert(region_ids(i)); // Insert each element of the subview column
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
    ndists(ndists), num_regions(num_regions),
    region_ids(region_ids_col), region_sizes(region_sizes_col)
{
    // check num_regions and num_districts inputs make sense
    if (ndists < 2) throw Rcpp::exception("Tried to create a plan with ndists < 2 regions!");
    if (region_sizes.n_elem != ndists) throw Rcpp::exception("The region dvals column passed in is not size ndists!");

    // set number of multidistricts, and V
    this->V = region_ids.n_elem;
    map_pop = 0;
    region_order_max = ndists+1;

    if(split_district_only){
        remainder_region = num_regions-1;
    }else{
        remainder_region = -1;
    }

    // check region sizes and labels if not just 1 region 
    if(num_regions > 1){
        check_inputted_region_ids();
        // this also sets the number of regions and multidistricts 
        check_inputted_region_sizes(split_district_only);
    }else{
        num_districts = 0;
        num_multidistricts = num_regions;
    }
    

    // now check these
    if (num_regions > ndists) throw Rcpp::exception("Tried to create a plan object with more regions than ndists!");
    if (num_districts > ndists) throw Rcpp::exception("Tried to create a plan object with more districts than ndists!");
    if (num_districts > num_regions) throw Rcpp::exception("Tried to create a plan object with more districts than total regions!");
    if (num_districts == num_regions && num_regions != ndists) throw Rcpp::exception("Tried to create a partial plan with only districts!");
    if (num_districts < 0 || num_regions < 0 || ndists < 0) throw Rcpp::exception("Tried to create a plan object with negative number of regions!");

    // Create other region-level information 
    region_added_order = std::vector<int>(ndists, -1);
    // fill first num_regions entries with 1,...,num_regions 
    std::iota(std::begin(region_added_order), std::begin(region_added_order) + num_regions, 1); 
    
    region_pops = std::vector<int>(ndists, 0);
    // compute the population for each of the regions 
    for (size_t v = 0; v < V; v++)
    {
        region_pops.at(region_ids(v)) += pop(v);
        map_pop += pop(v);
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
    int largest_value = -1 * ndists, second_largest_value = -1* ndists;

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
// After reordering it copies that over to the dummy plan as well
void Plan::reorder_plan_by_oldest_split(
    Plan &dummy_plan) {
    // Make dummy plan a shallow copy of the plan
    dummy_plan = *this;

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


    // First we relabel all the region vertex ids
    for (size_t i = 0; i < this->region_ids.n_elem; i++)
    {
        // Send region id i to old_region_id_to_new_vec[i]
        this->region_ids(i) = old_region_id_to_new_vec.at(dummy_plan.region_ids(i));
    }


    // Relabel the remainder region if needed
    if(this->remainder_region >= 0){
        this->remainder_region = old_region_id_to_new_vec.at(dummy_plan.remainder_region);
    }
    

    // Now we reorder the region dvals and population 
    for (size_t i = 0; i < this->num_regions; i++)
    {
        // Recall indices[i] == k means the old region with id k now has id i
        // so we want to set the value at the new region id `i` to the value it
        // had in the old one which is `indices[i]`
        int old_region_id = indices[i];
        int new_region_id = i;

        this->region_sizes(new_region_id) = dummy_plan.region_sizes(old_region_id);
        this->region_pops.at(new_region_id) = dummy_plan.region_pops.at(old_region_id);
        // Since the regions are in order of added reset the order added to just be 1,..., n
        this->region_added_order.at(i) = i+1;
    }

    // reset the max region counter 
    this->region_order_max = this->num_regions + 6;
    
    // copy the dummy plan over
    dummy_plan = *this;
}



// Prints our object using Rcout. Should be used in Rcpp call
void Plan::Rprint() const{
    RcppThread::Rcout << "Plan with " << num_regions << " regions, " << num_districts
                      << " districts, " << num_multidistricts << " multidistricts and "
                      << arma::sum(region_sizes) << " sum of sizes and "
                      << V << " Vertices.\n";


    RcppThread::Rcout << "Region Level Values:[";
    for(int region_id = 0; region_id < num_regions; region_id++){
        RcppThread::Rcout << "(Region " << region_id <<
            ", Size=" << region_sizes(region_id) << ", pop= " <<
            region_pops.at(region_id) <<"), ";
    }
    RcppThread::Rcout << "]\n";

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
    Graph const &g,
    std::vector<bool> const &check_adj_to_regions,
    std::vector<std::vector<bool>> const &valid_merge_pairs
) const{
    // 2D matrix for tracking if regions are adjacent 
    // To save space we will only ever access with [smaller id][larger id]
    // So it will be `num_regions-1` rows and row i will have num_regions-1
    // elements and the second element will need to be 1 indexed so 
    // subtract off 1 when accessing
    std::vector<std::vector<bool>> is_valid_adj_pair(num_regions-1);

    for (size_t i = 0; i < num_regions-1; i++)
    {
        // for smaller value i there are num_regions-(i+1) larger values 
        is_valid_adj_pair[i] = std::vector<bool>(num_regions-(i+1), false);
    }

    int num_valid_adj_pairs = 0;

    int const V = g.size();

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = region_ids(v);
        auto v_region_size = region_sizes(v_region_num);

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!check_adj_to_regions.at(v_region_size)){
            continue;
        }

        // get neighbors
        std::vector<int> nbors = g[v];

        // now iterate over its neighbors
        for (int v_nbor : nbors) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = region_ids(v_nbor);

            // ignore if they are in the same region
            if(v_region_num == v_nbor_region_num) continue;
            // ignore if pair can't be merged
            auto v_nbor_region_size = region_sizes(v_nbor_region_num);
            if(!valid_merge_pairs[v_region_size][v_nbor_region_size]) continue;

            // else they are different so regions are adj
            // check if we need to be aware of double counting. IE if both regions
            // are ones where check adjacent is true then its possible to double count edges
            bool double_counting_not_possible = !check_adj_to_regions[v_nbor_region_size];
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
                    ++num_valid_adj_pairs;
                    is_valid_adj_pair[smaller_region_id][bigger_region_id-1] = true;
                }
            }
        }
    }

    return num_valid_adj_pairs;
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
        ){


    if(num_multidistricts < 1){
        throw Rcpp::exception("ERROR: Trying to find multidistrict to split when there are none!\n");
    }
    if(num_multidistricts == 1){ // if just one multidistrict return that 
        int region_id_to_split = region_sizes.index_max();
        return region_id_to_split;
    }

    // make vectors with cumulative d value and region label for later
    std::vector<int> valid_region_ids, associated_region_sizes;
    valid_region_ids.reserve(num_multidistricts);
    associated_region_sizes.reserve(num_multidistricts);


    for(int region_id = 0 ; region_id < num_regions; region_id++) {
        auto region_size = region_sizes(region_id);

        // if valid then add id to vector 
        if(valid_region_sizes_to_split[region_size]){
            // add the count and label to vector
            valid_region_ids.push_back(region_id);
            associated_region_sizes.push_back(region_size);
        }
    }

    auto num_candidates = valid_region_ids.size();

    // pick index unif at random 
    // int idx = r_int(num_candidates);

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
        ignore[i] = region_ids(i) != region_to_draw_tree_on;
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
        EdgeCut cut_edge, bool split_district_only,
        const int split_region1_id, const int split_region2_id
){
    // Get information on the two new regions cut
    int split_region1_tree_root, split_region2_tree_root;
    int split_region1_size, split_region2_size;
    int split_region1_pop, split_region2_pop;

    cut_edge.get_split_regions_info(
        split_region1_tree_root, split_region1_size, split_region1_pop,
        split_region2_tree_root, split_region2_size, split_region2_pop
    );

    // update plan with new regions
    num_regions++; // increase region count by 1
    num_multidistricts--; // Decrease by one to avoid double counting later

    // Create info for two new districts

    // Set label and count depending on if district or multi district
    if(split_region1_size == 1){
        num_districts++;
    }else{
        num_multidistricts++;
    }

    // Now do it for second region
    if(split_region2_size == 1){
        num_districts++;
    }else{
        num_multidistricts++;
    }


    // make the first new region get max plus one
    region_order_max++;
    int new_region1_order_added_num = region_order_max;
    // Now make the second new region max plus one again
    region_order_max++;
    int new_region2_order_added_num = region_order_max;


    // Now update the region level information
    // updates the new region 1
    region_sizes(split_region1_id) = split_region1_size;
    region_added_order.at(split_region1_id) = new_region1_order_added_num;
    region_pops.at(split_region1_id) = split_region1_pop;

    // updates the new region 2
    // New region 2's id is the highest id number so push back
    region_sizes(split_region2_id) = split_region2_size;
    region_added_order.at(split_region2_id) = new_region2_order_added_num;
    region_pops.at(split_region2_id) = split_region2_pop;


    // If district split only then set the remainder region as well
    if(split_district_only){
        // update the remainder region value if needed
        remainder_region = split_region2_id;
    }
}



bool Plan::attempt_split(const MapParams &map_params, const SplittingSchedule &splitting_schedule,
                 Tree &ust, TreeSplitter &tree_splitter,
                 std::vector<bool> &visited, std::vector<bool> &ignore, 
                 RNGState &rng_state,
                 int const min_region_cut_size, int const max_region_cut_size, 
                 std::vector<int> const &smaller_cut_sizes_to_try,
                 const bool split_district_only, 
                 const int region_id_to_split, const int new_region_id)
{
    if(TREE_SPLITTING_DEBUG_VERBOSE){
    REprintf("Drawing tree on region %d which is size %d. Smallest/Biggest is (%d, %d)\n", 
    region_id_to_split, (int) region_sizes(region_id_to_split), min_region_cut_size,
    max_region_cut_size);
    int split_region_size = (int) region_sizes(region_id_to_split);
    REprintf("We can try cut sizes: ");
    for (auto const &small_size: smaller_cut_sizes_to_try)
    {
        REprintf("(%d, %d), ", small_size, split_region_size-small_size);
    }
    REprintf("\n");
    }
    

    // Now try to draw a tree on the region
    int root;
    
    // Try to draw a tree on region
    bool tree_drawn = draw_tree_on_region(map_params, region_id_to_split,
        ust, visited, ignore, root, rng_state);

    if(!tree_drawn) return false;

    // Now try to select an edge to cut
    std::pair<bool, EdgeCut> edge_search_result = tree_splitter.select_edge_to_cut(
        map_params, *this, ust, root, rng_state,
        min_region_cut_size, max_region_cut_size, 
        smaller_cut_sizes_to_try,
        region_id_to_split
    );

    // return false if unsuccessful
    if(!edge_search_result.first) return false;


    // If successful extract the edge cut info
    EdgeCut cut_edge = edge_search_result.second;
    // Now erase the cut edge in the tree
    erase_tree_edge(ust, cut_edge);


    // now update the region level information from the edge cut
    update_region_info_from_cut(
        cut_edge, split_district_only,
        region_id_to_split, new_region_id
    );

    // Now update the vertex level information
    update_vertex_info_from_cut(
        ust, cut_edge, 
        region_id_to_split, new_region_id, split_district_only
    );

    return true;
}
