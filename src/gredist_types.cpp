/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2024/08
 * Purpose: Sequential Monte Carlo redistricting sampler
 ********************************************************/

#include "gredist_types.h"


// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED


// TODO: Need to add option to pass region population

// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED
Plan::Plan(
    arma::subview_col<arma::uword> region_ids_col, 
    arma::subview_col<arma::uword> region_dvals_col, 
    int N, int total_map_pop, bool split_district_only,
    int num_regions, int num_districts,
    const arma::uvec &pop)
: region_ids(region_ids_col), region_dvals(region_dvals_col) // Initialize the reference
{
    // check num_regions and num_districts inputs make sense
    if (N < 2) throw Rcpp::exception("Tried to create a plan with N < 2 regions!");
    if (num_regions > N) throw Rcpp::exception("Tried to create a plan object with more regions than N!");
    if (num_districts > N) throw Rcpp::exception("Tried to create a plan object with more districts than N!");
    if (num_districts > num_regions) throw Rcpp::exception("Tried to create a plan object with more districts than total regions!");
    if (num_districts == num_regions && num_regions != N) throw Rcpp::exception("Tried to create a partial plan with only districts!");
    if (num_districts < 0 || num_regions < 0 || N < 0) throw Rcpp::exception("Tried to create a plan object with negative number of regions!");
    if (region_dvals.n_elem != N) throw Rcpp::exception("The region dvals column passed in is not size N!");

    // set number of regions, districts, multidistricts, map pop, and V
    this->num_regions = num_regions;
    this->num_districts = num_districts;
    num_multidistricts = num_regions - num_districts;
    this->V = region_ids.n_elem;
    map_pop = total_map_pop;
    this->N = N;
    region_order_max = N+1;

    // TODO add checks if number of regions > 1
    if(num_regions > 1){
        // check sum of first num_region elements is N and it matches expected
        // number of districts 
        int num_districts_implied_by_dval_mat = 0;
        int total_dvals_implied_by_dval_mat = 0;
        for (size_t i = 0; i < num_regions; i++)
        {
            // make sure each regions dval is non-zero
            if(region_dvals(i) <= 0) throw Rcpp::exception("Region dvals input for region is 0 or less!"); 

            total_dvals_implied_by_dval_mat += region_dvals(i);
            if(region_dvals(i) == 1) num_districts_implied_by_dval_mat++;

            // add check that if split district only then only last one has dval > 1
            if (split_district_only && i != num_regions-1)
            {
                if(region_dvals(i) != 1) throw Rcpp::exception("For partial plan the remainder does not have region id=num_regions-1!"); 
            }
        }

        if(num_districts_implied_by_dval_mat != num_districts) throw Rcpp::exception("Region dvals mat does not have the correct number of districts!"); 
        if(total_dvals_implied_by_dval_mat != N) throw Rcpp::exception("Sum of dvals in region dvals mat does equal N!"); 
        
        // If multiple regions it assumes that 0 is oldest region then 1 ... num_regions
    }else{
        if(region_dvals(0) != N) throw Rcpp::exception("Region dvals input for 1-region partial plan not N!");
    }


    // Create other region-level information 
    region_added_order = std::vector<int>(N, -1);
    // fill first num_regions entries with 1,...,num_regions 
    std::iota(std::begin(region_added_order), std::begin(region_added_order) + num_regions, 1); 
    region_pops = std::vector<int>(N, -1);
    
    if(num_regions > 1){
        // compute the population for each of the regions 
        for (size_t v = 0; v < V; v++)
        {
            region_pops.at(region_ids(v)) += pop(v);
        }
        if(split_district_only){
            remainder_region = num_regions-1;
        }else{
            remainder_region = -1;
        }
    }else{
        region_pops.at(0) = total_map_pop;
        if(split_district_only){
            remainder_region = 0;
        }else{
            remainder_region = -1;
        }
    }
    
}


Plan::Plan(const Plan& other)
    : region_ids(other.region_ids), region_dvals(other.region_dvals) // Share the same reference
{
    // Copy simple members
    N = other.N;
    V = other.V;
    num_regions = other.num_regions;
    num_districts = other.num_districts;
    num_multidistricts = other.num_multidistricts;
    map_pop = other.map_pop;
    remainder_region = other.remainder_region;

    // Deep copy std::vector members
    region_pops = other.region_pops;
    region_added_order = other.region_added_order;
    region_order_max = other.region_order_max;
}


// Performs a shallow copy of another plan meaning it 
// copies all its data but does not change the arma:: matrix it points at
void Plan::shallow_copy(const Plan& other){
    // Copy simple members
    N = other.N;
    V = other.V;
    num_regions = other.num_regions;
    num_districts = other.num_districts;
    num_multidistricts = other.num_multidistricts;
    map_pop = other.map_pop;
    remainder_region = other.remainder_region;

    // shallow copy
    region_ids = other.region_ids;
    region_dvals = other.region_dvals;

    // Deep copy std::vector members
    region_pops = other.region_pops;
    region_added_order = other.region_added_order;
    region_order_max = other.region_order_max;
}





// Takes a plan and reorders the regions according to the order the regions
// were split
// IT IS VERY IMPORTANT THAT THE TWO PLANS NOT POINT TO THE SAME arma::umat or everything breaks
// After reordering it copies that over to the dummy plan as well
void Plan::reorder_plan_by_oldest_split(
    Plan &dummy_plan) {
    // Make dummy plan a shallow copy of the plan
    dummy_plan.shallow_copy(*this);

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

        this->region_dvals(new_region_id) = dummy_plan.region_dvals(old_region_id);
        this->region_pops.at(new_region_id) = dummy_plan.region_pops.at(old_region_id);
        // Since the regions are in order of added reset the order added to just be 1,..., n
        this->region_added_order.at(i) = i+1;
    }

    // reset the max region counter 
    this->region_order_max = this->num_regions + 6;
    
    // copy the dummy plan over
    dummy_plan.shallow_copy(*this);
}



// Prints our object using Rcout. Should be used in Rcpp call
void Plan::Rprint() const{
    RcppThread::Rcout << "Plan with " << num_regions << " regions, " << num_districts
                      << " districts, " << num_multidistricts << " multidistricts and "
                      << arma::sum(region_dvals) << " sum of dnk and "
                      << V << " Vertices.\n";


    RcppThread::Rcout << "Region Level Values:[";
    for(int region_id = 0; region_id < num_regions; region_id++){
        RcppThread::Rcout << "(Region " << region_id <<
            ", dval=" << region_dvals(region_id) << ", pop= " <<
            region_pops.at(region_id) <<"), ";
    }
    RcppThread::Rcout << "]\n";

}


