/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2024/08
 * Purpose: Sequential Monte Carlo redistricting sampler
 ********************************************************/

#include "gredist_types.h"


// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED




// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED
Plan::Plan(
    arma::subview_col<arma::uword> region_ids_col, 
    arma::subview_col<arma::uword> region_dvals_col, 
    int N, double total_map_pop, bool split_district_only)
: region_ids(region_ids_col), region_dvals(region_dvals_col) // Initialize the reference
{
    // set number of regions, districts, multidistricts, map pop, and V
    num_regions = 1;
    num_districts = 0;
    num_multidistricts = num_regions - num_districts;
    this->V = region_ids.n_elem;
    map_pop = total_map_pop;
    this->N = N;


    // Set array which will map region id to d_nk, pop values, and str_label if tracked
    // need to assert sum is N
    // TODO add 

    region_added_order = std::vector<int>(N, -1);
    region_added_order.at(0) = 0;

    region_order_max = 0;

    region_pops = std::vector<double>(N, -1.0);
    region_pops.at(0) = total_map_pop;


    if(split_district_only){
        remainder_region = 0;
    }else{
        remainder_region = -1;
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







// Prints our object using Rcout. Should be used in Rcpp call
void Plan::Rprint() const{
    RcppThread::Rcout << "Plan with " << num_regions << " regions, " << num_districts
                      << " districts, " << num_multidistricts << " multidistricts and "
                      << sum(region_dvals) << " sum of dnk and"
                      << V << " Vertices.\n[";


    RcppThread::Rcout << "Region Level Values:";
    // for(int region_id = 0; region_id < num_regions; region_id++){
    //     RcppThread::Rcout << "( Region " << region_id <<
    //         ", dval=" << region_dvals(region_id) << ", " << region_pops.at(region_id) <<" ), ";
    // }
    for(int region_id = 0; region_id < region_dvals.n_elem; region_id++){
        RcppThread::Rcout << "( Region " << region_id <<
            ", dval=" << region_dvals(region_id) <<" ), ";
    }
    RcppThread::Rcout << "\n";
//
//     for(int i = 0; i<V; i++){
//         RcppThread::Rcout << region_labels[i] << " and " << region_num_ids[i] << '\n';
//     }
//
//     for (auto i: region_labels)
//         RcppThread::Rcout << i << ' ';
//     RcppThread::Rcout << "]\n";
}


