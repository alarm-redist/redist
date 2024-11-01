/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2024/08
 * Purpose: Sequential Monte Carlo redistricting sampler
 ********************************************************/

#include "redist_types.h"


// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED
Plan::Plan(int V, int N, double total_map_pop, bool split_district_only){
    // set number of regions, districts, multidistricts, map pop, and V
    num_regions = 1;
    num_districts = 0;
    num_multidistricts = num_regions - num_districts;
    this->V = V;
    map_pop = total_map_pop;
    this->N = N;


    // Set array which will map region id to d_nk, pop values, and str_label if tracked
    region_dvals.reserve(N); // Reserve N spots in memory
    region_dvals.push_back(N); // Make entry 0 the map population

    region_added_order.reserve(N);
    region_added_order.push_back(0);

    region_order_max = 0;

    region_pops.reserve(N);
    region_pops.push_back(total_map_pop);

    // Create array which maps vertices to their region id
    region_ids = std::vector<int>(V, 0);

    if(split_district_only){
        remainder_region = 0;
    }else{
        remainder_region = -1;
    }
    


}


// Prints our object using Rcout. Should be used in Rcpp call
void Plan::Rprint() const{
    RcppThread::Rcout << "Plan with " << num_regions << " regions, " << num_districts
                      << " districts, " << num_multidistricts << " multidistricts and "
                      << V << " Vertices.\n[";


    RcppThread::Rcout << "Region Level Values:";
    for(int region_id = 0; region_id < num_regions; region_id++){
        RcppThread::Rcout << "( Region " << region_id <<
            ", dval=" << region_dvals.at(region_id) << ", " << region_pops.at(region_id) <<" ), ";
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
