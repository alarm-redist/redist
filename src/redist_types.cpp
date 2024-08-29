/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2024/08
 * Purpose: Sequential Monte Carlo redistricting sampler
 ********************************************************/

#include "redist_types.h"


// Define the constructor template outside the class
// THIS ONLY CONSTRUCTS A ONE REGION MAP. ANYTHING ELSE MUST BE UPDATED
Plan::Plan(int V, int N, double total_map_pop){
    // Check number of regions and districts make sense
    /*
    if(num_regions < 1){
        RcppThread::Rcerr << "ERROR: Invalid number of regions!\n";
    }
    if(num_districts > num_regions){
        RcppThread::Rcout << "ERROR: Invalid number of districts!\n";
    }
     */

    // set number of regions, districts, multidistricts, map pop, and V
    num_regions = 1;
    num_districts = 0;
    num_multidistricts = num_regions - num_districts;
    this->V = V;
    map_pop = total_map_pop;
    this->N = N;




    // Set label to dnk and pop value map
    r_to_d_map = std::map<std::string, int> {{"R1", N}};
    r_to_pop_map = std::map<std::string, double> {{"R1", map_pop}};

    // Create region labels, pop, and dnk
    // THIS ONLY WORKS FOR A BLANK REGION RIGHT NOW
    region_labels = std::vector<std::string>(V, "R1");
    region_num_ids = std::vector<int>(V, 0);
    region_dval = std::vector<int>(V, r_to_d_map["R1"]);
    region_pop = std::vector<double>(V, r_to_pop_map["R1"]);

    // Set maps allowing going between integer and string region ids
    str_label_to_num_id_map = std::map<std::string, int> {{"R1", 0}};; // Maps region label values to integer id number
    num_id_to_str_label_map = std::map<int, std::string> {{0, "R1"}};; // Maps integer id number to region label values

}


// Prints our object using Rcout. Should be used in Rcpp call
void Plan::Rprint() const{
    RcppThread::Rcout << "Plan with " << num_regions << " regions, " << num_districts
                      << " districts and "<< V << " Vertices.\n[";


    RcppThread::Rcout << "Region Level Values:";
    for (const auto &[key, value]: r_to_d_map ) {
        RcppThread::Rcout << "(" << str_label_to_num_id_map.at(key) << " aka " << key << " also " <<
            num_id_to_str_label_map.at(str_label_to_num_id_map.at(key)) << ", " << value << ", " << r_to_pop_map.at(key) <<" ), ";
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
