/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2024/08
 * Purpose: Sequential Monte Carlo redistricting sampler
 ********************************************************/

#include "gredist_types.h"

/*
 * Convert zero-indxed R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const Rcpp::List &l) {
    int V = l.size();
    Graph g;
    for (int i = 0; i < V; i++) {
        g.push_back(Rcpp::as<std::vector<int>>((Rcpp::IntegerVector) l[i]));
    }
    return g;
}

/*
 * Initialize empty multigraph structure on graph with `V` vertices
 */
// TESTED
Multigraph init_multigraph(int V) {
    Multigraph g;
    for (int i = 0; i < V; i++) {
        std::vector<std::vector<int>> el;
        g.push_back(el);
    }
    return g;
}

/*
 * Make a county graph from a precinct graph and list of counties
 * County graph is list of list of 3: <cty of nbor, index of vtx, index of nbor>
 */
// TESTED
Multigraph county_graph(const Graph &g, const arma::uvec &counties) {
    int n_county = arma::max(counties);
    Multigraph cg = init_multigraph(n_county);

    int V = g.size();
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        int length = nbors.size();
        int county = counties.at(i) - 1;
        for (int j = 0; j < length; j++) {
            int nbor_cty = counties.at(nbors[j]) - 1;
            if (county == nbor_cty) continue;
            std::vector<int> el = {nbor_cty, i, nbors[j]};
            cg.at(county).push_back(el);
        }
    }

    return cg;
}


MapParams::MapParams(Rcpp::List adj_list, const arma::uvec &counties, const arma::uvec &pop,
        int ndists, double lower, double target, double upper) :
    counties(counties), pop(pop), ndists(ndists), lower(lower), target(target), upper(upper) {
        g = list_to_graph(adj_list);
        cg = county_graph(g, counties);
        V = static_cast<int>(g.size());
}




void EdgeCut::get_split_regions_info(
    int &split_region1_tree_root, int &split_region1_dval, int &split_region1_pop,
    int &split_region2_tree_root, int &split_region2_dval, int &split_region2_pop
){
    // Always make region 1 the smaller one by size (allowing for ties)

    if(cut_below_region_size <= cut_above_region_size){
        // if true then cut below is smalelr so make that region 1
        split_region1_tree_root = cut_vertex;
        split_region1_dval = cut_below_region_size;
        split_region1_pop = cut_below_pop;

        split_region2_tree_root = tree_root;
        split_region2_dval = cut_above_region_size;
        split_region2_pop = cut_above_pop;
    }else{
        // if false then cut above is smalelr so make that region 1
        split_region2_tree_root = cut_vertex;
        split_region2_dval = cut_below_region_size;
        split_region2_pop = cut_below_pop;

        split_region1_tree_root = tree_root;
        split_region1_dval = cut_above_region_size;
        split_region1_pop = cut_above_pop;
    }
};


std::array<double, 2> EdgeCut::compute_signed_pop_deviances(double target){
    // get the target populations for the regions 
    double cut_below_target = target*cut_below_region_size;
    double cut_above_target = target*cut_above_region_size;
    // get the deviation 
    double below_dev = (static_cast<double>(cut_below_pop) - cut_below_target)/cut_below_target;
    double above_dev = (static_cast<double>(cut_above_pop) - cut_above_target)/cut_above_target;
    
    std::array<double, 2> unsigned_devs = {below_dev, above_dev};

    return std::array<double, 2>{below_dev, above_dev};
}


std::array<double, 2> EdgeCut::compute_abs_pop_deviances(double target){
    // get the raw unsigned deviations
    std::array<double, 2> unsigned_devs = compute_signed_pop_deviances(target);
    // take the absolute value
    std::array<double, 2> signed_devs = {
        std::fabs(unsigned_devs.at(0)), std::fabs(unsigned_devs.at(1))
    };

    return signed_devs;
}



SplittingMethodType get_splitting_type(std::string const &splitting_type_str){
    // find the type or throw an error 
    if(splitting_type_str == "naive_top_k"){
        return SplittingMethodType::NaiveTopK;
    }else if(splitting_type_str == "uniform_valid_edge"){
        return SplittingMethodType::UnifValid;
    }else if(splitting_type_str == "expo_bigger_abs_dev"){
        return SplittingMethodType::ExpBiggerAbsDev;
    }else if(splitting_type_str == "expo_smaller_abs_dev"){
        return SplittingMethodType::ExpSmallerAbsDev;
    }else if(splitting_type_str == "experimental"){
        return SplittingMethodType::Experimental;
    }else{
        REprintf("Splitting Type %s is not a valid type!\n", 
        splitting_type_str.c_str());
        throw Rcpp::exception("Invalid splitting type passed");
    }
}

std::string splitting_method_to_str(SplittingMethodType splitting_method){
    if(splitting_method == SplittingMethodType::NaiveTopK){
        return "Naive Top K Splitter";
    }else if(splitting_method == SplittingMethodType::UnifValid){
        return "Uniform Valid Edge Splitter";
    }else if(splitting_method == SplittingMethodType::ExpBiggerAbsDev){
        return "Exponentially Weighted Absolute Bigger Deviance Splitter";
    }else if(splitting_method == SplittingMethodType::ExpSmallerAbsDev){
        return "Exponentially Weighted Absolute Smaller Deviance Splitter";
    }else if(splitting_method == SplittingMethodType::Experimental){
        return "Experimental Splitter";
    }else{
        REprintf("Splitting Type %c has no to str form!\n", 
        splitting_method);
        throw Rcpp::exception("Invalid splitting type passed to_str");
    }
}


SplittingSizeScheduleType get_splitting_size_regime(std::string const &splitting_size_regime_str){
    // find the type or throw an error 
    if(splitting_size_regime_str == "split_district_only"){
        return SplittingSizeScheduleType::DistrictOnly;
    }else if(splitting_size_regime_str == "any_valid_sizes"){
        return SplittingSizeScheduleType::AnyValidSize;
    }else if(splitting_size_regime_str == "custom"){
        return SplittingSizeScheduleType::CustomSizes;
    }else{
        REprintf("Splitting Size Regime %s is not a valid regime!\n", 
        splitting_size_regime_str.c_str());
        throw Rcpp::exception("Invalid splitting size regime passed");
    }
};



// init splitting schedule object 
SplittingSchedule::SplittingSchedule(
    const int num_splits, const int ndists, const int initial_num_regions, 
    SplittingSizeScheduleType const schedule_type,
    Rcpp::List const &control) :
    schedule_type(schedule_type), 
    ndists(ndists),
    valid_split_region_sizes_list(num_splits, std::vector<bool>(ndists+1, false)),
    all_regions_smaller_cut_sizes_to_try(ndists+1),
    all_regions_min_and_max_possible_cut_sizes(ndists+1){
    
    // If doing split district only for all sizes above 1 the only possible smaller size
    // is 1 so we can set this at the beginning and not touch it again.
    if(schedule_type == SplittingSizeScheduleType::DistrictOnly){
        for (int region_size = 2; region_size <= ndists; region_size++)
        {
            all_regions_smaller_cut_sizes_to_try[region_size] = std::vector<int> {1};
            // smallest possible cut size is 1 and biggest is region_size-1
            all_regions_min_and_max_possible_cut_sizes[region_size] = {1, region_size-1};
        }
    }else if(schedule_type == SplittingSizeScheduleType::AnyValidSize){
        // If doing any sizes then for `region_size` the possible split sizes are always 
        // if any sizes allowed then for 2, ..., ndists make it 1,...,floor(size/2)
        for (int region_size = 2; region_size <= ndists; region_size++)
        {
            // Get the lrgest smaller cut size we'll try which is the
            // floor of the size of region over 2
            int smaller_cut_size_max = std::floor(region_size / 2);
            all_regions_smaller_cut_sizes_to_try[region_size].resize(smaller_cut_size_max);
            // Now makes the entry for region size go 1,...smaller_cut_size_max
            std::iota(
                all_regions_smaller_cut_sizes_to_try[region_size].begin(), 
                all_regions_smaller_cut_sizes_to_try[region_size].end(), 
                1);
            // smallest possible cut size is 1 and biggest is region_size-1
            all_regions_min_and_max_possible_cut_sizes[region_size] = {1, region_size-1};
        }
    }else{
        throw Rcpp::exception("Not implemented yet!!");
    }

};

void SplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int num_regions
){
    // We only have to do this for custom, for any size and split district only 
    // it never changes 

    if(schedule_type != SplittingSizeScheduleType::CustomSizes){
        return;
    }

    // If split district only then make the remainder region size 1
    
}