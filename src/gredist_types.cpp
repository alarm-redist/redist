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


// MapParams::MapParams(Rcpp::List adj_list, const arma::uvec &counties, const arma::uvec &pop,
//         int ndists, double lower, double target, double upper) :
//      g(list_to_graph(adj_list)), counties(counties), cg(county_graph(g, counties)), pop(pop),
//      V(static_cast<int>(g.size())), ndists(ndists), lower(lower), target(target), upper(upper) {
//         // g = list_to_graph(adj_list);
//         // cg = county_graph(g, counties);
//         // V = static_cast<int>(g.size());
// }




void EdgeCut::get_split_regions_info(
    int &split_region1_tree_root, int &split_region1_dval, int &split_region1_pop,
    int &split_region2_tree_root, int &split_region2_dval, int &split_region2_pop
) const{
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


std::array<double, 2> EdgeCut::compute_signed_pop_deviances(double target) const{
    // get the target populations for the regions 
    double cut_below_target = target*cut_below_region_size;
    double cut_above_target = target*cut_above_region_size;
    // get the deviation 
    double below_dev = (static_cast<double>(cut_below_pop) - cut_below_target)/cut_below_target;
    double above_dev = (static_cast<double>(cut_above_pop) - cut_above_target)/cut_above_target;
    
    return std::array<double, 2>{below_dev, above_dev};
}


std::array<double, 2> EdgeCut::compute_abs_pop_deviances(double target) const{
    // get the raw unsigned deviations
    std::array<double, 2> unsigned_devs = compute_signed_pop_deviances(target);
    // take the absolute value
    std::array<double, 2> signed_devs = {
        std::fabs(unsigned_devs.at(0)), std::fabs(unsigned_devs.at(1))
    };

    return signed_devs;
}


// loads a sampling spaces type enum from a control string
SamplingSpace get_sampling_space(std::string const &sampling_space_str){
    // find the type or throw an error 
    if(sampling_space_str == "graph_plan_space"){
        return SamplingSpace::GraphSpace;
    }else if(sampling_space_str == "spanning_forest_space"){
        return SamplingSpace::ForestSpace;
    }else{
        REprintf("Splitting Type %s is not a valid sampling space!\n", 
            sampling_space_str.c_str());
        throw Rcpp::exception("Invalid sampling space passed");
    }
}

// Get convinient string representation
std::string sampling_space_to_str(SamplingSpace sampling_space){
    if(sampling_space == SamplingSpace::GraphSpace){
        return "Graph";
    }else if(sampling_space == SamplingSpace::ForestSpace){
        return "Forest";
    }else{
        REprintf("Sampling Space Type ?? has no to str form!\n");
        throw Rcpp::exception("Invalid splitting type passed to_str");
    }
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
        REprintf("Splitting Type ?? has no to str form!\n");
        throw Rcpp::exception("Invalid splitting type passed to_str");
    }
}


SplittingSizeScheduleType get_splitting_size_regime(std::string const &splitting_size_regime_str){
    // find the type or throw an error 
    if(splitting_size_regime_str == "split_district_only"){
        return SplittingSizeScheduleType::DistrictOnly;
    }else if(splitting_size_regime_str == "any_valid_sizes"){
        return SplittingSizeScheduleType::AnyValidSize;
    }else if(splitting_size_regime_str == "one_custom_size"){
        return SplittingSizeScheduleType::OneCustomSize;
    }else if(splitting_size_regime_str == "custom"){
        return SplittingSizeScheduleType::CustomSizes;
    }else{
        REprintf("Splitting Size Regime %s is not a valid regime!\n", 
        splitting_size_regime_str.c_str());
        throw Rcpp::exception("Invalid splitting size regime passed");
    }
};


