#pragma once
#ifndef SPLITTING_SCHEDULE_TYPES_H
#define SPLITTING_SCHEDULE_TYPES_H


#include "gredist_types.h"



// base class used to manage what sizes we allow for splits at each step 
// used for splitting and in computing the weights 
class SplittingSchedule {

    private:
        // this are both only actually tracked for custom splitting schedules 
        std::vector<std::vector<bool>> valid_split_region_sizes_list; // the allowable splitting sizes for each split
        std::vector<std::vector<bool>> valid_presplit_region_sizes_list; // the allowable presplit region sizes for each split
        
    public:
        // constructor 
        SplittingSchedule(
            const int num_splits, const int ndists, const int initial_num_regions, 
            SplittingSizeScheduleType const schedule_type,
            Rcpp::List const &control);
        
        SplittingSchedule(
            const int num_splits, const int ndists, const int initial_num_regions, 
            Rcpp::List const &control);
    
        SplittingSizeScheduleType schedule_type; // the splitting type 
        int ndists; // the number of districts 
    
        // This is a vector of vectors where entry 
        // for index `r` it is a vector of smaller cut sizes to try
        // ie (element, r-element)
        // For any split its just 1,...,floor(r/2)
        // for District splits its just 1
        std::vector<std::vector<int>> all_regions_smaller_cut_sizes_to_try; 
    
        // For each element its an array of the minimum and maximum possible cut sizes 
        // of regions cut from 
        std::vector<std::array<int, 2>> all_regions_min_and_max_possible_cut_sizes;
    
        // 1-indexed vector of if the region can be split
        std::vector<bool> valid_region_sizes_to_split;
        // 1-indexed vector of if that region size could have been the result of a valid split
        std::vector<bool> valid_split_region_sizes;
        // 1-indexed vector of if we should check adjacency of vertices in a region of that size 
        std::vector<bool>  check_adj_to_regions;

        // 1-indexed matrix for if two regions (size1, size2) can be merged 
        std::vector<std::vector<bool>> valid_merge_pair_sizes;
    
        // This sets it for that split
        void set_potential_cut_sizes_for_each_valid_size(int split_num, int num_regions);
        
};
    

/* 
 * Derived Class for one district split 
 */
class DistrictOnlySplittingSchedule : public SplittingSchedule {

    public:
        // constructor
        DistrictOnlySplittingSchedule(
            const int num_splits, const int ndists, const int initial_num_regions, 
            SplittingSizeScheduleType const schedule_type,
            Rcpp::List const &control);
        
        // void set_potential_cut_sizes_for_each_valid_size(int split_num, int num_regions) override;
    
};



/* 
 * Derived Class for one custom split schedule. This is where only one type
 * of split it allowed in each step 
 * 
 */
class OneCustomSplitSchedule : public SplittingSchedule {

    public:
        // constructor
        OneCustomSplitSchedule(
            const int num_splits, const int ndists, const int initial_num_regions, 
            SplittingSizeScheduleType const schedule_type,
            Rcpp::List const &control);
        
        // void set_potential_cut_sizes_for_each_valid_size(int split_num, int num_regions) override;
    
};




std::tuple<
    std::vector<bool>, 
    std::vector<std::vector<int>>, 
    std::vector<std::array<int, 2>>
> get_all_valid_split_NAME_NEEDED_STUFF(
    const std::vector<bool> &valid_split_region_sizes, const std::vector<bool> &valid_presplit_region_sizes
);
    
    


#endif