#pragma once
#ifndef SPLITTING_SCHEDULE_TYPES_H
#define SPLITTING_SCHEDULE_TYPES_H


#include "gredist_types.h"


/*
 * Abstract class for object used to manage what sizes are allowed to be split
 * at each step. 
 * 
 * It is used in all three splitting, merge-split, and weight computation steps.
 */
class SplittingSchedule {

    private:
        // this are both only actually tracked for custom splitting schedules 
        std::vector<std::vector<bool>> valid_split_region_sizes_list; // the allowable splitting sizes for each split
        std::vector<std::vector<bool>> valid_presplit_region_sizes_list; // the allowable presplit region sizes for each split
        
    public:
        virtual ~SplittingSchedule() = default;
        // constructor 
        SplittingSchedule(
            const int num_splits, const int ndists, const int initial_num_regions, 
            SplittingSizeScheduleType const schedule_type,
            Rcpp::List const &control);
        
        // Base constructor just creates vectors, does not set them
        SplittingSchedule(
            const SplittingSizeScheduleType schedule_type, const int num_splits, const int ndists, const int initial_num_regions
        ):
            schedule_type(schedule_type),
            ndists(ndists),
            all_regions_smaller_cut_sizes_to_try(ndists+1),
            all_regions_min_and_max_possible_cut_sizes(ndists+1),
            valid_region_sizes_to_split(ndists+1),
            valid_split_region_sizes(ndists+1),
            check_adj_to_regions(ndists+1),
            valid_merge_pair_sizes(ndists+1, std::vector<bool>(ndists+1, false)
            ) {};
    
        SplittingSizeScheduleType schedule_type; // the splitting type 
        int ndists; // the number of districts 
    

        /*
         * This is a vector of vectors where entry for index `r` it is a vector
         * of smaller cut sizes to try. ie (element, r-element)
         *      - For any split its just 1,...,floor(r/2)
         *      - for District splits its just 1
         */
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

        // 1-indexed 2D matrix for if two regions (size1, size2) can be merged 
        std::vector<std::vector<bool>> valid_merge_pair_sizes;
    
        // This sets it for that split
        virtual void set_potential_cut_sizes_for_each_valid_size(int split_num, int presplit_num_regions) = 0;
        // Sets all the boolean vectors and matrices related to splitting and merging to false
        void reset_splitting_and_merge_booleans();    
};
    

/* 
 * Derived Class for one district split schedule
 * 
 * This is the splitting schedule for splitting 1 district at a time
 */
class DistrictOnlySplittingSchedule : public SplittingSchedule {

    public:
        // constructor
        DistrictOnlySplittingSchedule(
            const int num_splits, const int ndists, const int initial_num_regions
        );


        void set_potential_cut_sizes_for_each_valid_size(int split_num, int presplit_num_regions) override;    
};


/* 
 * Derived Class for any region split schedule (full gsmc)
 * This is the splitting schedule for which allows for any size splits
 * to be made (where both split regions have to have at least size 1)
 */
class AnyRegionSplittingSchedule : public SplittingSchedule {

    public:
        // constructor
        AnyRegionSplittingSchedule(
            const int num_splits, const int ndists, const int initial_num_regions
        );

        void set_potential_cut_sizes_for_each_valid_size(int split_num, int presplit_num_regions) override;    
};



/* 
 * Derived Class for one custom split schedule. This is where only one type
 * of split is allowed in each step 
 * 
 * BE VERY CAREFUL WITH THIS! Not all configurations are guaranteed to be correct
 */
class OneCustomSplitSchedule : public SplittingSchedule {

    private:
        std::vector<std::array<int, 3>> split_array_list; // list of the split for each round
        // array is (presplit size, split size 1, split size 2)

    public:
        // constructor
        OneCustomSplitSchedule(
            const int num_splits, const int ndists, const int initial_num_regions, 
            Rcpp::List const &control);
        
        void set_potential_cut_sizes_for_each_valid_size(int split_num, int presplit_num_regions) override;
};


    
std::unique_ptr<SplittingSchedule> get_splitting_schedule(
    const int num_splits, const int ndists, const int initial_num_regions, 
    SplittingSizeScheduleType const schedule_type,
    Rcpp::List const &control
);



#endif