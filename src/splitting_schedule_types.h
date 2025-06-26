#pragma once
#ifndef SPLITTING_SCHEDULE_TYPES_H
#define SPLITTING_SCHEDULE_TYPES_H


#include "redist_types.h"


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
        
        // Base constructor just creates vectors, does not set them
        SplittingSchedule(
            const SplittingSizeScheduleType schedule_type, 
            const int ndists, const int total_seats, 
            std::vector<int> const &district_seat_sizes
        ):
            schedule_type(schedule_type),
            ndists(ndists),
            total_seats(total_seats),
            smallest_district_size(*min_element(district_seat_sizes.begin(), district_seat_sizes.end())),
            largest_district_size(*max_element(district_seat_sizes.begin(), district_seat_sizes.end())),
            district_seat_sizes(district_seat_sizes),
            is_district(total_seats + 1, false),
            all_regions_smaller_cut_sizes_to_try(total_seats+1),
            all_regions_min_and_max_possible_cut_sizes(total_seats+1),
            valid_region_sizes_to_split(total_seats+1),
            valid_split_region_sizes(total_seats+1),
            check_adj_to_regions(total_seats+1),
            valid_merge_pair_sizes(total_seats+1, std::vector<bool>(total_seats+1, false)
            ){
                // make sure the district sizes are ok 
                for (auto const &a_size: district_seat_sizes){
                    if(a_size < 0) throw Rcpp::exception("District Seat Sizes must be strictly positive!\n");
                    if(a_size >= total_seats)  throw Rcpp::exception("District Seat Sizes must be less than total seats!\n");
                    // mark this as a district size 
                    is_district[a_size] = true;
                }
            };
    
        SplittingSizeScheduleType schedule_type; // the splitting type 
        int const ndists; // the number of districts 
        int const total_seats;
        double multidistrict_alpha;
        int const smallest_district_size;
        int const largest_district_size;
        std::vector<int> const &district_seat_sizes;
        std::vector<bool> is_district;

        /*
         * This is a vector of vectors where entry for index `r` it is a vector
         * of smaller cut sizes to try. ie (element, r-element)
         *      - For any split its just 1,...,floor(r/2)
         *      - for District splits its just 1
         */
        std::vector<std::vector<int>> all_regions_smaller_cut_sizes_to_try; 
    
        // For each element its an array of the minimum and maximum possible cut sizes 
        // of regions cut from 
        std::vector<std::pair<int,int>> all_regions_min_and_max_possible_cut_sizes;
    
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
        // for district splits this allows to merge and split them during merge split 
        virtual void update_cut_sizes_for_mergesplit_step(int split_num, int num_regions){
            return;
        };

        // Sets all the boolean vectors and matrices related to splitting and merging to false
        void reset_splitting_and_merge_booleans();
        
        void print_current_step_splitting_info();
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
            const int num_splits, const int ndists
        );


        void set_potential_cut_sizes_for_each_valid_size(int split_num, int presplit_num_regions) override;
        void update_cut_sizes_for_mergesplit_step(int split_num, int num_regions) override;    
};


/* 
 * Derived Class for any region split schedule (full gsmc) for single member districts
 * This is the splitting schedule for which allows for any size splits
 * to be made (where both split regions have to have at least size 1)
 */
class AnyRegionSMDSplittingSchedule : public SplittingSchedule {

    public:
        // constructor
        AnyRegionSMDSplittingSchedule(
            const int num_splits, const int ndists
        );
        
        void set_potential_cut_sizes_for_each_valid_size(
            int split_num, int presplit_num_regions) override;    
};



/* 
 * Derived Class for district only split schedule for multi-member districts
 * This is the splitting schedule where at each step we split off a district
 * from the remainder.
 */
class DistrictOnlyMMDSplittingSchedule : public SplittingSchedule {

    public:
        // constructor
        DistrictOnlyMMDSplittingSchedule(
            const int num_splits, const int ndists,
            const int total_seats,
            std::vector<int> const &district_seat_sizes
        );

        void set_potential_cut_sizes_for_each_valid_size(
            int split_num, int presplit_num_regions) override;    
        void update_cut_sizes_for_mergesplit_step(int split_num, int num_regions) override; 
};


/* 
 * Derived Class for any region split schedule for multi-member districts
 * This is the splitting schedule for which allows for any size splits
 * to be made (where both split regions have to have at least size district_size_lb)
 */
class AnyRegionMMDSplittingSchedule : public SplittingSchedule {

    public:
        // constructor
        AnyRegionMMDSplittingSchedule(
            const int num_splits, const int ndists,
            const int total_seats,
            std::vector<int> const &district_seat_sizes
        );

        std::vector<bool> reachable_size;
        
        void set_potential_cut_sizes_for_each_valid_size(
            int split_num, int presplit_num_regions) override;    
};


/* 
 * Derived Class for pure merge split. So you can merge any valid district
 * sizes and create any district sizes. 
 */
class PureMSSplittingSchedule : public SplittingSchedule {

    public:
        // constructor
        PureMSSplittingSchedule(
            const int ndists, const int total_size,
            std::vector<int> const &district_seat_sizes
        );

        void set_potential_cut_sizes_for_each_valid_size(
            int split_num, int presplit_num_regions) override{
                throw Rcpp::exception("Dont call this method for pure MS!");
            };    
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
            const int num_splits, 
            const int ndists,  
            Rcpp::List const &control);
        
        void set_potential_cut_sizes_for_each_valid_size(int split_num, int presplit_num_regions) override;
};


    
std::unique_ptr<SplittingSchedule> get_splitting_schedule(
    const int num_splits, const int ndists, const int total_seats,
    std::vector<int> const &district_seat_sizes,
    SplittingSizeScheduleType const schedule_type,
    Rcpp::List const &control
);



#endif
