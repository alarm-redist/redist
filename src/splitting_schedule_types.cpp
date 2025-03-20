/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/01
 * Purpose: Implements Various Splitting Schedule Types
 ********************************************************/

 #include "splitting_schedule_types.h"


void SplittingSchedule::reset_splitting_and_merge_booleans(){
    // reset all regions that can be split to false 
    std::fill(valid_region_sizes_to_split.begin(), valid_region_sizes_to_split.end(), false);
    std::fill(valid_split_region_sizes.begin(), valid_split_region_sizes.end(), false);
    // reset merge pairs
    for (size_t region_size = 1; region_size <= ndists; region_size++)
    {
        std::fill(
            valid_merge_pair_sizes[region_size].begin(), 
            valid_merge_pair_sizes[region_size].end(), 
            false
        );
    }
    // reset check adj to 
    std::fill(check_adj_to_regions.begin(), check_adj_to_regions.end(), false);
}



/*
 * Constructor for district splits only splitting schedule
 */
DistrictOnlySplittingSchedule::DistrictOnlySplittingSchedule(
    const int num_splits, const int ndists
):
    SplittingSchedule(SplittingSizeScheduleType::DistrictOnly, ndists){
    // For each size above 2 make the min=1 and the max=size-1
    for (int region_size = 2; region_size <= ndists; region_size++)
    {
        all_regions_smaller_cut_sizes_to_try[region_size] = std::vector<int> {1};
        // smallest possible cut size is 1 and biggest is region_size-1
        all_regions_min_and_max_possible_cut_sizes[region_size] = {1, region_size-1};
    }
}


/*
 * Constructor for district splits only splitting schedule
 */
void DistrictOnlySplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions
){
    // The biggest size a region can be is if all but one region is a district
    int presplit_remainder_size = ndists - presplit_num_regions + 1;
    int after_split_remainder_size = presplit_remainder_size-1;

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    // we can only split the remainder region which has size presplit_remainder_size
    valid_region_sizes_to_split[presplit_remainder_size] = true;
    // the only valid split sizes are 1 and remainder_region size -1
    valid_split_region_sizes[1] = true;
    valid_split_region_sizes[after_split_remainder_size] = true;
    // we can only merge 1 and reminader region size -1 
    valid_merge_pair_sizes[1][after_split_remainder_size] = true;
    valid_merge_pair_sizes[after_split_remainder_size][1] = true;
    check_adj_to_regions[after_split_remainder_size]=true;

    return;
}


AnyRegionSplittingSchedule::AnyRegionSplittingSchedule(
    const int num_splits, const int ndists
):
    SplittingSchedule(SplittingSizeScheduleType::AnyValidSize, ndists){
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
}




void AnyRegionSplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions
){
    // The biggest size a region can be is if all but one region is a district
    int presplit_biggest_possible_size = ndists - presplit_num_regions + 1;

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    for (size_t region_size = 2; region_size <=  presplit_biggest_possible_size; region_size++)
    {
        // we can split any region between 2 and largest possible region 
        // size which is max_valid_region_size
        valid_region_sizes_to_split[region_size] = true;
        // a valid split size can be anything from 1 to largest possible region size minus 1
        valid_split_region_sizes[region_size] = true;
    }
    // make 1 a valid size
    valid_split_region_sizes[1] = true;
    // make the largest size false since you'll never get this from a split
    valid_split_region_sizes[presplit_biggest_possible_size] = false;
    check_adj_to_regions = valid_split_region_sizes;
    // valid merge pairs are anything that sums to less than or equal to the max valid size
    for (size_t region1_size = 1; region1_size < presplit_biggest_possible_size; region1_size++){
        for (size_t region2_size = 1; region2_size < presplit_biggest_possible_size; region2_size++){
            if(region1_size+region2_size <= presplit_biggest_possible_size){
                valid_merge_pair_sizes[region1_size][region2_size] = true;
            }
        }
    }
    
    return;
}



PureMSSplittingSchedule::PureMSSplittingSchedule(
    const int ndists, int const district_size_lb, int const district_size_ub
):
    SplittingSchedule(SplittingSizeScheduleType::PureMergeSplitSize, ndists){
    
    // check bounds make sense
    if(district_size_lb < 1){
        throw Rcpp::exception("Pure MS District Size Lower Bound must be at least 1!\n");
    }
    if(district_size_lb > district_size_ub){
        throw Rcpp::exception("Pure MS District Size Lower Bound must not be greater than the upper bound!\n");
    }
    if(district_size_ub > ndists){
        throw Rcpp::exception("Pure MS District Size Upper Bound must not be greater than the total number of districts in the map!\n");
    }

    for(int district_size1 = district_size_lb; 
            district_size1 <= district_size_ub;
            district_size1++){
        // mark each district size as a valid split region size
        valid_split_region_sizes[district_size1] = true;

        // don't support irregular district sizes like this right now
        if(district_size1*2 > ndists){
            throw Rcpp::exception("Pure MS not supported when one district bigger than half total seats\n");
        }

        for(int district_size2 = district_size_lb; 
            district_size2 <= district_size_ub;
            district_size2++){
            
            if(district_size1 + district_size2 <= ndists){
                valid_merge_pair_sizes[district_size1][district_size2] = true;
                valid_merge_pair_sizes[district_size2][district_size1] = true;
                valid_region_sizes_to_split[district_size1+district_size2] = true;
                // just add we'll sort and make unique later 
                all_regions_smaller_cut_sizes_to_try[district_size1+district_size2].emplace_back(
                    std::min(district_size1, district_size2)
                );
            }
        }
    }

    check_adj_to_regions = valid_split_region_sizes;
    // make smaller cut sizes unique and sort
    // get min and max sizes 

    // NOT TESTED
    for (size_t region_size = 0; region_size <= ndists; region_size++)
    {
        if(!valid_region_sizes_to_split[region_size]) continue;
        // sort and remove duplicates 
        // 1) Sort the vector
        std::sort(
            all_regions_smaller_cut_sizes_to_try[region_size].begin(), 
            all_regions_smaller_cut_sizes_to_try[region_size].end());
        
        all_regions_smaller_cut_sizes_to_try[region_size].erase( 
            std::unique(
                all_regions_smaller_cut_sizes_to_try[region_size].begin(), 
                all_regions_smaller_cut_sizes_to_try[region_size].end())
            , all_regions_smaller_cut_sizes_to_try[region_size].end() );

        int smallest_possible_size = all_regions_smaller_cut_sizes_to_try[region_size][0];

        // set biggest and smallest possible size
        all_regions_min_and_max_possible_cut_sizes[region_size] = {
            smallest_possible_size, 
            region_size-smallest_possible_size
        };
    }

}


OneCustomSplitSchedule::OneCustomSplitSchedule(
    const int num_splits, const int ndists,
    Rcpp::List const &control
): SplittingSchedule(SplittingSizeScheduleType::OneCustomSize, ndists){

    // Ensure "custom_size_split_list" exists and is a list
    if (!control.containsElementNamed("custom_size_split_list")) {
        Rcpp::stop("Error: 'permitted_split_region_sizes_list' is missing from control.");
    }

    // extract
    Rcpp::List custom_size_split_list = control["custom_size_split_list"];
    std::vector<std::vector<int>> custom_splits_vector =
        Rcpp::as<std::vector<std::vector<int>>>(custom_size_split_list);
    // ensure its the same size as the number of splits 
    if(num_splits != custom_splits_vector.size()){
        Rcpp::stop("Error: 'custom_size_split_list' must be the same length as number of splits.");
    }

    split_array_list.reserve(custom_splits_vector.size());
    // convert to vector of arrays of size 3
    for(auto split_vec: custom_splits_vector){
        // check it is length 3
        if(split_vec.size() != 3){
            Rcpp::stop("Error: All vectors in 'custom_size_split_list' must be of length 3.");
        }
        // check the first element is sum of second two
        if(split_vec[0] != split_vec[1]+split_vec[2]){
            Rcpp::stop("Error: First element in each vector in 'custom_size_split_list' must be sum of other two elements.");
        }
        // check first element is above 0
        if(split_vec[0] <= 0 || split_vec[2] <= 0 || split_vec[2] <= 0){
            Rcpp::stop("Error: All elements in each vector in 'custom_size_split_list' must be greater than zero.");
        }
        split_array_list.push_back({split_vec[0], split_vec[1], split_vec[2]});
    }
}


void OneCustomSplitSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions
){

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    // get the (presplit size, split size 1, split size 2)
    auto the_custom_split = split_array_list[split_num];

    int valid_region_size_to_split = the_custom_split[0];
    int valid_split_region1_size = the_custom_split[1];
    int valid_split_region2_size = the_custom_split[2];
    int smaller_valid_split_size = std::min(valid_split_region1_size, valid_split_region2_size);
    int bigger_valid_split_size = std::max(valid_split_region1_size, valid_split_region2_size);


    // we can only split the one valid size to split 
    valid_region_sizes_to_split[valid_region_size_to_split] = true;
    // only two valid split sizes 
    valid_split_region_sizes[valid_split_region1_size] = true;
    valid_split_region_sizes[valid_split_region2_size] = true;
    check_adj_to_regions = valid_split_region_sizes;
    // we can only merge those two valid split sizes 
    valid_merge_pair_sizes[valid_split_region1_size][valid_split_region2_size] = true;
    valid_merge_pair_sizes[valid_split_region2_size][valid_split_region1_size] = true;

    // Now set the splitting size 
    all_regions_smaller_cut_sizes_to_try[valid_region_size_to_split] = {smaller_valid_split_size};
    // set min and max 
    all_regions_min_and_max_possible_cut_sizes[valid_region_size_to_split] = {smaller_valid_split_size, bigger_valid_split_size};

    return;
}


// Returns a splitting schedule depending on the schedule type 
std::unique_ptr<SplittingSchedule> get_splitting_schedule(
    const int num_splits, const int ndists, 
    SplittingSizeScheduleType const schedule_type,
    Rcpp::List const &control
){
    if(schedule_type == SplittingSizeScheduleType::DistrictOnly){
        return std::make_unique<DistrictOnlySplittingSchedule>(num_splits, ndists);
    }else if(schedule_type == SplittingSizeScheduleType::AnyValidSize){
        return std::make_unique<AnyRegionSplittingSchedule>(num_splits, ndists);
    }else if(schedule_type == SplittingSizeScheduleType::OneCustomSize){
        return std::make_unique<OneCustomSplitSchedule>(num_splits, ndists, control);
    }else if(schedule_type == SplittingSizeScheduleType::CustomSizes){
        Rprintf("Not implemented!");
        throw Rcpp::exception("Schedule not impliemented yet!");
    }else{
        throw Rcpp::exception("Schedule not implemented yet!");
    }
}