/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/01
 * Purpose: Implements Various Splitting Schedule Types
 ********************************************************/

 #include "splitting_schedule_types.h"


// init splitting schedule object 
SplittingSchedule::SplittingSchedule(
    const int num_splits, const int ndists, const int initial_num_regions, 
    SplittingSizeScheduleType const schedule_type,
    Rcpp::List const &control) :
    schedule_type(schedule_type), 
    ndists(ndists),
    all_regions_smaller_cut_sizes_to_try(ndists+1),
    all_regions_min_and_max_possible_cut_sizes(ndists+1),
    valid_region_sizes_to_split(ndists+1),
    valid_split_region_sizes(ndists+1),
    valid_merge_pair_sizes(ndists+1, std::vector<bool>(ndists+1, false))
    {
    
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
    }else if(schedule_type == SplittingSizeScheduleType::CustomSizes){
        // Get the permitted split sizes 
        Rcpp::List permitted_split_region_sizes_list = control["permitted_split_region_sizes_list"];
        std::vector<std::vector<int>> permitted_split_region_sizes =
            Rcpp::as<std::vector<std::vector<int>>>(permitted_split_region_sizes_list);
        // make sure the length is the same number of splits 
        if(permitted_split_region_sizes.size() != num_splits){
            throw Rcpp::exception("Number of splits not equal to permitted split region sizes list!\n");
        }
        // Get the permitted split sizes 
        Rcpp::List permitted_presplit_region_sizes_list = control["permitted_presplit_region_sizes_list"];
        std::vector<std::vector<int>> permitted_presplit_region_sizes =
            Rcpp::as<std::vector<std::vector<int>>>(permitted_presplit_region_sizes_list);
        // make sure the length is the same number of splits 
        if(permitted_presplit_region_sizes.size() != num_splits){
            throw Rcpp::exception("Number of splits not equal to permitted split region sizes list!\n");
        }


        // set this
        valid_split_region_sizes_list = std::vector<std::vector<bool>>(num_splits, std::vector<bool>(ndists+1, false));
        valid_presplit_region_sizes_list = std::vector<std::vector<bool>>(num_splits, std::vector<bool>(ndists+1, false));
        // for each split we go through and set things
        for (size_t split_num = 0; split_num < num_splits; split_num++)
        {   
            // sanity check that this is always at least one
            int num_valid_split_regions = 0;
            // we go through each valid split size and mark as true 
            for(auto const valid_split_region_size:permitted_split_region_sizes[split_num]){
                valid_split_region_sizes_list[split_num][valid_split_region_size] = true;
                ++num_valid_split_regions;
            }
            if(num_valid_split_regions < 1) throw Rcpp::exception("NO bad custom size!\n");
            // sanity check that this is always at least one
            int num_valid_presplit_regions = 0;
            // we go through each valid split size and mark as true 
            for(auto const valid_presplit_region_size:permitted_presplit_region_sizes[split_num]){
                valid_presplit_region_sizes_list[split_num][valid_presplit_region_size] = true;
                ++num_valid_presplit_regions;
            }
            if(num_valid_presplit_regions < 1) throw Rcpp::exception("No bad custom size!\n");
        }
    }else{
        throw Rcpp::exception("Not implemented yet!!");
    }

};



std::tuple<std::vector<bool>, std::vector<std::vector<int>>, std::vector<std::array<int, 2>>> get_all_valid_split_NAME_NEEDED_STUFF(
    const std::vector<bool> &valid_split_region_sizes, const std::vector<bool> &valid_presplit_region_sizes
){
    int ndists = static_cast<int>(valid_split_region_sizes.size())-1;
    int num_valid_split_regions = std::accumulate(valid_split_region_sizes.begin(), valid_split_region_sizes.end(), 0);
    // make the ouput vector 
    std::vector<std::vector<int>> all_regions_smaller_cut_sizes_to_try(
        ndists+1, std::vector<int>(0)
    );
    std::vector<bool> valid_region_sizes_to_split(ndists+1, false);


    // make a vector of sets to avoid possible duplicates
    // I don't think duplicates are possible but too lazy to prove
    std::vector<std::set<int>> all_regions_smaller_cut_sizes_to_try_sets(
            ndists+1
        );

    // get the actual sizes and put them in a vector to iterate over later
    std::vector<int> valid_split_region_sizes_vals;
    valid_split_region_sizes_vals.reserve(num_valid_split_regions);
    for (int split_region_size = 1; split_region_size <= ndists; split_region_size++)
    {
        if(valid_split_region_sizes[split_region_size]){
            valid_split_region_sizes_vals.push_back(split_region_size);
            int doubled_region_size = 2*split_region_size;
            // see if we can split into two pieces of this size
            if(doubled_region_size <= ndists && valid_presplit_region_sizes[doubled_region_size]){
                // Then we mark the sum of the regions as valid 
                valid_region_sizes_to_split[2*split_region_size] = true;
                // add the smaller of the two to the set of valid cut sizes 
                all_regions_smaller_cut_sizes_to_try_sets[2*split_region_size].insert(
                    split_region_size
                );
            }
        }
    }


    // Now iterate over all possible pairs of different valid sizes
    for (size_t i = 0; i < num_valid_split_regions; i++){
        for (size_t j = 0; j < i; j++)
        {
            int region1_size = valid_split_region_sizes_vals[i];
            int region2_size = valid_split_region_sizes_vals[j];
            int combined_region_size = region1_size+region2_size;

            // check if the merged region has both a valid size and its
            // a valid presplit region size 
            if(combined_region_size <= ndists && valid_presplit_region_sizes[combined_region_size]){
                // Then we mark the sum of the regions as valid 
                valid_region_sizes_to_split[region1_size + region2_size] = true;
                // add the smaller of the two to the set of valid cut sizes 
                all_regions_smaller_cut_sizes_to_try_sets[region1_size + region2_size].insert(
                    std::min(region1_size, region2_size)
                );
            }
        }
    }

    std::vector<std::array<int, 2>> all_regions_min_and_max_possible_cut_sizes(ndists+1);
    
    // now we convert the sets to vectors 
    for (int region_size = 1; region_size <= ndists; region_size++)
    {
    all_regions_smaller_cut_sizes_to_try[region_size] = std::vector<int>(
            all_regions_smaller_cut_sizes_to_try_sets[region_size].begin(),
            all_regions_smaller_cut_sizes_to_try_sets[region_size].end()
        );
        // set this too
        if(valid_region_sizes_to_split[region_size]){
            int smallest_possible_split_region_size = all_regions_smaller_cut_sizes_to_try[region_size][0];
            int biggest_possible_split_region_size = region_size-smallest_possible_split_region_size;
            all_regions_min_and_max_possible_cut_sizes[region_size][0] = smallest_possible_split_region_size;
            all_regions_min_and_max_possible_cut_sizes[region_size][1] = biggest_possible_split_region_size;
        }
    
    }

    
    return std::make_tuple(valid_region_sizes_to_split, all_regions_smaller_cut_sizes_to_try, all_regions_min_and_max_possible_cut_sizes);

}

void SplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int num_regions
){
    // The biggest size a region can be is if all but one region is a district
    int max_valid_region_size = ndists - num_regions + 1;
    // for custom one we can just use results from initialization
    if(schedule_type == SplittingSizeScheduleType::CustomSizes){
        // Use previous result for this step 
        valid_split_region_sizes = valid_split_region_sizes_list[split_num];
        // helper function 
        auto helper_output = get_all_valid_split_NAME_NEEDED_STUFF(
            valid_split_region_sizes, 
            valid_presplit_region_sizes_list[split_num]
        );

        valid_region_sizes_to_split = std::get<0>(helper_output); 
        all_regions_smaller_cut_sizes_to_try = std::get<1>(helper_output);
        all_regions_min_and_max_possible_cut_sizes = std::get<2>(helper_output);
        return;
    }

    // for split district only and any size we can compute from scratch easily 
    // reset all regions that can be split to false 
    std::fill(valid_region_sizes_to_split.begin(), valid_region_sizes_to_split.end(), false);
    std::fill(valid_split_region_sizes.begin(), valid_split_region_sizes.end(), false);
    
    if(schedule_type == SplittingSizeScheduleType::DistrictOnly){
        // we can only split the remainder region which has size max_valid_region_size
        valid_region_sizes_to_split[max_valid_region_size] = true;
        // the only valid split sizes are 1 and remainder_region size -1
        valid_split_region_sizes[1] = true;
        valid_split_region_sizes[ndists-num_regions] = true;
    }else if(schedule_type == SplittingSizeScheduleType::AnyValidSize){
        for (size_t region_size = 2; region_size <=  max_valid_region_size; region_size++)
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
        valid_split_region_sizes[max_valid_region_size] = false;
    }else{
        throw Rcpp::exception("Not implemented yet!!");
    }

}