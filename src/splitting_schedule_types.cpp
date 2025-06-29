/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/01
 * Purpose: Implements Various Splitting Schedule Types
 ********************************************************/

 #include "splitting_schedule_types.h"

constexpr bool DEBUG_SPLITTING_SCHEDULES_VERBOSE = false;

void SplittingSchedule::reset_splitting_and_merge_booleans(){
    // reset all regions that can be split to false 
    std::fill(valid_region_sizes_to_split.begin(), valid_region_sizes_to_split.end(), false);
    std::fill(valid_split_region_sizes.begin(), valid_split_region_sizes.end(), false);
    // reset merge pairs and split sizes to try
    for (size_t region_size = 1; region_size <= total_seats; region_size++)
    {
        std::fill(
            valid_merge_pair_sizes[region_size].begin(), 
            valid_merge_pair_sizes[region_size].end(), 
            false
        );
        // clear split sizes to try 
        all_regions_smaller_cut_sizes_to_try[region_size].clear();

    }
    // reset check adj to 
    std::fill(check_adj_to_regions.begin(), check_adj_to_regions.end(), false);
}

void SplittingSchedule::print_current_step_splitting_info(){
    Rprintf("The following region sizes can split:\n");
    for (size_t i = 1; i <= total_seats; i++)
    {                
        if(valid_region_sizes_to_split[i]){
        Rprintf("\tRegion Size %d | ", (int) i);
        Rprintf("min/max (%d, %d)", 
            all_regions_min_and_max_possible_cut_sizes[i].first,
            all_regions_min_and_max_possible_cut_sizes[i].second);
        Rprintf(" | possible split sizes: ");
        for(auto smaller_size: all_regions_smaller_cut_sizes_to_try[i]){
            Rprintf("(%d, %d) ", (int) smaller_size, static_cast<int>(i-smaller_size));
        }
        Rprintf("\n");
        }
        
    }
    Rprintf("All Possible Split Sizes are:(");
    for (int i = 1; i <= total_seats; i++)
    {
        if(valid_split_region_sizes[i]){
            Rprintf("%d, ", i);
        }
    }
    Rprintf(")\n");
    // Rprintf("Now doing check_adj_to:\n");
    // for (size_t i = 1; i <= total_seats; i++)
    // {
    //     std::string t_str = check_adj_to_regions[i] ?  "true" : "false";
    //     Rprintf("%d: %s\n", (int) i, t_str.c_str());
    // }

    // Rprintf("The following valid merge pairs are:\n");
    // std::cout << "_ |";
    // for (int i = 1; i <= total_seats; i++)
    // {
    //     if(!valid_region_sizes_to_split[i] && !valid_split_region_sizes[i]) continue;
    //     // std::cout << " " << i;
    //     std::cout << " " << "\033[4m" << i << "\033[0m";
    // }
    // std::cout << std::endl;
    // for (int i = 1; i <= total_seats; i++)
    // {
    //     if(!valid_region_sizes_to_split[i] && !valid_split_region_sizes[i]) continue;
    //     std::cout << i << " | ";
    //     for (int j = 1; j <= total_seats; j++)
    //     {
    //         if(!valid_region_sizes_to_split[j] && !valid_split_region_sizes[j]) continue;
    //         std::cout << (valid_merge_pair_sizes[i][j] ? "1 " : "0 ");
    //     }
    //     std::cout << std::endl;
    // }
    
    

    // int row_num = 0;
    // for (const auto& row : valid_merge_pair_sizes) {  // Iterate over rows
    //     if(row_num == 0){
    //         row_num++;
    //         continue;
    //     }
    //     std::cout << row_num++ << " | ";
    //     int col_num = 0;
    //     for (bool val : row) {        // Iterate over elements in the row
    //         if(col_num == 0){
    //             col_num++;
    //             continue;
    //         }
    //         std::cout << (val ? "1 " : "0 ");  // Print 1 for true, 0 for false
    //     }
    //     std::cout << "\n";  // Newline after each row
    // }
}



/*
 * Constructor for district splits only splitting schedule
 */
DistrictOnlySplittingSchedule::DistrictOnlySplittingSchedule(
    const int num_splits, const int ndists
):
    SplittingSchedule(SplittingSizeScheduleType::DistrictOnlySMD, ndists, ndists, std::vector<int> {1}){
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

// makes it so districts can be merged during merge split
void DistrictOnlySplittingSchedule::update_cut_sizes_for_mergesplit_step(int split_num, int num_regions){
    // if only one region then do nothing because only 1 district (if ndists=2 then we don't need to worry)
    if(num_regions == 2) return;
    // make it so we can merge districts 
    valid_merge_pair_sizes[1][1] = true;
    // make it so we can split a region of size 2
    valid_region_sizes_to_split[2] = true;
    return;
}


AnyRegionSMDSplittingSchedule::AnyRegionSMDSplittingSchedule(
    const int num_splits, const int ndists
):
    SplittingSchedule(SplittingSizeScheduleType::AnyValidSizeSMD, ndists, ndists, std::vector<int> {1}){
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


void AnyRegionSMDSplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
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


DistrictOnlyMMDSplittingSchedule::DistrictOnlyMMDSplittingSchedule(
            const int num_splits, const int ndists,
            const int total_seats, std::vector<int> const &district_seat_sizes
):
    SplittingSchedule(SplittingSizeScheduleType::DistrictOnlyMMD, ndists, total_seats, district_seat_sizes){
    // make sure the district sizes are ok 
    for (auto const &a_size: district_seat_sizes){
        if(a_size < 0) throw Rcpp::exception("District Seat Sizes must be strictly positive!\n");
        if(a_size >= total_seats)  throw Rcpp::exception("District Seat Sizes must be less than total seats!\n");
        // mark this as a district size 
        is_district[a_size] = true;
    }

    // right now only ranges are supported 
    int cur_index = 0;
    for (int i = smallest_district_size; i < largest_district_size; i++)
    {
        if(district_seat_sizes[cur_index] != i){
            throw Rcpp::exception("For MMD only ranges of sizes are supported!\n");
        }
        cur_index++;
    }
    if(largest_district_size - smallest_district_size + 1 != district_seat_sizes.size()){
        throw Rcpp::exception("For MMD only ranges of sizes are supported!\n");
    }
    
    // For a range it must be smallest_district_size <= total_seats/ndists <= largest_district_size
    if(smallest_district_size > total_seats / static_cast<double>(ndists)){
        throw Rcpp::exception("It is not possible to split this map with the given district sizes! The district sizes are too large\n");
    }
    if(largest_district_size < total_seats / static_cast<double>(ndists)){
        throw Rcpp::exception("It is not possible to split this map with the given district sizes! There's not enough regions\n");
    }
}

void DistrictOnlyMMDSplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions
){
    // get current number of districts before splitting
    int presplit_ndists = presplit_num_regions - 1;
    // get the number of districts the remainder will need to be split into 
    // which is the number of districts minus the number of districts currently
    int remainder_ndists = ndists - presplit_ndists;
    // Define largest and smallest possible remainder sizes

    // the biggest remainder size assumes we split smallest district each time 
    int presplit_biggest_possible_size = total_seats - presplit_ndists * smallest_district_size;

    if(DEBUG_SPLITTING_SCHEDULES_VERBOSE){
    REprintf("%d Remaining Districts!\n", remainder_ndists);
    REprintf("Total Seats: %d\n Presplit Districts: %d\n Smallest District Size: %d\n Largest District Size: %d\n", 
        total_seats, presplit_ndists, smallest_district_size, largest_district_size);
    REprintf("There are %d presplit regions! Starting with Biggest presplit size %d\n", 
        presplit_num_regions, presplit_biggest_possible_size);
    }

    // If this isn't a splittable number then keep increasing by one until we get something 
    while(
        smallest_district_size > presplit_biggest_possible_size / static_cast<double>(remainder_ndists) ||
        largest_district_size < presplit_biggest_possible_size / static_cast<double>(remainder_ndists)
    ){
        if(DEBUG_SPLITTING_SCHEDULES_VERBOSE) REprintf("%d\n", presplit_biggest_possible_size);
        presplit_biggest_possible_size--;
        if(presplit_biggest_possible_size <= 0 || presplit_biggest_possible_size >= total_seats){
            REprintf("WE'RE BREAKING FREE!\n");
            throw Rcpp::exception("We got MMD remainder size issue with the presplit biggest possible size!\n");
        }
     }
    // REprintf("Biggest Remainder Size: %d\n", presplit_biggest_possible_size);


     // start off with 
    int presplit_smallest_possible_size = std::max(0, total_seats - presplit_ndists * largest_district_size);
    if(DEBUG_SPLITTING_SCHEDULES_VERBOSE){
    REprintf("There are %d presplit regions! Starting with Smallest presplit size %d\n", 
        presplit_num_regions, presplit_smallest_possible_size);
    }

    while(
        smallest_district_size > presplit_smallest_possible_size / static_cast<double>(remainder_ndists) ||
        largest_district_size < presplit_smallest_possible_size / static_cast<double>(remainder_ndists)
    ){
        if(DEBUG_SPLITTING_SCHEDULES_VERBOSE) REprintf("%d\n", presplit_smallest_possible_size);
        presplit_smallest_possible_size++;
        if(presplit_smallest_possible_size >= total_seats){
            REprintf("WE'RE BREAKING FREE!\n");
            throw Rcpp::exception("We got MMD remainder size issue!\n");
        }
    }
    // REprintf("Smallest Remainder Size: %d\n", presplit_smallest_possible_size);

    // Don't think this should be possible but just adding a flag 
    if(presplit_biggest_possible_size <= largest_district_size){
        REprintf("BIG PROBLEM: Biggest remainder size is %d but that is less than largest district %d!\n", 
        presplit_biggest_possible_size, largest_district_size);
        throw Rcpp::exception("Error in MMD Split district only!!\n");
    }

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    if(DEBUG_SPLITTING_SCHEDULES_VERBOSE){
    REprintf("We're looking at remainders between %d and %d\n",
    presplit_biggest_possible_size ,presplit_smallest_possible_size);
    }
    // now look at each possible remainder size
    for (int remainder_size = presplit_smallest_possible_size; 
         remainder_size <= presplit_biggest_possible_size; 
         remainder_size++
    ){
        // check if it is possible to reach 
        // if not smallest_district_size <= remainder_size/remainder_ndists <= largest_district_size
        if(
            smallest_district_size > remainder_size/static_cast<double>(remainder_ndists) ||
            largest_district_size < remainder_size/static_cast<double>(remainder_ndists)
        ){
            continue;
        }
        // Now check for each district size if it can be split off 
        for (int district_size = smallest_district_size; 
             district_size <= largest_district_size; 
             district_size++)
        {
            if(district_size >= remainder_size) continue;
            int new_remainder_size = remainder_size - district_size;

            double new_piece_ratio = (new_remainder_size) /static_cast<double>(remainder_ndists - 1);
            // if smallest_district_size <= (remainder_size-district_size)/(remainder_ndists-1) <= largest_district_size then it is ok
            if(smallest_district_size <= new_piece_ratio && new_piece_ratio <= largest_district_size){
                // If so then add this size as an option 
                // Add the smaller of the two cut pieces
                all_regions_smaller_cut_sizes_to_try[remainder_size].push_back(
                    std::min(new_remainder_size, district_size)
                );
                // Now also add this district as a possible split size 
                valid_split_region_sizes[district_size] = true;
                valid_split_region_sizes[new_remainder_size] = true;
                // The district and the new remainder are also valid merge pairs 
                valid_merge_pair_sizes[new_remainder_size][district_size] = true;
                valid_merge_pair_sizes[district_size][new_remainder_size] = true;
                // We also add these as check adjacent to just because things aren't as 
                // easy as SMD case 
                check_adj_to_regions[new_remainder_size] = true;
                check_adj_to_regions[district_size] = true;
            }
        }
        
        if(all_regions_smaller_cut_sizes_to_try[remainder_size].size() == 0) continue;
        valid_region_sizes_to_split[remainder_size] = true;
        // Now sort the vector and make it unique 
        std::sort(
            all_regions_smaller_cut_sizes_to_try[remainder_size].begin(), 
            all_regions_smaller_cut_sizes_to_try[remainder_size].end()
        );
        auto last = std::unique(
            all_regions_smaller_cut_sizes_to_try[remainder_size].begin(), 
            all_regions_smaller_cut_sizes_to_try[remainder_size].end()
        );
        all_regions_smaller_cut_sizes_to_try[remainder_size].erase(last, all_regions_smaller_cut_sizes_to_try[remainder_size].end());
        // now find the min and max values 
        all_regions_min_and_max_possible_cut_sizes[remainder_size] = {
            all_regions_smaller_cut_sizes_to_try[remainder_size][0],
            remainder_size - all_regions_smaller_cut_sizes_to_try[remainder_size][0]
        };

    }

    std::fill(
        check_adj_to_regions.begin(),
        check_adj_to_regions.end(),
        true
    );
    
}


void DistrictOnlyMMDSplittingSchedule::update_cut_sizes_for_mergesplit_step(int split_num, int num_regions){
    // make it 
    return;
}


AnyRegionMMDSplittingSchedule::AnyRegionMMDSplittingSchedule(
            const int num_splits, const int ndists,
            const int total_seats, std::vector<int> const &district_seat_sizes
):
    SplittingSchedule(SplittingSizeScheduleType::AnyValidSizeMMD, ndists, total_seats, district_seat_sizes),
    reachable_size(total_seats + 1, false){


    // the smallest size divided by total seats must be at least ndists or its not even
    // possible to split this map. 
    // For example: total_seats=6, ndists=4, district_seat_sizes = {2,3}
    if(total_seats / static_cast<double>(smallest_district_size) < ndists){
        throw Rcpp::exception("It is not possible to split this map with the given district sizes! There's not enough regions\n");
    }

    // Now we define the initial batch of reachable sizes as sizes that are a sum 
    for (auto const &district_size1: district_seat_sizes){
        for (auto const &district_size2: district_seat_sizes){
            if(district_size1 + district_size2 <= total_seats){
                reachable_size[district_size1 + district_size2] = true;
                // REprintf("%d is reachable!\n", district_size1 + district_size2);
            }
        }
    }


    REprintf("Smallest size is %d!\n", smallest_district_size);
    // Now for each possible region size we figure out 
    //  - the size it can be split into 
    //  - the min and max smallest possible sizes 

    for (int region_size = smallest_district_size+1; region_size <= total_seats; region_size++){
        // REprintf("Starting size %d!\n", region_size);
        // first check that this sum is reachable, if not then skip 
        if(!reachable_size[region_size]) continue; 
        // Now we know this size is reachable so first we update the reachable sizes 
        for (int potential_size = smallest_district_size; potential_size <= total_seats; potential_size++)
        {
            // skip if this size is not reachable or a district 
            if(!reachable_size[potential_size] && !is_district[potential_size]) continue; 
            // Else if its reachable and the sum is leq total seats then that size is reachable
            auto merged_size = region_size + potential_size;
            if(merged_size <= total_seats) reachable_size[merged_size] = true;
            // Now check if the equivalent 
            auto other_split_size = region_size - potential_size;
            // skip if the other size isn't possible 
        }
        
        // Now we want to figure out what reachable sizes can be combined to get this one
        // which is equivalent to what splits are ok. We'll only go up to the floor of 
        // the size over two which is the largest any smaller cut size can be 
        int smaller_cut_size_max = std::floor(region_size / 2);
        int smallest_possible_cut_size = total_seats;
        int largest_possible_cut_size = 0;

        // REprintf("Found all reachable ones! Now going from %d to %d\n", 
        // smallest_district_size, smaller_cut_size_max);


        for (int potential_size = smallest_district_size; potential_size <= smaller_cut_size_max; potential_size++){
            // REprintf("ALERT: Trying size %d: ", NEW_potential_size);
            // skip if this size is not reachable 
            // if(!reachable_size[NEW_potential_size] && !is_district[NEW_potential_size]) REprintf("Skipping %d!\n", NEW_potential_size);
            if(!reachable_size[potential_size] && !is_district[potential_size]) continue; 
            // skip if the other size is not reachable
            int other_split_size = region_size - potential_size;
            // if(!reachable_size[NEW_potential_size] && !is_district[NEW_potential_size]) REprintf("Skipping %d!\n", NEW_potential_size);
            if(!reachable_size[other_split_size] && !is_district[other_split_size]) continue; 
            // // If both sizes are reachable then this is an acceptable split 
            all_regions_smaller_cut_sizes_to_try[region_size].push_back(potential_size);
            // REprintf("Adding %d and %d!\n", NEW_potential_size, other_split_size);
            // // update the smallest and largest possible cut sizes 
            smallest_possible_cut_size = std::min(smallest_possible_cut_size, potential_size);
            largest_possible_cut_size = std::max(largest_possible_cut_size, other_split_size);
        }

        // Now we set these
        all_regions_min_and_max_possible_cut_sizes[region_size] = {
            smallest_possible_cut_size, largest_possible_cut_size
        };
        // REprintf("Done with size %d!\n", region_size);
        
    }
    // If the whole map isn't reachable then stop
    if(!reachable_size[total_seats]){
        throw Rcpp::exception("It is not possible to split this map with the given district sizes!\n");
    }
    REprintf("All done!\n");

}


void AnyRegionMMDSplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions
){
    // The biggest size a region can be is if all but one region is a district and
    // all districts are the smallest. Note this might not be possible but its just an UB

    // TODO might need to make this more precise!
    int presplit_biggest_possible_size = total_seats - (presplit_num_regions - 1) * smallest_district_size;
    REprintf("Its %d!\n", presplit_biggest_possible_size);

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    // MAKE ALL DISTRICT SIZES VALID SPLIT REGION SIZES 
    for (auto const &district_size: district_seat_sizes){
        valid_split_region_sizes[district_size] = true;
    }

    for (size_t region_size = smallest_district_size; region_size <=  presplit_biggest_possible_size; region_size++)
    {
        // We can split any reachable region 
        if(!reachable_size[region_size]) continue;
        valid_region_sizes_to_split[region_size] = true;
        // a valid split size is any reachable size or district where its partner is also district or
        // reachable and they sum to this region size so we can just iterate over the split sizes
        for(auto const &split_size1: all_regions_smaller_cut_sizes_to_try[region_size]){
            valid_split_region_sizes[split_size1] = true;
            valid_split_region_sizes[region_size - split_size1] = true;
        }
    }
    
    check_adj_to_regions = valid_split_region_sizes;

    
    // valid merge pairs are anything where both sizes are reachable 
    // TODO this
    for (size_t region1_size = smallest_district_size; region1_size < presplit_biggest_possible_size; region1_size++){
        if(!reachable_size[region1_size] && !is_district[region1_size]) continue;
        for (size_t region2_size = smallest_district_size; region2_size < presplit_biggest_possible_size; region2_size++){
            if(!reachable_size[region2_size] && !is_district[region2_size]){
                continue;
            }else{
                // since its a double loop we don't need to set twice since condition being checked is symmetric
                valid_merge_pair_sizes[region1_size][region2_size] = true;
            }
        }
    }
    
    return;
}

PureMSSplittingSchedule::PureMSSplittingSchedule(
    const int ndists, const int total_size,
    std::vector<int> const &district_seat_sizes
):
    SplittingSchedule(SplittingSizeScheduleType::PureMergeSplitSize, ndists, total_size, district_seat_sizes){
    

    for(int district_size1 = smallest_district_size;
            district_size1 <= largest_district_size;
            district_size1++){
        // mark each district size as a valid split region size
        valid_split_region_sizes[district_size1] = true;

        // don't support irregular district sizes like this right now
        if(district_size1*2 > ndists){
            throw Rcpp::exception("Pure MS not supported when one district bigger than half total seats\n");
        }

        for(int district_size2 = smallest_district_size; 
            district_size2 <= largest_district_size;
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
    for (size_t region_size = 1; region_size <= total_seats; region_size++)
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
): SplittingSchedule(SplittingSizeScheduleType::OneCustomSize, ndists, ndists, std::vector<int> {1}){

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
    const int total_seats,
    std::vector<int> const &district_seat_sizes,
    SplittingSizeScheduleType const schedule_type,
    Rcpp::List const &control
){
    if(schedule_type == SplittingSizeScheduleType::DistrictOnlySMD){
        return std::make_unique<DistrictOnlySplittingSchedule>(num_splits, ndists);
    }else if(schedule_type == SplittingSizeScheduleType::AnyValidSizeSMD){
        return std::make_unique<AnyRegionSMDSplittingSchedule>(num_splits, ndists);
    }else if(schedule_type == SplittingSizeScheduleType::OneCustomSize){
        return std::make_unique<OneCustomSplitSchedule>(num_splits, ndists, control);
    }else if(schedule_type == SplittingSizeScheduleType::DistrictOnlyMMD){
        return std::make_unique<DistrictOnlyMMDSplittingSchedule>(num_splits, ndists, total_seats, district_seat_sizes);
    }else if(schedule_type == SplittingSizeScheduleType::AnyValidSizeMMD){
        return std::make_unique<AnyRegionMMDSplittingSchedule>(num_splits, ndists, total_seats, district_seat_sizes);
    }else if(schedule_type == SplittingSizeScheduleType::CustomSizes){
        Rprintf("Not implemented!");
        throw Rcpp::exception("Schedule not impliemented yet!");
    }else{
        throw Rcpp::exception("Schedule not implemented yet!");
    }
}