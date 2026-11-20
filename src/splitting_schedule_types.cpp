/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/01
 * Purpose: Implements Various Splitting Schedule Types
 ********************************************************/

#include "splitting_schedule_types.h"
#include "base_plan_type.h" // Plan
#include "redist_constants.h" // perf_config

// for `Rprintf` / `REprintf`
#include <R_ext/Print.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>

constexpr bool DEBUG_SPLITTING_SCHEDULES_VERBOSE = false;

DistrictSizeSums::DistrictSizeSums(int const total_seats, int const ndists,
                                   std::vector<int> const &district_seat_sizes)
    : total_seats(total_seats), ndists(ndists),
      smallest_district_size(
          *std::min_element(district_seat_sizes.begin(), district_seat_sizes.end())),
      largest_district_size(
          *std::max_element(district_seat_sizes.begin(), district_seat_sizes.end())),
      seat_sizes_contiguous([&district_seat_sizes]() {
          // contiguous if the distinct sizes fill the range from smallest to largest
          std::vector<int> sizes(district_seat_sizes);
          std::sort(sizes.begin(), sizes.end());
          sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
          return sizes.back() - sizes.front() + 1 == static_cast<int>(sizes.size());
      }()) {}

bool DistrictSizeSums::table_is_sum_of(int const region_size, int const num_districts) const {
    throw std::runtime_error(
        "DistrictSizeSums::table_is_sum_of: district seat sizes that are not a contiguous "
        "range are not supported yet.\n");
}

void SplittingSchedule::check_context(Plan const &plan, ScheduleContext const &ctx) const {
    if (ctx.num_regions != plan.num_regions) {
        std::ostringstream oss;
        oss << "ScheduleContext was built for a plan with " << ctx.num_regions
            << " regions but was used with a plan with " << plan.num_regions
            << " regions. Rebuild the context with `build_context` after the plan changes.\n";
        throw std::runtime_error(oss.str());
    }
}

void SplittingSchedule::build_context(Plan const &plan, ScheduleContext &ctx) const {
    ctx.num_regions = plan.num_regions;
    // schedules that only use the size tables have nothing to compute
    if (!plan_context_needed) return;
    fill_context(plan, ctx);
}

/*
 * Whether some district size is itself a sum of two or more district sizes
 * (eg sizes 3 to 6, where 6 = 3 + 3). When none is, the district only MMD
 * schedule can be answered from the size tables alone; otherwise which
 * regions can be split depends on the rest of the plan.
 *
 * For a contiguous range this is just largest >= 2 * smallest. Non-contiguous
 * sizes conservatively return true, since building a context that isn't
 * needed is always safe.
 */
static bool some_district_size_is_sum_of_districts(std::vector<int> const &district_seat_sizes) {
    DistrictSizeSums const size_sums(0, 0, district_seat_sizes);
    if (!size_sums.seat_sizes_contiguous) return true;
    return size_sums.largest_district_size >= 2 * size_sums.smallest_district_size;
}

bool SplittingSchedule::region_is_splittable(Plan const &plan, ScheduleContext const &ctx,
                                             int const region_id) const {
    if constexpr (perf_config::object_integrity_checking) check_context(plan, ctx);
    return region_is_splittable_impl(plan, ctx, region_id);
}

CutSizeOptions SplittingSchedule::region_cut_sizes(Plan const &plan, ScheduleContext &ctx,
                                                   int const region_id) const {
    if constexpr (perf_config::object_integrity_checking) check_context(plan, ctx);
    return region_cut_sizes_impl(plan, ctx, region_id);
}

CutSizeOptions SplittingSchedule::merged_region_cut_sizes(Plan const &plan,
                                                          ScheduleContext &ctx,
                                                          int const region1_id,
                                                          int const region2_id) const {
    if constexpr (perf_config::object_integrity_checking) check_context(plan, ctx);
    return merged_region_cut_sizes_impl(plan, ctx, region1_id, region2_id);
}

bool SplittingSchedule::is_valid_merge_pair(Plan const &plan, ScheduleContext const &ctx,
                                            int const region1_id, int const region2_id) const {
    if constexpr (perf_config::object_integrity_checking) check_context(plan, ctx);
    return is_valid_merge_pair_impl(plan, ctx, region1_id, region2_id);
}

double SplittingSchedule::log_merged_region_selection_prob(Plan const &plan,
                                                           ScheduleContext const &ctx,
                                                           int const region1_id,
                                                           int const region2_id,
                                                           double const selection_alpha) const {
    if constexpr (perf_config::object_integrity_checking) check_context(plan, ctx);
    return log_merged_region_selection_prob_impl(plan, ctx, region1_id, region2_id,
                                                 selection_alpha);
}

bool SplittingSchedule::region_is_splittable_impl(Plan const &plan, ScheduleContext const &ctx,
                                                  int const region_id) const {
    return valid_region_sizes_to_split[plan.region_sizes[region_id]];
}

CutSizeOptions SplittingSchedule::region_cut_sizes_impl(Plan const &plan, ScheduleContext &ctx,
                                                        int const region_id) const {
    return size_cut_options(plan.region_sizes[region_id]);
}

CutSizeOptions SplittingSchedule::merged_region_cut_sizes_impl(Plan const &plan,
                                                               ScheduleContext &ctx,
                                                               int const region1_id,
                                                               int const region2_id) const {
    return size_cut_options(plan.region_sizes[region1_id] + plan.region_sizes[region2_id]);
}

bool SplittingSchedule::is_valid_merge_pair_impl(Plan const &plan, ScheduleContext const &ctx,
                                                 int const region1_id,
                                                 int const region2_id) const {
    return valid_merge_pair_sizes[plan.region_sizes[region1_id]]
                                 [plan.region_sizes[region2_id]];
}

/*
 * Regions are picked to split with probability proportional to size^alpha
 * among the splittable regions. This is the probability the region made by
 * merging `region1_id` and `region2_id` would have been picked in the plan
 * where they are merged.
 */
double SplittingSchedule::log_merged_region_selection_prob_impl(
    Plan const &plan, ScheduleContext const &ctx, int const region1_id, int const region2_id,
    double const selection_alpha) const {
    // weight of the merged region
    auto merged_region_size = plan.region_sizes[region1_id] + plan.region_sizes[region2_id];
    double merged_region_wgt = std::pow(merged_region_size, selection_alpha);
    double wgt_sum = merged_region_wgt;

    for (int region_id = 0; region_id < plan.num_regions; region_id++) {
        auto region_size = plan.region_sizes[region_id];
        // add if splittable and not one of the two being merged
        if (valid_region_sizes_to_split[region_size] && region_id != region1_id &&
            region_id != region2_id) {
            wgt_sum += std::pow(region_size, selection_alpha);
        }
    }

    return std::log(merged_region_wgt) - std::log(wgt_sum);
}

void SplittingSchedule::reset_splitting_and_merge_booleans() {
    // reset all regions that can be split to false
    std::fill(valid_region_sizes_to_split.begin(), valid_region_sizes_to_split.end(), false);
    std::fill(valid_split_region_sizes.begin(), valid_split_region_sizes.end(), false);
    // reset merge pairs and split sizes to try
    for (size_t region_size = 1; region_size <= total_seats; region_size++) {
        std::fill(valid_merge_pair_sizes[region_size].begin(),
                  valid_merge_pair_sizes[region_size].end(), false);
    }
}

void SplittingSchedule::print_current_step_splitting_info() {
    Rprintf("District Seat sizes are: (");
    for (auto const &v : district_seat_sizes) {
        Rprintf("%d, ", v);
    }
    Rprintf(")\n");

    Rprintf("The following region sizes can split:\n");
    for (size_t i = 1; i <= total_seats; i++) {
        if (valid_region_sizes_to_split[i]) {
            Rprintf("\tRegion Size %d | ", (int)i);
            Rprintf("min/max (%d, %d)", all_regions_min_and_max_possible_cut_sizes[i].first,
                    all_regions_min_and_max_possible_cut_sizes[i].second);
            Rprintf(" | possible split sizes: ");
            for (auto smaller_size : all_regions_smaller_cut_sizes_to_try[i]) {
                Rprintf("(%d, %d) ", (int)smaller_size, static_cast<int>(i - smaller_size));
            }
            Rprintf("\n");
        }
    }
    Rprintf("All Possible Split Sizes are:(");
    for (int i = 1; i <= total_seats; i++) {
        if (valid_split_region_sizes[i]) {
            Rprintf("%d, ", i);
        }
    }
    Rprintf(")\n");

    Rprintf("The following valid merge pairs are:\n");
    Rprintf("_ |");
    for (int i = 1; i <= total_seats; i++) {
        if (!valid_region_sizes_to_split[i] && !valid_split_region_sizes[i])
            continue;
        // \033[4m .. \033[0m underlines the column header
        Rprintf(" \033[4m%d\033[0m", i);
    }
    Rprintf("\n");
    for (int i = 1; i <= total_seats; i++) {
        if (!valid_region_sizes_to_split[i] && !valid_split_region_sizes[i])
            continue;
        Rprintf("%d | ", i);
        for (int j = 1; j <= total_seats; j++) {
            if (!valid_region_sizes_to_split[j] && !valid_split_region_sizes[j])
                continue;
            Rprintf("%s", valid_merge_pair_sizes[i][j] ? "1 " : "0 ");
        }
        Rprintf("\n");
    }
}

/*
 * Constructor for district splits only splitting schedule
 */
DistrictOnlySplittingSchedule::DistrictOnlySplittingSchedule(const int num_splits,
                                                             const int ndists)
    : SplittingSchedule(SplittingSizeScheduleType::DistrictOnlySMD, ndists, ndists,
                        std::vector<int>{1}, false) {
    // For each size above 2 make the min=1 and the max=size-1
    for (int region_size = 2; region_size <= ndists; region_size++) {
        all_regions_smaller_cut_sizes_to_try[region_size] = std::vector<int>{1};
        // smallest possible cut size is 1 and biggest is region_size-1
        all_regions_min_and_max_possible_cut_sizes[region_size] = {1, region_size - 1};
    }
}

/*
 * Constructor for district splits only splitting schedule
 */
void DistrictOnlySplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions) {
    // The biggest size a region can be is if all but one region is a district
    int presplit_remainder_size = ndists - presplit_num_regions + 1;
    int after_split_remainder_size = presplit_remainder_size - 1;

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

    return;
}

// makes it so districts can be merged during merge split
void DistrictOnlySplittingSchedule::update_cut_sizes_for_mergesplit_step(int split_num,
                                                                         int num_regions) {
    // if only one region then do nothing because only 1 district (if ndists=2 then we don't
    // need to worry)
    if (num_regions == 2)
        return;
    // make it so we can merge districts
    valid_merge_pair_sizes[1][1] = true;
    // make it so we can split a region of size 2
    valid_region_sizes_to_split[2] = true;
    return;
}

AnyRegionSMDSplittingSchedule::AnyRegionSMDSplittingSchedule(const int num_splits,
                                                             const int ndists)
    : SplittingSchedule(SplittingSizeScheduleType::AnyValidSizeSMD, ndists, ndists,
                        std::vector<int>{1}, false) {
    // If doing any sizes then for `region_size` the possible split sizes are always
    // if any sizes allowed then for 2, ..., ndists make it 1,...,floor(size/2)
    for (int region_size = 2; region_size <= ndists; region_size++) {
        // Get the lrgest smaller cut size we'll try which is the
        // floor of the size of region over 2
        int smaller_cut_size_max = std::floor(region_size / 2);
        all_regions_smaller_cut_sizes_to_try[region_size].resize(smaller_cut_size_max);
        // Now makes the entry for region size go 1,...smaller_cut_size_max
        std::iota(all_regions_smaller_cut_sizes_to_try[region_size].begin(),
                  all_regions_smaller_cut_sizes_to_try[region_size].end(), 1);
        // smallest possible cut size is 1 and biggest is region_size-1
        all_regions_min_and_max_possible_cut_sizes[region_size] = {1, region_size - 1};
    }
}

void AnyRegionSMDSplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions) {
    // The biggest size a region can be is if all but one region is a district
    int presplit_biggest_possible_size = ndists - presplit_num_regions + 1;

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    for (size_t region_size = 2; region_size <= presplit_biggest_possible_size; region_size++) {
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
    // valid merge pairs are anything that sums to less than or equal to the max valid size
    for (size_t region1_size = 1; region1_size < presplit_biggest_possible_size;
         region1_size++) {
        for (size_t region2_size = 1; region2_size < presplit_biggest_possible_size;
             region2_size++) {
            if (region1_size + region2_size <= presplit_biggest_possible_size) {
                valid_merge_pair_sizes[region1_size][region2_size] = true;
            }
        }
    }

    return;
}

DistrictOnlyMMDSplittingSchedule::DistrictOnlyMMDSplittingSchedule(
    const int num_splits, const int ndists, const int total_seats,
    std::vector<int> const &district_seat_sizes)
    : SplittingSchedule(SplittingSizeScheduleType::DistrictOnlyMMD, ndists, total_seats,
                        district_seat_sizes,
                        some_district_size_is_sum_of_districts(district_seat_sizes)) {

    // right now only ranges are supported
    if (!size_sums.seat_sizes_contiguous) {
        std::ostringstream oss;
        oss << "For multi-member districting with district only splits the district seat "
            << "sizes must be a consecutive range of sizes with no gaps, but got sizes: ";
        for (int const size : district_seat_sizes) oss << size << " ";
        oss << "\n";
        throw std::runtime_error(oss.str());
    }

    // For a range it must be smallest_district_size <= total_seats/ndists <=
    // largest_district_size
    if (smallest_district_size > total_seats / static_cast<double>(ndists)) {
        std::ostringstream oss;
        oss << "It is not possible to split " << total_seats << " seats into " << ndists
            << " districts: the smallest district size (" << smallest_district_size
            << " seats) is larger than the average of " << total_seats / static_cast<double>(ndists)
            << " seats per district.\n";
        throw std::runtime_error(oss.str());
    }
    if (largest_district_size < total_seats / static_cast<double>(ndists)) {
        std::ostringstream oss;
        oss << "It is not possible to split " << total_seats << " seats into " << ndists
            << " districts: the largest district size (" << largest_district_size
            << " seats) is smaller than the average of " << total_seats / static_cast<double>(ndists)
            << " seats per district.\n";
        throw std::runtime_error(oss.str());
    }

    // For merge split: the ways to split two merged districts back into two
    // districts. These only depend on the merged size, not the step.
    district_pair_smaller_cut_sizes.resize(total_seats + 1);
    district_pair_min_and_max_cut_sizes.resize(total_seats + 1, {0, 0});
    for (int const district1_size : district_seat_sizes) {
        for (int const district2_size : district_seat_sizes) {
            int const merged_size = district1_size + district2_size;
            if (merged_size > total_seats)
                continue;
            district_pair_smaller_cut_sizes[merged_size].push_back(
                std::min(district1_size, district2_size));
        }
    }
    for (int merged_size = 1; merged_size <= total_seats; merged_size++) {
        auto &cut_sizes = district_pair_smaller_cut_sizes[merged_size];
        if (cut_sizes.empty())
            continue;
        std::sort(cut_sizes.begin(), cut_sizes.end());
        cut_sizes.erase(std::unique(cut_sizes.begin(), cut_sizes.end()), cut_sizes.end());
        district_pair_min_and_max_cut_sizes[merged_size] = {cut_sizes[0],
                                                            merged_size - cut_sizes[0]};
    }
}

void DistrictOnlyMMDSplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions) {
    // get current number of districts before splitting
    int presplit_ndists = presplit_num_regions - 1;
    // get the number of districts the remainder will need to be split into
    // which is the number of districts minus the number of districts currently
    int remainder_ndists = ndists - presplit_ndists;
    // Define largest and smallest possible remainder sizes

    // the context shared by the errors below
    auto schedule_error = [&](std::string const &problem) {
        std::ostringstream oss;
        oss << "DistrictOnlyMMDSplittingSchedule::set_potential_cut_sizes_for_each_valid_size: "
            << problem << "\n"
            << "split_num=" << split_num << ", presplit_num_regions=" << presplit_num_regions
            << ", total_seats=" << total_seats << ", ndists=" << ndists
            << ", district sizes " << smallest_district_size << " to " << largest_district_size
            << ", remainder region must hold " << remainder_ndists << " districts.\n";
        return std::runtime_error(oss.str());
    };

    if (presplit_num_regions < 1 || presplit_num_regions >= ndists) {
        throw schedule_error("the number of presplit regions must be between 1 and ndists - 1.");
    }

    // district pair merges are only for merge split rounds
    mergesplit_district_pairs = false;

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    // A remainder of size s is possible if s can be split into `remainder_ndists`
    // districts and the other `presplit_ndists` districts can hold the rest of the seats.
    bool any_remainder_size = false;
    for (int remainder_size = 1; remainder_size <= total_seats; remainder_size++) {
        if (!size_sums.is_sum_of(remainder_size, remainder_ndists) ||
            !size_sums.is_sum_of(total_seats - remainder_size, presplit_ndists)) {
            continue;
        }
        any_remainder_size = true;
        // clear split sizes to try
        all_regions_smaller_cut_sizes_to_try[remainder_size].clear();
        // We can split off any district that leaves a remainder that can still be
        // split into the remaining districts
        for (int const district_size : district_seat_sizes) {
            if (district_size >= remainder_size)
                continue;
            int new_remainder_size = remainder_size - district_size;
            if (!size_sums.is_sum_of(new_remainder_size, remainder_ndists - 1))
                continue;

            // Add the smaller of the two cut pieces
            all_regions_smaller_cut_sizes_to_try[remainder_size].push_back(
                std::min(new_remainder_size, district_size));
            // Now also add this district as a possible split size
            valid_split_region_sizes[district_size] = true;
            valid_split_region_sizes[new_remainder_size] = true;
            // The district and the new remainder are also valid merge pairs
            valid_merge_pair_sizes[new_remainder_size][district_size] = true;
            valid_merge_pair_sizes[district_size][new_remainder_size] = true;
        }

        if (all_regions_smaller_cut_sizes_to_try[remainder_size].size() == 0)
            continue;
        valid_region_sizes_to_split[remainder_size] = true;
        // Now sort the vector and make it unique
        std::sort(all_regions_smaller_cut_sizes_to_try[remainder_size].begin(),
                  all_regions_smaller_cut_sizes_to_try[remainder_size].end());
        auto last = std::unique(all_regions_smaller_cut_sizes_to_try[remainder_size].begin(),
                                all_regions_smaller_cut_sizes_to_try[remainder_size].end());
        all_regions_smaller_cut_sizes_to_try[remainder_size].erase(
            last, all_regions_smaller_cut_sizes_to_try[remainder_size].end());
        // now find the min and max values
        all_regions_min_and_max_possible_cut_sizes[remainder_size] = {
            all_regions_smaller_cut_sizes_to_try[remainder_size][0],
            remainder_size - all_regions_smaller_cut_sizes_to_try[remainder_size][0]};
    }

    if (!any_remainder_size) {
        throw schedule_error("no remainder region size can be split into the remaining "
                             "districts.");
    }
}

void DistrictOnlyMMDSplittingSchedule::update_cut_sizes_for_mergesplit_step(int split_num,
                                                                            int num_regions) {
    // let adjacent districts be merged and resplit during merge split rounds
    mergesplit_district_pairs = true;
}

bool DistrictOnlyMMDSplittingSchedule::is_district_pair_merge(Plan const &plan,
                                                              ScheduleContext const &ctx,
                                                              int const region1_id,
                                                              int const region2_id) const {
    if (!mergesplit_district_pairs)
        return false;
    if (!is_district[plan.region_sizes[region1_id]] ||
        !is_district[plan.region_sizes[region2_id]])
        return false;
    // The plan needs a non district sized remainder, which then isn't one of the
    // two districts. If no district size is a sum of other district sizes
    // every plan short of a full plan has one.
    if (plan_context_needed)
        return ctx.non_district_region >= 0;
    return plan.num_regions < ndists;
}

CutSizeOptions DistrictOnlyMMDSplittingSchedule::merged_region_cut_sizes_impl(
    Plan const &plan, ScheduleContext &ctx, int const region1_id, int const region2_id) const {
    if (is_district_pair_merge(plan, ctx, region1_id, region2_id)) {
        // two districts can only be resplit into two districts
        int const merged_size = plan.region_sizes[region1_id] + plan.region_sizes[region2_id];
        return {district_pair_smaller_cut_sizes[merged_size],
                district_pair_min_and_max_cut_sizes[merged_size].first,
                district_pair_min_and_max_cut_sizes[merged_size].second};
    }
    return SplittingSchedule::merged_region_cut_sizes_impl(plan, ctx, region1_id, region2_id);
}

void DistrictOnlyMMDSplittingSchedule::fill_context(Plan const &plan,
                                                    ScheduleContext &ctx) const {
    // find the one region that isn't district sized, if there is one
    ctx.non_district_region = -1;
    for (int region_id = 0; region_id < plan.num_regions; region_id++) {
        if (is_district[plan.region_sizes[region_id]])
            continue;
        if (ctx.non_district_region >= 0) {
            std::ostringstream oss;
            oss << "DistrictOnlyMMDSplittingSchedule: plan has more than one region that is not "
                << "district sized (regions " << ctx.non_district_region << " and " << region_id
                << ") so it can't be completed by splitting off one district at a time.\n"
                << plan.debug_string(true);
            throw std::runtime_error(oss.str());
        }
        ctx.non_district_region = region_id;
    }
}

bool DistrictOnlyMMDSplittingSchedule::region_is_splittable_impl(Plan const &plan,
                                                                 ScheduleContext const &ctx,
                                                                 int const region_id) const {
    bool const size_splittable =
        SplittingSchedule::region_is_splittable_impl(plan, ctx, region_id);
    // if every region is district sized any splittable size can be the remainder
    if (!plan_context_needed || ctx.non_district_region < 0)
        return size_splittable;
    // otherwise the non district sized region is the remainder and the only
    // region that can be split
    return region_id == ctx.non_district_region && size_splittable;
}

bool DistrictOnlyMMDSplittingSchedule::is_valid_merge_pair_impl(Plan const &plan,
                                                                ScheduleContext const &ctx,
                                                                int const region1_id,
                                                                int const region2_id) const {
    // during merge split two districts can be merged if a remainder exists elsewhere
    if (is_district_pair_merge(plan, ctx, region1_id, region2_id))
        return true;

    bool const sizes_ok =
        SplittingSchedule::is_valid_merge_pair_impl(plan, ctx, region1_id, region2_id);
    if (!plan_context_needed || !sizes_ok)
        return sizes_ok;
    // If one region isn't district sized it is the remainder, and a split can
    // only have made it from the previous remainder, so it must be in the pair
    int const non_district_region = ctx.non_district_region;
    return non_district_region < 0 || region1_id == non_district_region ||
           region2_id == non_district_region;
}

double DistrictOnlyMMDSplittingSchedule::log_merged_region_selection_prob_impl(
    Plan const &plan, ScheduleContext const &ctx, int const region1_id, int const region2_id,
    double const selection_alpha) const {
    if (plan_context_needed) {
        // If the merged region isn't district sized it was the only region
        // that could be split, so it was picked with probability 1
        int const merged_region_size =
            plan.region_sizes[region1_id] + plan.region_sizes[region2_id];
        if (!is_district[merged_region_size])
            return 0.0;
    }
    // otherwise every region was district sized and it was picked from all
    // the regions whose size is splittable
    return SplittingSchedule::log_merged_region_selection_prob_impl(
        plan, ctx, region1_id, region2_id, selection_alpha);
}

void DistrictOnlyMMDSplittingSchedule::check_initial_plan(Plan const &plan) const {
    int const num_regions = plan.num_regions;
    // the number of districts the remainder still has to be split into
    int const remainder_ndists = ndists - num_regions + 1;

    auto plan_error = [&](std::string const &problem) {
        std::ostringstream oss;
        oss << "The initial partial plan can't be completed by splitting off one district at a "
            << "time: " << problem << "\n"
            << "num_regions=" << num_regions << ", ndists=" << ndists
            << ", total_seats=" << total_seats << ", district sizes " << smallest_district_size
            << " to " << largest_district_size << ".\n"
            << plan.debug_string(true);
        return std::runtime_error(oss.str());
    };

    if (num_regions < 1 || num_regions > ndists) {
        throw plan_error("the number of regions must be between 1 and ndists.");
    }

    int seat_sum = 0;
    int non_district_region = -1;
    bool any_possible_remainder = false;
    for (int region_id = 0; region_id < num_regions; region_id++) {
        int const region_size = plan.region_sizes[region_id];
        if (region_size < 1 || region_size > total_seats) {
            throw plan_error("region " + std::to_string(region_id) + " has " +
                             std::to_string(region_size) + " seats.");
        }
        seat_sum += region_size;
        if (size_sums.is_sum_of(region_size, remainder_ndists))
            any_possible_remainder = true;
        if (is_district[region_size])
            continue;
        if (non_district_region >= 0) {
            throw plan_error("more than one region is not district sized.");
        }
        non_district_region = region_id;
    }

    if (seat_sum != total_seats) {
        throw plan_error("the region seats sum to " + std::to_string(seat_sum) + ".");
    }
    if (non_district_region >= 0) {
        // the non district sized region has to be the remainder
        int const remainder_size = plan.region_sizes[non_district_region];
        if (!size_sums.is_sum_of(remainder_size, remainder_ndists)) {
            throw plan_error("region " + std::to_string(non_district_region) + " with " +
                             std::to_string(remainder_size) + " seats can't be split into " +
                             std::to_string(remainder_ndists) + " districts.");
        }
    } else if (!any_possible_remainder) {
        throw plan_error("no region can be split into the remaining " +
                         std::to_string(remainder_ndists) + " districts.");
    }
}

AnyRegionMMDSplittingSchedule::AnyRegionMMDSplittingSchedule(
    const int num_splits, const int ndists, const int total_seats,
    std::vector<int> const &district_seat_sizes)
    : SplittingSchedule(SplittingSizeScheduleType::AnyValidSizeMMD, ndists, total_seats,
                        district_seat_sizes, true),
      reachable_size(total_seats + 1, false) {

    // the smallest size divided by total seats must be at least ndists or its not even
    // possible to split this map.
    // For example: total_seats=6, ndists=4, district_seat_sizes = {2,3}
    if (total_seats / static_cast<double>(smallest_district_size) < ndists) {
        throw std::runtime_error("It is not possible to split this map with the given district "
                              "sizes! There's not enough regions\n");
    }

    // Now we define the initial batch of reachable sizes as sizes that are a sum
    for (auto const &district_size1 : district_seat_sizes) {
        for (auto const &district_size2 : district_seat_sizes) {
            if (district_size1 + district_size2 <= total_seats) {
                reachable_size[district_size1 + district_size2] = true;
                // REprintf("%d is reachable!\n", district_size1 + district_size2);
            }
        }
    }

    REprintf("Smallest size is %d!\n", smallest_district_size);
    // Now for each possible region size we figure out
    //  - the size it can be split into
    //  - the min and max smallest possible sizes

    for (int region_size = smallest_district_size + 1; region_size <= total_seats;
         region_size++) {
        // REprintf("Starting size %d!\n", region_size);
        // first check that this sum is reachable, if not then skip
        if (!reachable_size[region_size])
            continue;
        // Now we know this size is reachable so first we update the reachable sizes
        for (int potential_size = smallest_district_size; potential_size <= total_seats;
             potential_size++) {
            // skip if this size is not reachable or a district
            if (!reachable_size[potential_size] && !is_district[potential_size])
                continue;
            // Else if its reachable and the sum is leq total seats then that size is reachable
            auto merged_size = region_size + potential_size;
            if (merged_size <= total_seats)
                reachable_size[merged_size] = true;
        }

        // Now we want to figure out what reachable sizes can be combined to get this one
        // which is equivalent to what splits are ok. We'll only go up to the floor of
        // the size over two which is the largest any smaller cut size can be
        int smaller_cut_size_max = std::floor(region_size / 2);
        int smallest_possible_cut_size = total_seats;
        int largest_possible_cut_size = 0;

        // REprintf("Found all reachable ones! Now going from %d to %d\n",
        // smallest_district_size, smaller_cut_size_max);

        for (int potential_size = smallest_district_size;
             potential_size <= smaller_cut_size_max; potential_size++) {
            // REprintf("ALERT: Trying size %d: ", NEW_potential_size);
            // skip if this size is not reachable
            // if(!reachable_size[NEW_potential_size] && !is_district[NEW_potential_size])
            // REprintf("Skipping %d!\n", NEW_potential_size);
            if (!reachable_size[potential_size] && !is_district[potential_size])
                continue;
            // skip if the other size is not reachable
            int other_split_size = region_size - potential_size;
            // if(!reachable_size[NEW_potential_size] && !is_district[NEW_potential_size])
            // REprintf("Skipping %d!\n", NEW_potential_size);
            if (!reachable_size[other_split_size] && !is_district[other_split_size])
                continue;
            // // If both sizes are reachable then this is an acceptable split
            all_regions_smaller_cut_sizes_to_try[region_size].push_back(potential_size);
            // REprintf("Adding %d and %d!\n", NEW_potential_size, other_split_size);
            // // update the smallest and largest possible cut sizes
            smallest_possible_cut_size = std::min(smallest_possible_cut_size, potential_size);
            largest_possible_cut_size = std::max(largest_possible_cut_size, other_split_size);
        }

        // Now we set these
        all_regions_min_and_max_possible_cut_sizes[region_size] = {smallest_possible_cut_size,
                                                                   largest_possible_cut_size};
        // REprintf("Done with size %d!\n", region_size);
    }
    // If the whole map isn't reachable then stop
    if (!reachable_size[total_seats]) {
        throw std::runtime_error(
            "It is not possible to split this map with the given district sizes!\n");
    }
    REprintf("All done!\n");
}

void AnyRegionMMDSplittingSchedule::set_potential_cut_sizes_for_each_valid_size(
    int split_num, int presplit_num_regions) {
    // The biggest size a region can be is if all but one region is a district and
    // all districts are the smallest. Note this might not be possible but its just an UB

    // TODO might need to make this more precise!
    int presplit_biggest_possible_size =
        total_seats - (presplit_num_regions - 1) * smallest_district_size;
    REprintf("Its %d!\n", presplit_biggest_possible_size);

    // reset all regions split and merge booleans
    reset_splitting_and_merge_booleans();

    // MAKE ALL DISTRICT SIZES VALID SPLIT REGION SIZES
    for (auto const &district_size : district_seat_sizes) {
        valid_split_region_sizes[district_size] = true;
    }

    for (size_t region_size = smallest_district_size;
         region_size <= presplit_biggest_possible_size; region_size++) {
        // We can split any reachable region
        if (!reachable_size[region_size])
            continue;
        valid_region_sizes_to_split[region_size] = true;
        // a valid split size is any reachable size or district where its partner is also
        // district or reachable and they sum to this region size so we can just iterate over
        // the split sizes
        for (auto const &split_size1 : all_regions_smaller_cut_sizes_to_try[region_size]) {
            valid_split_region_sizes[split_size1] = true;
            valid_split_region_sizes[region_size - split_size1] = true;
        }
    }

    // valid merge pairs are anything where both sizes are reachable
    // TODO this
    for (size_t region1_size = smallest_district_size;
         region1_size < presplit_biggest_possible_size; region1_size++) {
        if (!reachable_size[region1_size] && !is_district[region1_size])
            continue;
        for (size_t region2_size = smallest_district_size;
             region2_size < presplit_biggest_possible_size; region2_size++) {
            if (!reachable_size[region2_size] && !is_district[region2_size]) {
                continue;
            } else {
                // since its a double loop we don't need to set twice since condition being
                // checked is symmetric
                valid_merge_pair_sizes[region1_size][region2_size] = true;
            }
        }
    }

    return;
}

PureMSSplittingSchedule::PureMSSplittingSchedule(const int ndists, const int total_size,
                                                 std::vector<int> const &district_seat_sizes)
    : SplittingSchedule(SplittingSizeScheduleType::PureMergeSplitSize, ndists, total_size,
                        district_seat_sizes, false) {

    for (int const district_size1: district_seat_sizes) {
        for (int const district_size2: district_seat_sizes) {

            if (district_size1 + district_size2 <= total_seats) {
                valid_split_region_sizes[district_size1] = true;
                valid_split_region_sizes[district_size2] = true;
                valid_merge_pair_sizes[district_size1][district_size2] = true;
                valid_merge_pair_sizes[district_size2][district_size1] = true;
                valid_region_sizes_to_split[district_size1 + district_size2] = true;
                // just add we'll sort and make unique later
                all_regions_smaller_cut_sizes_to_try[district_size1 + district_size2]
                    .emplace_back(std::min(district_size1, district_size2));
            }
        }
    }

    // make smaller cut sizes unique and sort
    // get min and max sizes


    for (size_t region_size = 1; region_size <= total_seats; region_size++) {
        if (!valid_region_sizes_to_split[region_size])
            continue;
        // sort and remove duplicates
        // 1) Sort the vector
        std::sort(all_regions_smaller_cut_sizes_to_try[region_size].begin(),
                  all_regions_smaller_cut_sizes_to_try[region_size].end());
        // remove duplicates
        all_regions_smaller_cut_sizes_to_try[region_size].erase(
            std::unique(all_regions_smaller_cut_sizes_to_try[region_size].begin(),
                        all_regions_smaller_cut_sizes_to_try[region_size].end()),
            all_regions_smaller_cut_sizes_to_try[region_size].end());
        // now the smallest possible size is the first one
        int smallest_possible_size = all_regions_smaller_cut_sizes_to_try[region_size][0];

        // set biggest and smallest possible size
        all_regions_min_and_max_possible_cut_sizes[region_size] = {
            smallest_possible_size, region_size - smallest_possible_size};
    }
}




// Returns a splitting schedule depending on the schedule type
std::unique_ptr<SplittingSchedule>
get_splitting_schedule(const int num_splits, const int ndists, const int total_seats,
                       std::vector<int> const &district_seat_sizes,
                       SplittingSizeScheduleType const schedule_type
                       ) {
    if (schedule_type == SplittingSizeScheduleType::DistrictOnlySMD) {
        return std::make_unique<DistrictOnlySplittingSchedule>(num_splits, ndists);
    } else if (schedule_type == SplittingSizeScheduleType::AnyValidSizeSMD) {
        return std::make_unique<AnyRegionSMDSplittingSchedule>(num_splits, ndists);
    } else if (schedule_type == SplittingSizeScheduleType::DistrictOnlyMMD) {
        return std::make_unique<DistrictOnlyMMDSplittingSchedule>(
            num_splits, ndists, total_seats, district_seat_sizes);
    } else if (schedule_type == SplittingSizeScheduleType::AnyValidSizeMMD) {
        return std::make_unique<AnyRegionMMDSplittingSchedule>(num_splits, ndists, total_seats,
                                                               district_seat_sizes);
    } else if (schedule_type == SplittingSizeScheduleType::CustomSizes) {
        throw std::runtime_error("The custom sizes splitting schedule is not implemented yet.");
    } else {
        throw std::runtime_error("No splitting schedule is implemented for splitting schedule "
                                 "type " +
                                 std::to_string(static_cast<int>(schedule_type)) + ".");
    }
}
