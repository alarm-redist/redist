#pragma once
#ifndef SPLITTING_SCHEDULE_TYPES_H
#define SPLITTING_SCHEDULE_TYPES_H

#include <memory>
#include "advanced_types.h"

class Plan;

/*
 * Answers whether a region with `region_size` seats can be divided into
 * exactly `num_districts` districts whose sizes are all allowed district
 * sizes, ie whether `region_size` is a sum of exactly `num_districts`
 * district sizes (with repetition).
 *
 * When the district sizes are a contiguous range [smallest, largest] this
 * is just `num_districts * smallest <= region_size <= num_districts * largest`.
 * Non-contiguous sizes can have gaps (eg sizes {3, 7} can make 21 from 3 or
 * 7 districts but not 4, 5, or 6) so they need a subset sum table instead.
 */
class DistrictSizeSums {
  public:
    DistrictSizeSums(int const total_seats, int const ndists,
                     std::vector<int> const &district_seat_sizes);

    int const total_seats;
    int const ndists;
    int const smallest_district_size;
    int const largest_district_size;
    // whether the district sizes are a contiguous range with no gaps
    bool const seat_sizes_contiguous;

    // whether `region_size` is a sum of exactly `num_districts` district sizes
    bool is_sum_of(int const region_size, int const num_districts) const {
        if (seat_sizes_contiguous) {
            return contiguous_is_sum_of(region_size, num_districts);
        } else {
            return table_is_sum_of(region_size, num_districts);
        }
    }

  private:
    bool contiguous_is_sum_of(int const region_size, int const num_districts) const {
        return num_districts * smallest_district_size <= region_size &&
               region_size <= num_districts * largest_district_size;
    }
    bool table_is_sum_of(int const region_size, int const num_districts) const;
};

/*
 * The allowed ways to split one region, which is either an existing region
 * or the region made by merging two regions.
 *
 *  - `smaller_cut_sizes` is the sorted list of sizes the smaller of the two
 *    new regions may have.
 *  - `min_cut_size` and `max_cut_size` are the smallest and largest sizes
 *    either new region may have, used for population bounds when drawing
 *    trees and searching for cuts.
 *
 * `smaller_cut_sizes` refers to storage owned by the schedule's size tables
 * or by the ScheduleContext passed to the query, so it is only valid until
 * that context is next queried or rebuilt.
 */
struct CutSizeOptions {
    std::vector<int> const &smaller_cut_sizes;
    int min_cut_size;
    int max_cut_size;
};

/*
 * Per plan information a splitting schedule caches while answering
 * questions about one plan.
 *
 * The schedule itself is shared read only by every thread, so anything that
 * depends on a particular plan (rather than just the current step) lives
 * here instead. Use one context per plan: build it with
 * `SplittingSchedule::build_context` and pass it with the plan to every
 * query. Rebuild it whenever the plan's regions change.
 */
class ScheduleContext {
  public:
    // number of regions in the plan the context was built for, used to catch
    // a context that wasn't rebuilt after its plan changed
    int num_regions = -1;
    // storage for plan dependent cut size lists returned in CutSizeOptions
    std::vector<int> cut_sizes_scratch;
    // District only MMD: the id of the one region whose size is not a district
    // size (the remainder), or -1 if every region is district sized
    int non_district_region = -1;
};

/*
 * Abstract class for object used to manage what sizes are allowed to be split
 * at each step.
 *
 * It is used in all three splitting, merge-split, and weight computation steps.
 *
 * The size tables below describe what is allowed for a region of a given
 * size at the current step. Questions about a particular plan should go
 * through the plan level queries instead, which by default answer from the
 * size tables but are overridden by schedules where the answer depends on
 * the rest of the plan.
 */
class SplittingSchedule {

  private:
    // this are both only actually tracked for custom splitting schedules
    std::vector<std::vector<bool>>
        valid_split_region_sizes_list; // the allowable splitting sizes for each split
    std::vector<std::vector<bool>>
        valid_presplit_region_sizes_list; // the allowable presplit region sizes for each split

  public:
    virtual ~SplittingSchedule() = default;

    // Base constructor just creates vectors, does not set them
    SplittingSchedule(const SplittingSizeScheduleType schedule_type, const int ndists,
                      const int total_seats, std::vector<int> const &district_seat_sizes,
                      bool const plan_context_needed)
        : schedule_type(schedule_type), ndists(ndists), total_seats(total_seats),
          plan_context_needed(plan_context_needed),
          smallest_district_size(
              *min_element(district_seat_sizes.begin(), district_seat_sizes.end())),
          largest_district_size(
              *max_element(district_seat_sizes.begin(), district_seat_sizes.end())),
          district_seat_sizes(district_seat_sizes),
          is_district([ndists, total_seats, &district_seat_sizes]() {
              // vector where index i is true iff i seats is a district
              std::vector<bool> is_district_vec(total_seats + 1, false);
              if (ndists > total_seats) {
                  throw std::runtime_error("The number of distrcts must be less than or equal to "
                                        "the total number of seats!\n");
              } else if (ndists != total_seats) {
                  for (auto const &a_size : district_seat_sizes) {
                      if (a_size <= 0)
                          throw std::runtime_error(
                              "District Seat Sizes must be strictly positive!\n");
                      if (a_size >= total_seats)
                          throw std::runtime_error(
                              "District Seat Sizes must be less than total seats!\n");
                      // mark this as a district size
                      is_district_vec[a_size] = true;
                  }

              } else {
                  is_district_vec[1] = true;
              }
              return is_district_vec;
          }()),
          is_mmd(ndists != total_seats),
          size_sums(total_seats, ndists, district_seat_sizes),
          all_regions_smaller_cut_sizes_to_try(total_seats + 1),
          all_regions_min_and_max_possible_cut_sizes(total_seats + 1),
          valid_region_sizes_to_split(total_seats + 1),
          valid_split_region_sizes(total_seats + 1),
          valid_merge_pair_sizes(total_seats + 1, std::vector<bool>(total_seats + 1, false)) {};

    SplittingSizeScheduleType schedule_type; // the splitting type
    int const ndists;                        // the number of districts
    int const total_seats;
    // Whether the plan level queries depend on more than the size tables. If
    // false `build_context` skips computing any per plan information.
    bool const plan_context_needed;
    double multidistrict_alpha;
    int const smallest_district_size;
    int const largest_district_size;
    std::vector<int> const district_seat_sizes;
    std::vector<bool> const is_district;
    bool const is_mmd;
    DistrictSizeSums const size_sums;

    /*
     * This is a vector of vectors where entry for index `r` it is a vector
     * of smaller cut sizes to try. ie (element, r-element)
     *      - For any split its just 1,...,floor(r/2)
     *      - for District splits its just 1
     */
    std::vector<std::vector<int>> all_regions_smaller_cut_sizes_to_try;

    // For each element its an array of the minimum and maximum possible cut sizes
    // of regions cut from
    std::vector<std::pair<int, int>> all_regions_min_and_max_possible_cut_sizes;

    // 1-indexed vector of if the region can be split
    std::vector<bool> valid_region_sizes_to_split;
    // 1-indexed vector of if that region size could have been the result of a valid split
    std::vector<bool> valid_split_region_sizes;

    // 1-indexed 2D matrix for if two regions (size1, size2) can be merged
    std::vector<std::vector<bool>> valid_merge_pair_sizes;

    // The size table entries for splitting a region of `region_size`
    CutSizeOptions size_cut_options(int const region_size) const {
        return {all_regions_smaller_cut_sizes_to_try[region_size],
                all_regions_min_and_max_possible_cut_sizes[region_size].first,
                all_regions_min_and_max_possible_cut_sizes[region_size].second};
    }

    /*
     * Plan level queries. Each takes the plan and the context built for it
     * by `build_context`.
     */
    // fills `ctx` with the per plan information this schedule needs
    void build_context(Plan const &plan, ScheduleContext &ctx) const;
    // whether `region_id` may be split
    bool region_is_splittable(Plan const &plan, ScheduleContext const &ctx,
                              int const region_id) const;
    // the allowed ways to split `region_id`
    CutSizeOptions region_cut_sizes(Plan const &plan, ScheduleContext &ctx,
                                    int const region_id) const;
    // the allowed ways to split the region made by merging the two regions
    CutSizeOptions merged_region_cut_sizes(Plan const &plan, ScheduleContext &ctx,
                                           int const region1_id, int const region2_id) const;
    // whether merging the two regions is a valid reversal of a split
    bool is_valid_merge_pair(Plan const &plan, ScheduleContext const &ctx,
                             int const region1_id, int const region2_id) const;
    // log probability the region made by merging the two regions would have
    // been picked to split, when regions are picked with weight size^alpha
    double log_merged_region_selection_prob(Plan const &plan, ScheduleContext const &ctx,
                                            int const region1_id, int const region2_id,
                                            double const selection_alpha) const;

    // This sets it for that split
    virtual void set_potential_cut_sizes_for_each_valid_size(int split_num,
                                                             int presplit_num_regions) = 0;
    // for district splits this allows to merge and split them during merge split
    virtual void update_cut_sizes_for_mergesplit_step(int split_num, int num_regions) {
        return;
    };

    // Throws if `plan`, a partial plan passed in to start sampling from, can't
    // be completed into a full plan under this schedule. Does nothing by default.
    virtual void check_initial_plan(Plan const &plan) const { return; };

    // Sets all the boolean vectors and matrices related to splitting and merging to false
    void reset_splitting_and_merge_booleans();

    void print_current_step_splitting_info();

  protected:
    /*
     * Implementations of the plan level queries. The defaults answer from
     * the size tables alone. Schedules where the answer depends on the rest
     * of the plan override these, and `fill_context` to cache whatever per
     * plan information they need.
     */
    virtual void fill_context(Plan const &plan, ScheduleContext &ctx) const { return; };
    virtual bool region_is_splittable_impl(Plan const &plan, ScheduleContext const &ctx,
                                           int const region_id) const;
    virtual CutSizeOptions region_cut_sizes_impl(Plan const &plan, ScheduleContext &ctx,
                                                 int const region_id) const;
    virtual CutSizeOptions merged_region_cut_sizes_impl(Plan const &plan, ScheduleContext &ctx,
                                                        int const region1_id,
                                                        int const region2_id) const;
    virtual bool is_valid_merge_pair_impl(Plan const &plan, ScheduleContext const &ctx,
                                          int const region1_id, int const region2_id) const;
    virtual double log_merged_region_selection_prob_impl(Plan const &plan,
                                                         ScheduleContext const &ctx,
                                                         int const region1_id,
                                                         int const region2_id,
                                                         double const selection_alpha) const;

  private:
    // throws if `ctx` wasn't built for a plan with `plan`'s number of regions
    void check_context(Plan const &plan, ScheduleContext const &ctx) const;
};

/*
 * Derived Class for one district split schedule
 *
 * This is the splitting schedule for splitting 1 district at a time
 */
class DistrictOnlySplittingSchedule : public SplittingSchedule {

  public:
    // constructor
    DistrictOnlySplittingSchedule(const int num_splits, const int ndists);

    void set_potential_cut_sizes_for_each_valid_size(int split_num,
                                                     int presplit_num_regions) override;
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
    AnyRegionSMDSplittingSchedule(const int num_splits, const int ndists);

    void set_potential_cut_sizes_for_each_valid_size(int split_num,
                                                     int presplit_num_regions) override;
};

/*
 * Derived Class for district only split schedule for multi-member districts
 * This is the splitting schedule where at each step we split off a district
 * from the remainder.
 *
 * With r regions every region but one (the remainder) must be a district and
 * the remainder must still be splittable into the remaining D - r + 1
 * districts. Which splits of the remainder are allowed only depends on the
 * step and its size, so the size tables handle that.
 *
 * Which region is the remainder is not always decided by size alone. If some
 * district size is a sum of other district sizes (eg sizes 3 to 6, 6 = 3 + 3)
 * a district sized region can also be splittable, so with a remainder that
 * isn't district sized, size alone would let that district be split and
 * leave a plan that can't be completed. In that case (`plan_context_needed`)
 * the context records the non district sized region, and when there is one
 * it is the only region that can be split.
 */
class DistrictOnlyMMDSplittingSchedule : public SplittingSchedule {

  public:
    // constructor
    DistrictOnlyMMDSplittingSchedule(const int num_splits, const int ndists,
                                     const int total_seats,
                                     std::vector<int> const &district_seat_sizes);

    void set_potential_cut_sizes_for_each_valid_size(int split_num,
                                                     int presplit_num_regions) override;
    void update_cut_sizes_for_mergesplit_step(int split_num, int num_regions) override;

    void check_initial_plan(Plan const &plan) const override;

  protected:
    void fill_context(Plan const &plan, ScheduleContext &ctx) const override;
    bool region_is_splittable_impl(Plan const &plan, ScheduleContext const &ctx,
                                   int const region_id) const override;
    CutSizeOptions merged_region_cut_sizes_impl(Plan const &plan, ScheduleContext &ctx,
                                                int const region1_id,
                                                int const region2_id) const override;
    bool is_valid_merge_pair_impl(Plan const &plan, ScheduleContext const &ctx,
                                  int const region1_id, int const region2_id) const override;
    double log_merged_region_selection_prob_impl(Plan const &plan, ScheduleContext const &ctx,
                                                 int const region1_id, int const region2_id,
                                                 double const selection_alpha) const override;

  private:
    /*
     * During merge split rounds two adjacent districts can also be merged and
     * resplit into two districts, as long as the plan has a non district sized
     * remainder that isn't one of them. That leaves the remainder unchanged so
     * the plan stays completable. (When every region is district sized this
     * isn't always safe, so it isn't allowed.)
     *
     * `mergesplit_district_pairs` is turned on by
     * `update_cut_sizes_for_mergesplit_step` and off again by the next
     * `set_potential_cut_sizes_for_each_valid_size`, so SMC steps never see it.
     */
    bool mergesplit_district_pairs = false;
    // For each merged size, the ways to split it into two districts
    std::vector<std::vector<int>> district_pair_smaller_cut_sizes;
    std::vector<std::pair<int, int>> district_pair_min_and_max_cut_sizes;

    // whether merging the two regions is a district pair merge as described above
    bool is_district_pair_merge(Plan const &plan, ScheduleContext const &ctx,
                                int const region1_id, int const region2_id) const;
};

/*
 * Derived Class for any region split schedule for multi-member districts
 * This is the splitting schedule for which allows for any size splits
 * to be made (where both split regions have to have at least size district_size_lb)
 */
class AnyRegionMMDSplittingSchedule : public SplittingSchedule {

  public:
    // constructor
    AnyRegionMMDSplittingSchedule(const int num_splits, const int ndists, const int total_seats,
                                  std::vector<int> const &district_seat_sizes);

    std::vector<bool> reachable_size;

    void set_potential_cut_sizes_for_each_valid_size(int split_num,
                                                     int presplit_num_regions) override;
};

/*
 * Derived Class for pure merge split. So you can merge any valid district
 * sizes and create any district sizes.
 */
class PureMSSplittingSchedule : public SplittingSchedule {

  public:
    // constructor
    PureMSSplittingSchedule(const int ndists, const int total_size,
                            std::vector<int> const &district_seat_sizes);

    void set_potential_cut_sizes_for_each_valid_size(int split_num,
                                                     int presplit_num_regions) override {
        throw std::runtime_error(
            "PureMSSplittingSchedule::set_potential_cut_sizes_for_each_valid_size was called, "
            "but a pure mergesplit schedule never splits plans with SMC so has no cut sizes.");
    };
};


std::unique_ptr<SplittingSchedule>
get_splitting_schedule(const int num_splits, const int ndists, const int total_seats,
                       std::vector<int> const &district_seat_sizes,
                       SplittingSizeScheduleType const schedule_type);

#endif
