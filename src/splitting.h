#pragma once
#ifndef SPLITTING_H
#define SPLITTING_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <iostream>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "redist_types.h"


//' Selects a multidistrict with probability proportional to its d_nk value and
//' returns the log probability of the selected region
//'
//' Given a plan object with at least one multidistrict this function randomly
//' selects a multidistrict with probability proporitional to its d_nk value
//' (relative to all multidistricts) and returns the log of the probability that
//' region was chosen.
//'
//'
//' @title Choose multidistrict to split
//'
//' @param plan A plan object
//' @param region_to_split an integer that will be updated by reference with the
//' id number of the region selected to split
//'
//' @details No modifications to inputs made
//'
//' @return the region level graph
//'
double choose_multidistrict_to_split(
        Plan const&plan, int &region_id_to_split);



//' Attempts to cut one region into two from a spanning tree and if successful
//' returns information on what the two new regions would be. Does not actually
//' update the plan
//'
//' Takes a spanning tree `ust` drawn on a specific region and attempts to cut
//' it to produce two new regions using the generalized splitting procedure
//' outlined <PAPER HERE>. This function is based on `cut_districts` in `smc.cpp`
//' however the crucial difference is even if a cut is successful it does not
//' update the plan. Instead it just returns the information on the two new
//' regions if successful and the vertices to use to update the plans.
//'
//' Depending on the value of split_district_only will only attempt to split off
//' a single district or allows for more general splits.
//'
//' By convention the first new region (`new_region1`) will always be the region
//' with the smaller d-value (although they can be equal).
//'
//' @title Attempt to Cut Region Tree into Two New Regions
//'
//' @param ust A directed spanning tree passed by reference
//' @param root The root vertex of the spanning tree
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param split_district_only If true then only tries to split a district, if false allows for
//' arbitrary region splits
//' @param pop A vector of the population associated with each vertex in `g`
//' @param plan A plan object
//' @param region_id_to_split The id of the region in the plan object we're attempting to split
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param new_region1_tree_root The index of the root of tree associated with
//' the first new region (if the tree cut was successful)
//' @param new_region1_dval The d-value of the first new region (if the tree cut
//'  was successful)
//' @param new_region1_pop The population of the first new region (if the tree cut
//' was successful)
//' @param new_region2_tree_root The index of the root of tree associated with
//' the second new region (if the tree cut was successful)
//' @param new_region2_dval The d-value of the second new region (if the tree cut
//'  was successful)
//' @param new_region2_pop The population of the second new region (if the tree cut
//' was successful)
//'
//' @details Modifications
//'    - If two new valid regions are split then the tree `ust` is cut into two
//'    distjoint pieces
//'    - If two new valid regions are split then the 6 `new_region` inputs are all
//'    updated by reference with the values associated with the new regions
//'
//' @return True if two valid regions were successfully split, false otherwise
//'
bool get_edge_to_cut(Tree &ust, int root,
                     int k_param, bool split_district_only,
                     const uvec &pop, const Plan &plan, const int region_id_to_split,
                     const double lower, const double upper, const double target,
                     int &new_region1_tree_root, int &new_region1_dval, double &new_region1_pop,
                     int &new_region2_tree_root, int &new_region2_dval, double &new_region2_pop
);



//' Updates a `Plan` object using a cut tree
//'
//' Takes a cut spanning tree `ust` and variables on the two new regions
//' induced by the cuts and updates `plan` to add those two new regions.
//'
//'
//' @title Update plan regions from cut tree
//'
//' @param ust A cut (ie has two partition pieces) directed spanning tree
//' passed by reference
//' @param plan A plan object
//' @param old_region_id The id of old (split) region
//' @param new_region1_tree_root The vertex of the root of one piece of the cut
//' tree. This always corresponds to the region with the smaller dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region1_dval The dval associated with the new region 1
//' @param new_region1_pop The population associated with the new region 1
//' @param new_region2_tree_root The vertex of the root of other piece of the cut
//' tree. This always corresponds to the region with the bigger dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region2_dval The dval associated with the new region 2
//' @param new_region2_pop The population associated with the new region 2
//' @param new_region1_id The id the new region 1 was assigned in the plan
//' @param new_region2_id The id the new region 2 was assigned in the plan
//'
//' @details Modifications
//'    - `plan` is updated in place with the two new regions and the old region
//'    is removed
//'    - `new_region1_id` and `new_region2_id` are updated by reference to what
//'    the values of the two new region ids were set to
//'
void update_plan_from_cut(
        Tree &ust, Plan &plan,
        const int old_region_id,
        const int new_region1_tree_root, const int new_region1_dval, const double new_region1_pop,
        const int new_region2_tree_root, const int new_region2_dval, const double new_region2_pop,
        int &new_region1_id,  int &new_region2_id
);



//' Splits a multidistrict in all of the plans
//'
//' Using the procedure outlined in <PAPER HERE> this function attempts to split
//' a multidistrict in a previous steps plan until M successful splits have been made. This
//' is based on the `split_maps` function in smc.cpp
//'
//' @title Split all the maps
//'
//' @param g A graph (adjacency list) passed by reference
//' @param counties Vector of county labels of each vertex in `g`
//' @param cg County level multigraph
//' @param pop A vector of the population associated with each vertex in `g`
//' @param old_plans_vec A vector of plans from the previous step
//' @param new_plans_vec A vector which will be filled with plans that had a
//' multidistrict split to make them
//' @param original_ancestor_vec A vector used to track which original ancestor
//' the new plans descended from. The value  of `original_ancestor_vec[i]`
//' is the index of the original ancestor the new plan `new_plans_vec[i]` is
//' descended from.
//' @param parent_vec A vector used to track the index of the previous plan
//' sampled that was successfully split. The value of `parent_vec[i]` is the
//' index of the old plan from which the new plan `new_plans_vec[i]` was
//' successfully split from. In other words `new_plans_vec[i]` is equal to
//' `attempt_region_split(old_plans_vec[parent_vec[i]], ...)`
//' @param prev_ancestor_vec A vector used to track the index of the original
//' ancestor of the previous plans. The value of `prev_ancestor_vec[i]` is the
//' index of the original ancestor of `old_plans_vec[i]`
//' @param unnormalized_sampling_weights A vector of weights used to sample indices
//' of the `old_plans_vec`. The value of `unnormalized_sampling_weights[i]` is
//' the unnormalized probability that index i is selected
//' @param normalized_weights_to_fill_in A vector which will be filled with the
//' normalized weights the index sampler uses. The value of
//' `normalized_weights_to_fill_in[i]` is the probability that index i is selected
//' @param draw_tries_vec A vector used to keep track of how many plan split
//' attempts were made for index i. The value `draw_tries_vec[i]` represents how
//' many split attempts were made for the i-th new plan (including the successful
//' split). For example, `draw_tries_vec[i] = 1` means that the first split
//' attempt was successful.
//' @param parent_unsuccessful_tries_vec A vector used to keep track of how many times the
//' previous rounds plans were sampled and unsuccessfully split. The value
//' `parent_unsuccessful_tries_vec[i]` represents how many times `old_plans_vec[i]` was sampled
//' and then unsuccessfully split while creating all `M` of the new plans.
//' THIS MAY NOT BE THREAD SAFE
//' @param accept_rate The number of accepted splits over the total number of
//' attempted splits. This is equal to `sum(draw_tries_vec)/M`
//' @param n_unique_parent_indices The number of unique parent indices, ie the
//' number of previous plans that had at least one descendant amongst the new
//' plans. This is equal to `unique(parent_vec)`
//' @param n_unique_original_ancestors The number of unique original ancestors,
//' in the new plans. This is equal to `unique(original_ancestor_vec)`
//' @param ancestors Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
//' WHAT IT IS DOING
//' @param lags Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
//' WHAT IT IS DOING
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param k_param The top edges to pick parameter for the region splitting
//' algorithm
//' @param split_district_only Whether or not to only allow for single district
//' splits. If set to `true` will only attempt to split off one district at a
//' time
//' @param pool A threadpool for multithreading
//' @param verbosity A parameter controlling the amount of detail printed out
//' during the algorithms running
//'
//' @details Modifications
//'    - The `new_plans_vec` is updated with all the newly split plans
//'    - The `old_plans_vec` is updated with all the newly split plans as well.
//'    Note that the reason both this and `new_plans_vec` are updated is because
//'    of the nature of the code you need both vectors and so both are passed by
//'    reference to save memory.
//'    - The `original_ancestor_vec` is updated to contain the indices of the
//'    original ancestors of the new plans
//'    - The `parent_vec` is updated to contain the indices of the parents of the
 //'    new plans
//'    - If two new valid regions are split then the new_region_ids is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'    - The `normalized_weights_to_fill_in` is updated to contain the normalized
//'    probabilities the index sampler used. This is only collected for diagnostics
//'    at this point and should just be equal to `unnormalized_sampling_weights`
//'    divided by `sum(unnormalized_sampling_weights)`
//'    - The `draw_tries_vec` is updated to contain the number of tries for each
//'    of the new plans
//'    - The `parent_unsuccessful_tries_vec` is updated to contain the number of unsuccessful
//'    samples of the old plans
//'    - The `accept_rate` is updated to contain the average acceptance rate for
//'    this iteration
//'    - `n_unique_parent_indices` and `n_unique_original_ancestors` are updated
//'    with the unique number of parents and original ancestors for all the new
//'    plans respectively
//'    - `ancestors` is updated to something. THIS IS FROM ORIGINAL SMC CODE,
//'    I DO NOT KNOW WHAT IT MEANS
//'
//' @return nothing
//'
void generalized_split_maps(
        const Graph &g, const uvec &counties, Multigraph &cg, const uvec &pop,
        std::vector<Plan> &old_plans_vec, std::vector<Plan> &new_plans_vec,
        std::vector<int> &original_ancestor_vec,
        std::vector<int> &parent_vec,
        const std::vector<int> &prev_ancestor_vec,
        const std::vector<double> &unnormalized_sampling_weights,
        std::vector<double> &normalized_weights_to_fill_in,
        std::vector<int> &draw_tries_vec,
        std::vector<int> &parent_unsuccessful_tries_vec,
        double &accept_rate,
        int &n_unique_parent_indices,
        int &n_unique_original_ancestors,
        umat &ancestors, const std::vector<int> &lags,
        double lower, double upper, double target,
        int k_param, bool split_district_only,
        RcppThread::ThreadPool &pool,
        int verbosity
);


#endif
