/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/02
 * Purpose: Implement constraint and scoring function
 ********************************************************/

 #include "scoring.h"


 /********************************************************
 * Constraints Supported
 *   - Population tempering (PopTemperConstraint)
 *      This one enoucrages smaller deviations 
 *   - Group Hinge (GroupHingeConstraint)
 *  
 ********************************************************/



 double compute_log_pop_temper(
    double const target, double const pop_temper, int const ndists,
    int const region_pop, int const region_size
){
    double region_target = target*region_size;
    // get population deviation 
    double const pop_dev = std::fabs(
        static_cast<double>(region_pop) - region_target
    )/region_target;

    double const pop_pen = std::sqrt(static_cast<double>(ndists) - 2) * std::log(1e-12 + pop_dev);

    // now return the values for the old region minus the two new ones
    return pop_pen * pop_temper;
}


double PopTemperConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    int const region_pop = plan.region_pops[region_id];
    int const region_size = plan.region_sizes[region_id];
    return compute_log_pop_temper(target, pop_temper, ndists, region_pop, region_size);
}


double PopTemperConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    int const merged_region_pop = plan.region_pops[region1_id] + plan.region_pops[region2_id];
    int const merged_region_size = plan.region_sizes[region1_id] + plan.region_sizes[region2_id];
    return compute_log_pop_temper(target, pop_temper, ndists, merged_region_pop, merged_region_size);
}



/*
 * Compute the new, hinge group penalty for district `distr`
 *
 * MODIFIED FROM eval_grp_hinge to just take idx instead of distr
 * idxs is  uvec idxs = find(region_ids == region_id);
 */
double eval_grp_hinge_gsmc_version(
    PlanVector const &region_ids, 
    const arma::vec &tgts_grp, const arma::uvec &grp_pop, const arma::uvec &total_pop,
    RegionID const &region1_id, RegionID const &possible_region2_id = -42) {
    RegionID const region2_id = possible_region2_id < 0 ? region1_id : possible_region2_id;
    double subsetted_grp_pop_sum = 0.0;
    double subsetted_total_pop_sum = 0.0;
    // get the sum of the two columns in region 1 or 2
    for (size_t i = 0; i < region_ids.size(); i++)
    {
        auto const region_id = region_ids[i];
        if(region_id == region1_id || region_id == region2_id){
            subsetted_grp_pop_sum += grp_pop(i);
            subsetted_total_pop_sum += total_pop(i);
        }
    }
    // do subsetted_grp_pop_sum/subsetted_total_pop_sum
    double frac = std::exp(
        std::log(subsetted_grp_pop_sum) - std::log(subsetted_total_pop_sum)
    );
    // figure out which to compare it to
    double target;
    double diff = 1;
    int n_tgt = tgts_grp.size();
    for (int i = 0; i < n_tgt; i++) {
        double new_diff = std::fabs(tgts_grp[i] - frac);
        if (new_diff <= diff) {
            diff = new_diff;
            target = tgts_grp[i];
        }
    }

    return std::sqrt(std::max(0.0, target - frac));
}



double GroupHingeConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_grp_hinge_gsmc_version(
        plan.region_ids, 
        tgts_group, group_pop, total_pop,
        region_id);

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double GroupHingeConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    double raw_score = eval_grp_hinge_gsmc_version(
        plan.region_ids,  
        tgts_group, group_pop, total_pop,
        region1_id, region2_id);
    return strength * raw_score;
}



std::pair<bool, double> ValidDistrictsConstraint::compute_plan_constraint_score(const Plan &plan) const{
    // first check no region is smaller than the smallest district size 
    for (size_t i = 0; i < plan.num_regions; i++)
    {
        auto region_size = plan.region_sizes[i];
        if(region_size < map_params.smallest_district_size){
            return std::make_pair(false, 0);
        }
        // if the plan is all districts, ie ndists = num regions then check every region is a district
        if(plan.num_regions == map_params.ndists && !map_params.is_district[region_size]){
            // if not a district return false
            return std::make_pair(false, 0);
        }
    }
    // Now return ok 
    return std::make_pair(true, 0);
    
    
}


// Scoring function
ScoringFunction::ScoringFunction(
    MapParams const &map_params,
    Rcpp::List const &constraints, double const pop_temper
):
map_params(map_params), num_non_final_soft_constraints(0), num_final_soft_constraints(0), all_rounds_soft_constraints(0), 
total_soft_constraints(0), num_hard_constraints(0){
    // First add region constraints 
    // add pop temper if doing that 
    if(pop_temper != 0){
        non_final_region_constraint_ptrs.emplace_back(
            std::make_unique<PopTemperConstraint>(
                map_params.target, map_params.ndists, 
                map_params, pop_temper
            ));
        num_non_final_soft_constraints++;
    }

    // add constraints if in the list 
    if (constraints.containsElementNamed("pop_dev")) {
        Rcpp::stop("Error: 'pop_dev' constraint not implemented yet!");
        Rcpp::List constr = constraints["pop_dev"];

        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            // check if the 
            if (strength != 0) {
                // Add constraint here 
            }
        }
        // create here
        // lp[i] += add_constraint("pop_dev", constraints,
        //     [&] (List l) -> double {
        //         return eval_pop_dev(districts.col(i), j,
        //                                pop, parity);
        //     });
        // all_rounds_soft_constraints++;
    }
    if (constraints.containsElementNamed("status_quo")) {
        throw Rcpp::exception("Error: 'status_quo' constraint not implemented yet!");
        // create here
        // lp[i] += add_constraint("status_quo", constraints,
        // [&] (List l) -> double {
        // return eval_sq_entropy(districts.col(i), as<uvec>(l["current"]),
        //          j, pop, n_distr,
        //          as<int>(l["n_current"]), V);
        // });
        // all_rounds_soft_constraints++;
    }
    if (constraints.containsElementNamed("segregation")) {
        throw Rcpp::exception("Error: 'segregation' constraint not implemented yet!");
        // create here
        // lp[i] += add_constraint("segregation", constraints,
        //     [&] (List l) -> double {
        //         return eval_segregation(districts.col(i), j,
        //                                 as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]));
        //     });
        // all_rounds_soft_constraints++;
    }
    if (constraints.containsElementNamed("grp_pow")) {
        throw Rcpp::exception("Error: 'grp_pow' constraint not implemented yet!");
        // create here
        // [&] (List l) -> double {
        // return eval_grp_pow(districts.col(i), j,
        //       as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]),
        //       as<double>(l["tgt_group"]), as<double>(l["tgt_other"]),
        //       as<double>(l["pow"]));
        // });
        // all_rounds_soft_constraints++;
    }
    if (constraints.containsElementNamed("compet")) {
        throw Rcpp::exception("Error: 'compet' constraint not implemented yet!");
        // lp[i] += add_constraint("compet", constraints,
        // [&] (List l) -> double {
        // uvec dvote = l["dvote"];
        // uvec total = dvote + as<uvec>(l["rvote"]);
        // return eval_grp_pow(districts.col(i), j,
        //       dvote, total, 0.5, 0.5, as<double>(l["pow"]));
        // });
        // all_rounds_soft_constraints++;
    }
    if (constraints.containsElementNamed("grp_hinge")) {
        // create constraint objects
        Rcpp::List constr = constraints["grp_hinge"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<GroupHingeConstraint>(
                        strength, as<arma::vec>(constr_inst["tgts_group"]),
                        as<arma::uvec>(constr_inst["group_pop"]), as<arma::uvec>(constr_inst["total_pop"]),
                        constr_score_districts_only
                    ));
                // REprintf("Added and its %d!\n", region_constraint_ptrs[0]->score_districts_only);
                // increase constraints applied at every split count
                all_rounds_soft_constraints++;
            }
        }
    }
    if (constraints.containsElementNamed("grp_inv_hinge")) {
        // NOTE: As far as I can tell grp_inv_hinge is same code as grp_inv_hinge
        // create constraint objects
        Rcpp::List constr = constraints["grp_inv_hinge"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<GroupHingeConstraint>(
                        strength, as<arma::vec>(constr_inst["tgts_group"]),
                        as<arma::uvec>(constr_inst["group_pop"]), as<arma::uvec>(constr_inst["total_pop"]),
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_constraints++;
            }
        }
    }

    // Now add plan constraints 
    if (constraints.containsElementNamed("plan_valid_district_sizes")) {
        // Add this check 
        hard_plan_constraint_ptrs.emplace_back(
            std::make_unique<ValidDistrictsConstraint>(map_params)
        );
    }

    total_soft_constraints = num_non_final_soft_constraints+num_final_soft_constraints+all_rounds_soft_constraints;
    any_soft_constraints = total_soft_constraints != 0; 

}




// lp[i] += add_constraint("incumbency", constraints,
// [&] (List l) -> double {
// return eval_inc(districts.col(i), j, as<uvec>(l["incumbents"]));
// });

// lp[i] += add_constraint("splits", constraints,
// [&] (List l) -> double {
// return eval_splits(districts.col(i), j, as<uvec>(l["admin"]), l["n"], true);
// });

// lp[i] += add_constraint("multisplits", constraints,
// [&] (List l) -> double {
// return eval_multisplits(districts.col(i), j, as<uvec>(l["admin"]), l["n"], true);
// });

// lp[i] += add_constraint("total_splits", constraints,
// [&] (List l) -> double {
// return eval_total_splits(districts.col(i), j, as<uvec>(l["admin"]), l["n"], true);
// });

// lp[i] += add_constraint("polsby", constraints,
//     [&] (List l) -> double {
//         return eval_polsby(districts.col(i), j,
//                            as<ivec>(l["from"]),
//                            as<ivec>(l["to"]), as<vec>(l["area"]),
//                            as<vec>(l["perimeter"]));
//     });

// lp[i] += add_constraint("fry_hold", constraints,
//     [&] (List l) -> double {
//         return eval_fry_hold(districts.col(i), j,
//                              as<uvec>(l["total_pop"]),
//                              as<mat>(l["ssdmat"]),
//                              as<double>(l["denominator"]));
//     });

// lp[i] += add_constraint("qps", constraints,
//     [&] (List l) -> double {
//         return eval_qps(districts.col(i), j,
//                         as<uvec>(l["total_pop"]),
//                         as<uvec>(l["cities"]), as<int>(l["n_city"]),
//                         n_distr);
//     });

// lp[i] += add_constraint("custom", constraints,
// [&] (List l) -> double {
// Function fn = l["fn"];
// return as<NumericVector>(fn(districts.col(i), j))[0];
// });



double ScoringFunction::compute_region_score(const Plan &plan, int const region_id, bool const is_final) 
    const{    
    // start out with score of zero
    double region_score = 0.0;

    // check if its a multidistrict 
    bool const is_multidistrict = !map_params.district_seat_sizes[plan.region_sizes[region_id]];

    // add the score from each constraint 
    for(auto const &constraint_ptr: region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;
        region_score += constraint_ptr->compute_region_constraint_score(plan, region_id);
    }
    // if final then return 
    if(is_final) return region_score;

    // if not final then add those constraints (for now just pop_temper)
    for(auto const &constraint_ptr: non_final_region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;
        region_score += constraint_ptr->compute_region_constraint_score(plan, region_id);
    }

    return region_score;
}


double ScoringFunction::compute_merged_region_score(const Plan &plan, 
    int const region1_id, int const region2_id, bool const is_final) 
    const{
    // start out with log score of zero
    double region_score = 0.0;

    // check if its a multidistrict 
    bool const is_multidistrict = !map_params.district_seat_sizes[
        plan.region_sizes[region1_id] + plan.region_sizes[region2_id]
    ];

    // get the log constraint 
    for(auto const &constraint_ptr: region_constraint_ptrs){
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;
        region_score += constraint_ptr->compute_merged_region_constraint_score(plan, region1_id, region2_id);
    }
    // if final then return 
    if(is_final) return region_score;

    // if not final then add those constraints 
    for(auto const &constraint_ptr: non_final_region_constraint_ptrs){
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;
        region_score += constraint_ptr->compute_merged_region_constraint_score(plan, region1_id, region2_id);
    }

    return region_score;
}

std::pair<bool, double> ScoringFunction::compute_hard_plan_constraints_score(const Plan &plan) const{
    double plan_score = 0.0;
 
    for(auto const &constraint_ptr: hard_plan_constraint_ptrs){
        auto result = constraint_ptr->compute_plan_constraint_score(plan);
        // if false then failed to satisfy constraint so return false
        if(!result.first){
            return std::make_pair(false, 0.0);
        }else{
            plan_score += result.second;
        }
    }

    return std::make_pair(true, plan_score);
}