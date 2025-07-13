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
 *   - Incumbency (IncumbentConstraint)
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



double PopDevConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_pop_dev(
        plan.region_ids, 
        region_id, region_id,
        total_pop, parity
    );

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double PopDevConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    double raw_score = eval_pop_dev(
        plan.region_ids, 
        region1_id, region2_id,
        total_pop, parity
    );
    return strength * raw_score;
}




double GroupHingeConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_grp_hinge(
        plan.region_ids, V, region_id, region_id,
        tgts_group, group_pop, total_pop
    );

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double GroupHingeConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    double raw_score = eval_grp_hinge(
        plan.region_ids, V, region1_id, region2_id,
        tgts_group, group_pop, total_pop
    );
    return strength * raw_score;
}



double IncumbentConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_inc(
        plan.region_ids, 
        region_id, region_id,
        incumbents);

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double IncumbentConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    double raw_score = eval_inc(
        plan.region_ids, 
        region1_id, region2_id,
        incumbents);
    return strength * raw_score;
}


double StatusQuoConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_sq_entropy(
        plan.region_ids, current,
        region_id, region_id,
        pop, 
        ndists, n_current, V
    );

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double StatusQuoConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    double raw_score = eval_sq_entropy(
        plan.region_ids, current,
        region1_id, region2_id,
        pop, 
        ndists, n_current, V
    );
    return strength * raw_score;
}


double SplitsConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_splits(
        plan.region_ids, region_id,
        admin_units, n_admin_units, smc
    );

    
    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double SplitsConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    // make a copy of the plans 
    std::vector<RegionID> dummy_merged_vec(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );
    // now make all instance of region2 region1
    for (size_t i = 0; i < dummy_merged_vec.size(); i++)
    {
        if(dummy_merged_vec[i] == region2_id){
            dummy_merged_vec[i] = region1_id;
        }
    }
    // very crude, if a bottleneck can probably just make eval splits support
    // implicit merging 
    double raw_score = eval_splits(
        dummy_merged_vec, region1_id,
        admin_units, n_admin_units, smc
    );
    return strength * raw_score;
}

double MultisplitsConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_multisplits(
        plan.region_ids, region_id,
        admin_units, n_admin_units, smc
    );

    
    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double MultisplitsConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    // make a copy of the plans 
    std::vector<RegionID> dummy_merged_vec(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );
    // now make all instance of region2 region1
    for (size_t i = 0; i < dummy_merged_vec.size(); i++)
    {
        if(dummy_merged_vec[i] == region2_id){
            dummy_merged_vec[i] = region1_id;
        }
    }
    // very crude, if a bottleneck can probably just make eval splits support
    // implicit merging 
    double raw_score = eval_multisplits(
        dummy_merged_vec, region1_id,
        admin_units, n_admin_units, smc
    );
    return strength * raw_score;
}


double TotalSplitsConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_total_splits(
        plan.region_ids, region_id,
        admin_units, n_admin_units, smc
    );

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double TotalSplitsConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    // make a copy of the plans 
    std::vector<RegionID> dummy_merged_vec(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );
    // now make all instance of region2 region1
    for (size_t i = 0; i < dummy_merged_vec.size(); i++)
    {
        if(dummy_merged_vec[i] == region2_id){
            dummy_merged_vec[i] = region1_id;
        }
    }
    // very crude, if a bottleneck can probably just make eval splits support
    // implicit merging 
    double raw_score = eval_total_splits(
        dummy_merged_vec, region1_id,
        admin_units, n_admin_units, smc
    );
    return strength * raw_score;
}


double PolsbyConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    double raw_score = eval_polsby(
        plan.region_ids, region_id, region_id, V,
        from, to, area, perimeter
    );

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double PolsbyConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    double raw_score = eval_polsby(
        plan.region_ids, region1_id, region2_id, V,
        from, to, area, perimeter
    );

    return strength * raw_score;
}


double CustomRegionConstraint::compute_region_constraint_score(const Plan &plan, int const region_id) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );

    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, region_id))[0]
    );

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double CustomRegionConstraint::compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );

    // now make all instance of region2 region1
    for (size_t i = 0; i < plan.region_ids.size(); i++)
    {
        if(rcpp_plan_wrap[i] == region2_id){
            rcpp_plan_wrap[i] = region1_id;
        }
    }
    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, region1_id))[0]
    );

    return strength * raw_score;
}

double CustomPlanConstraint::compute_plan_constraint_score(const Plan &plan) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );

    Rcpp::IntegerVector rcpp_sizes_wrap(
        plan.region_sizes.begin(),
        plan.region_sizes.end()
    );
    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, rcpp_sizes_wrap))[0]
    );

    return strength * raw_score;
}
// log constraint for region made by merging region 1 and 2
double CustomPlanConstraint::compute_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );

    Rcpp::IntegerVector rcpp_sizes_wrap(
        plan.region_sizes.begin(),
        plan.region_sizes.end()
    );

    // now make all instance of region2 region1
    for (size_t i = 0; i < plan.region_ids.size(); i++)
    {
        if(rcpp_plan_wrap[i] == region2_id){
            rcpp_plan_wrap[i] = region1_id;
        }
    }

    // adjust the sizes as well so region 1 has all the sizes
    rcpp_sizes_wrap[region1_id] += rcpp_sizes_wrap[region2_id];
    rcpp_sizes_wrap[region2_id] = 0;
    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, rcpp_sizes_wrap))[0]
    );

    return strength * raw_score;
}


std::pair<bool, double> CustomHardPlanConstraint::compute_plan_constraint_score(const Plan &plan) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );

    Rcpp::IntegerVector rcpp_sizes_wrap(
        plan.region_sizes.begin(),
        plan.region_sizes.end()
    );

    bool status = as<bool>(fn(rcpp_plan_wrap, rcpp_sizes_wrap));
    
    // Now return result
    return std::make_pair(status, 0);
}



std::pair<bool, double> CustomHardPlanConstraint::compute_merged_plan_constraint_score(
            const Plan &plan, int const region1_id, int const region2_id
) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );

    Rcpp::IntegerVector rcpp_sizes_wrap(
        plan.region_sizes.begin(),
        plan.region_sizes.end()
    );

    // now make all instance of region2 region1
    for (size_t i = 0; i < plan.region_ids.size(); i++)
    {
        if(rcpp_plan_wrap[i] == region2_id){
            rcpp_plan_wrap[i] = region1_id;
        }
    }

    // adjust the sizes as well so region 1 has all the sizes
    rcpp_sizes_wrap[region1_id] += rcpp_sizes_wrap[region2_id];
    rcpp_sizes_wrap[region2_id] = 0;

    bool status = as<bool>(fn(rcpp_plan_wrap, rcpp_sizes_wrap));
    
    // Now return result
    return std::make_pair(status, 0);
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
    Rcpp::List const &constraints, double const pop_temper, bool const smc
):
map_params(map_params), num_non_final_soft_region_constraints(0), num_final_soft_region_constraints(0), all_rounds_soft_region_constraints(0), 
num_non_final_soft_plan_constraints(0), num_final_soft_plan_constraints(0), all_rounds_soft_plan_constraints(0), 
num_hard_plan_constraints(0),
total_soft_constraints(0),
any_soft_custom_constraints(false), any_hard_custom_constraints(false){
    // First add region constraints 
    // add pop temper if doing that 
    if(pop_temper != 0){
        non_final_region_constraint_ptrs.emplace_back(
            std::make_unique<PopTemperConstraint>(
                map_params.target, map_params.ndists, 
                map_params, pop_temper
            ));
        num_non_final_soft_region_constraints++;
    }

    // Add region constraints 
    // add constraints if in the list 
    if (constraints.containsElementNamed("pop_dev")) {
        Rcpp::List constr = constraints["pop_dev"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<PopDevConstraint>(
                        strength, 
                        map_params.target, map_params.pop,
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
    }
    if (constraints.containsElementNamed("status_quo")) {
        Rcpp::List constr = constraints["status_quo"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<StatusQuoConstraint>(
                        strength, 
                        as<arma::uvec>(constr_inst["current"]), map_params.pop,
                        map_params.ndists, as<int>(constr_inst["n_current"]), map_params.V,
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
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
                        strength, map_params.V, as<arma::vec>(constr_inst["tgts_group"]),
                        as<arma::uvec>(constr_inst["group_pop"]), as<arma::uvec>(constr_inst["total_pop"]),
                        constr_score_districts_only
                    ));
                // REprintf("Added and its %d!\n", region_constraint_ptrs[0]->score_districts_only);
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
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
                        strength, map_params.V, as<arma::vec>(constr_inst["tgts_group"]),
                        as<arma::uvec>(constr_inst["group_pop"]), as<arma::uvec>(constr_inst["total_pop"]),
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
    }
    if (constraints.containsElementNamed("incumbency")) {
        Rcpp::List constr = constraints["incumbency"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<IncumbentConstraint>(
                        strength, as<arma::uvec>(constr_inst["incumbents"]),
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
    }
    if (constraints.containsElementNamed("splits")) {
        Rcpp::List constr = constraints["splits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<SplitsConstraint>(
                        strength, 
                        as<arma::uvec>(constr_inst["admin"]), as<int>(constr_inst["n"]),
                        smc,
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
    }
    if (constraints.containsElementNamed("multisplits")) {
        Rcpp::List constr = constraints["multisplits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<MultisplitsConstraint>(
                        strength, 
                        as<arma::uvec>(constr_inst["admin"]), as<int>(constr_inst["n"]),
                        smc,
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
    }
    if (constraints.containsElementNamed("total_splits")) {
        Rcpp::List constr = constraints["total_splits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<TotalSplitsConstraint>(
                        strength, 
                        as<arma::uvec>(constr_inst["admin"]), as<int>(constr_inst["n"]),
                        smc,
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
    }
    if (constraints.containsElementNamed("polsby")) {
        Rcpp::List constr = constraints["polsby"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                region_constraint_ptrs.emplace_back(
                    std::make_unique<PolsbyConstraint>(
                        strength, map_params.V,
                        as<arma::ivec>(constr_inst["from"]), 
                        as<arma::ivec>(constr_inst["to"]),
                        as<arma::vec>(constr_inst["area"]),
                        as<arma::vec>(constr_inst["perimeter"]),
                        constr_score_districts_only
                    ));
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
            }
        }
    }

    if (constraints.containsElementNamed("custom")) {
        Rcpp::List constr = constraints["custom"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                bool constr_score_districts_only = false;
                if (constr_inst.containsElementNamed("score_districts_only")){
                    constr_score_districts_only = as<bool>(constr_inst["score_districts_only"]);
                }
                // Function fn = constr_inst["fn"];
                region_constraint_ptrs.emplace_back(
                    std::make_unique<CustomRegionConstraint>(
                        strength, 
                        as<Rcpp::Function>(constr_inst["fn"]), 
                        constr_score_districts_only
                    )
                );
                // increase constraints applied at every split count
                all_rounds_soft_region_constraints++;
                // mark custom R constraints as true 
                any_soft_custom_constraints = true; 
            }
        }
    }

    // Now add plan constraints 
    // Add soft plan constraints 
    if (constraints.containsElementNamed("custom_plan")) {
        Rcpp::List constr = constraints["custom_plan"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            if (strength != 0) {
                plan_constraint_ptrs.emplace_back(
                    std::make_unique<CustomPlanConstraint>(
                        strength, 
                        as<Rcpp::Function>(constr_inst["fn"])
                    )
                );
                // increase constraints applied at every split count
                all_rounds_soft_plan_constraints++;
                any_soft_custom_constraints = true; 
            }
        }
    }

    // Now add hard constraints 
    if (constraints.containsElementNamed("custom_hard_plan")) {
        Rcpp::List constr = constraints["custom_hard_plan"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = 1;
            if (strength != 0) {
                hard_plan_constraint_ptrs.emplace_back(
                    std::make_unique<CustomHardPlanConstraint>(
                        as<Rcpp::Function>(constr_inst["fn"])
                    )
                );
                // increase constraints applied at every split count
                num_hard_plan_constraints++;
                any_hard_custom_constraints = true; 
            }
        }
    }

    if (constraints.containsElementNamed("plan_valid_district_sizes")) {
        // Add this check 
        hard_plan_constraint_ptrs.emplace_back(
            std::make_unique<ValidDistrictsConstraint>(map_params)
        );
        num_hard_plan_constraints++;
    }

    total_soft_region_constraints = num_non_final_soft_region_constraints+num_final_soft_region_constraints+all_rounds_soft_region_constraints;
    total_soft_plan_constraints = num_non_final_soft_plan_constraints + num_final_soft_plan_constraints + all_rounds_soft_plan_constraints;

    total_soft_constraints = total_soft_region_constraints + total_soft_plan_constraints;

    any_soft_region_constraints = total_soft_region_constraints > 0; 
    any_soft_plan_constraints = total_soft_plan_constraints > 0;
    any_hard_plan_constraints = num_hard_plan_constraints > 0;

}




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



double ScoringFunction::compute_plan_score(const Plan &plan, bool const is_final) const{
    double plan_score = 0.0;

    // add the score from each constraint 
    for(auto const &constraint_ptr: plan_constraint_ptrs){
        plan_score += constraint_ptr->compute_plan_constraint_score(plan);
    }
    // if final then return 
    if(is_final) return plan_score;

    // if not final then add those constraints (for now just pop_temper)
    for(auto const &constraint_ptr: non_plan_constraint_ptrs){
        plan_score += constraint_ptr->compute_plan_constraint_score(plan);
    }

    return plan_score;
}


double ScoringFunction::compute_merged_plan_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const{
    double plan_score = 0.0;

    // add the score from each constraint 
    for(auto const &constraint_ptr: plan_constraint_ptrs){
        plan_score += constraint_ptr->compute_merged_plan_constraint_score(plan, region1_id, region2_id);
    }
    // if final then return 
    if(is_final) return plan_score;

    // if not final then add those constraints (for now just pop_temper)
    for(auto const &constraint_ptr: non_plan_constraint_ptrs){
        plan_score += constraint_ptr->compute_merged_plan_constraint_score(plan, region1_id, region2_id);
    }

    return plan_score;
}

std::pair<bool, double> ScoringFunction::compute_hard_plan_constraints_score(const Plan &plan) const{
    if(!any_hard_plan_constraints) return std::make_pair(true, 0.0);
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


std::pair<bool, double> ScoringFunction::compute_hard_merged_plan_constraints_score(
    const Plan &plan, int const region1_id, int const region2_id
) const{
    if(!any_hard_plan_constraints) return std::make_pair(true, 0.0);

    double plan_score = 0.0;
 
    for(auto const &constraint_ptr: hard_plan_constraint_ptrs){
        auto result = constraint_ptr->compute_merged_plan_constraint_score(plan, region1_id, region2_id);
        // if false then failed to satisfy constraint so return false
        if(!result.first){
            return std::make_pair(false, 0.0);
        }else{
            plan_score += result.second;
        }
    }

    return std::make_pair(true, plan_score);
}