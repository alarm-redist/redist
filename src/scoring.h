#pragma once
#ifndef SCORING_H
#define SCORING_H

#include "base_plan_type.h"
#include "redist_types.h"
#include "map_calc.h"


class Plan;

/****************
 * Simply Hard Constraint Functions
 ********************/

// Checks if plan is connected and if not returns first disconnected region
template <typename PlanID>
std::pair<bool, int> is_plan_connected(
    Graph const &g, CircularQueue<int> &vertex_queue,
    std::vector<bool> &regions_visited, 
    const PlanID &region_ids
) {

    return std::make_pair(true, -1);
}

/****************
 * Constraint Functions
 ********************/



/****************
 * Constraint Classes
 ********************/

class RegionConstraint {
    public:
        RegionConstraint(bool const score_districts_only): score_districts_only(score_districts_only) {};  
        virtual ~RegionConstraint() = default;
        // whether or not to score districts only 
        bool score_districts_only;
        // log constraint for one region
        virtual double compute_region_constraint_score(const Plan &plan, int region_id) const = 0;
        // log constraint for region made by merging region 1 and 2
        virtual double compute_merged_region_constraint_score(const Plan &plan, int region1_id, int region2_id) const = 0;
};
    

// population tempering constraint 
class PopTemperConstraint : public RegionConstraint {
    private:
        double const target;
        int const ndists;
        double const pop_temper;

    public:
        PopTemperConstraint(
            double const target, int const ndists, const MapParams &map_params, 
            const double pop_temper, bool const score_districts_only = false) :
            RegionConstraint(score_districts_only),
            target(target),
            ndists(ndists),
            pop_temper(pop_temper) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};


class PopDevConstraint : public RegionConstraint {
    private:
        double const strength;
        double const parity;
        arma::uvec const total_pop;


    public:
        PopDevConstraint(
            double const strength, 
            double const parity, arma::uvec const &total_pop,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            parity(parity),
            total_pop(total_pop) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};




class StatusQuoConstraint : public RegionConstraint {
    private:
        double const strength;
        arma::uvec const current;
        arma::uvec const pop;
        int const ndists;
        int const n_current;
        int const V;

    public:
        StatusQuoConstraint(
            double const strength, 
            arma::uvec const &current, arma::uvec const &pop,
            int const ndists, int const n_current, int const V,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            current(current),
            pop(pop),
            ndists(ndists),
            n_current(n_current),
            V(V) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};



class SegregationConstraint : public RegionConstraint {
    private:
        double const strength;
        const arma::uvec grp_pop;
        const arma::uvec total_pop;
        int const V;

    public:
        SegregationConstraint(
            double const strength, 
            arma::uvec const &grp_pop, arma::uvec const &total_pop,
            int const V,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            grp_pop(grp_pop),
            total_pop(total_pop),
            V(V) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};

class GroupPowerConstraint : public RegionConstraint {
    private:
        double const strength;
        int const V;
        arma::uvec const grp_pop;
        arma::uvec const total_pop;
        double const tgt_grp;
        double const tgt_other;
        double const pow;


    public:
        GroupPowerConstraint(
            double const strength, 
            int const V,
            arma::uvec const &grp_pop, arma::uvec const &total_pop,
            double const tgt_grp, double const tgt_other, double const pow,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            V(V),
            grp_pop(grp_pop),
            total_pop(total_pop),
            tgt_grp(tgt_grp), tgt_other(tgt_other), pow(pow)
            {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};

class GroupHingeConstraint : public RegionConstraint {
    private:
        double const strength;
        int const V;
        arma::vec const tgts_group;
        arma::uvec const group_pop;
        arma::uvec const total_pop;

    public:
        GroupHingeConstraint(double const strength, int const V, arma::vec const &tgts_group, arma::uvec const &group_pop,
            arma::uvec const &total_pop, bool const score_districts_only) :
            RegionConstraint(score_districts_only),
            strength(strength),
            V(V), 
            tgts_group(tgts_group),
            group_pop(group_pop), 
            total_pop(total_pop) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};



class IncumbentConstraint : public RegionConstraint {
    private:
        double const strength;
        arma::uvec const incumbents;

    public:
        IncumbentConstraint(double const strength, const arma::uvec &incumbents, bool const score_districts_only) :
            RegionConstraint(score_districts_only),
            strength(strength),
            incumbents(incumbents) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};



class SplitsConstraint : public RegionConstraint {
    private:
        double const strength;
        arma::uvec const admin_units;
        int const n_admin_units;
        bool const smc;

    public:
        SplitsConstraint(
            double const strength, 
            arma::uvec const &admin_units, int const n_admin_units, bool const smc,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            admin_units(admin_units),
            n_admin_units(n_admin_units),
            smc(smc) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};


class MultisplitsConstraint : public RegionConstraint {
    private:
        double const strength;
        arma::uvec const admin_units;
        int const n_admin_units;
        bool const smc;

    public:
        MultisplitsConstraint(
            double const strength, 
            arma::uvec const &admin_units, int const n_admin_units, bool const smc,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            admin_units(admin_units),
            n_admin_units(n_admin_units),
            smc(smc) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};


class TotalSplitsConstraint : public RegionConstraint {
    private:
        double const strength;
        arma::uvec const admin_units;
        int const n_admin_units;
        bool const smc;

    public:
        TotalSplitsConstraint(
            double const strength, 
            arma::uvec const &admin_units, int const n_admin_units, bool const smc,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            admin_units(admin_units),
            n_admin_units(n_admin_units),
            smc(smc) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};


class PolsbyConstraint : public RegionConstraint {
    private:
        double const strength;
        int const V;
        arma::ivec const from;
        arma::ivec const to;
        arma::vec const area; 
        arma::vec const perimeter;
        

    public:
        PolsbyConstraint(
            double const strength, int const V,
            arma::ivec const &from, arma::ivec const &to, arma::vec const &area,
            arma::vec const &perimeter, 
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            V(V),
            from(from), to(to), area(area), perimeter(perimeter)
            {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};


class CustomRegionConstraint : public RegionConstraint {
    private:
        double const strength;
        Rcpp::Function const fn;        

    public:
        CustomRegionConstraint(
            double const strength, Rcpp::Function fn,
            bool const score_districts_only
        ) :
            RegionConstraint(score_districts_only),
            strength(strength),
            fn(Rcpp::clone(fn))
            {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;
};


// Soft plan constraint 
// These constraints should take an entire plan and return a finite score  
class SoftPlanConstraint {
    public:
        virtual ~SoftPlanConstraint() = default;
        // computes score for a plan 
        virtual double compute_plan_constraint_score(const Plan &plan) const = 0;
        virtual double compute_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const = 0;

};

// custom plan constraint. Assumes an Rcpp::Function is passed in
// which takes a vector 
class CustomPlanConstraint : public SoftPlanConstraint {
    private:
        double const strength;
        Rcpp::Function const fn; 

    public:
        CustomPlanConstraint(double const strength, Rcpp::Function fn):
            strength(strength),
            fn(Rcpp::clone(fn))
            {};
        // computes score for a plan 
        double compute_plan_constraint_score(const Plan &plan) const override;
        double compute_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;

};

// Hard plan constraint 
// These constraints should take an entire plan and potentially return a value of infinity 
// which means the plan has zero target probability  
class HardPlanConstraint {
    public:
        virtual ~HardPlanConstraint() = default;
        // computes score for a plan potentially allowing for infinity to be returned 
        // A value of false is equivalent to returning infinity so e^-score = 0
        virtual std::pair<bool, double> compute_plan_constraint_score(const Plan &plan) const = 0;
        virtual std::pair<bool, double> compute_merged_plan_constraint_score(
            const Plan &plan, int const region1_id, int const region2_id
        ) const = 0;

};


// Custom hard constraint
// Assumes that fn is a function which takes plan assignments and region sizes
// and returns a boolean 
class CustomHardPlanConstraint : public HardPlanConstraint {
    private:
        Rcpp::Function const fn;

    public:
        CustomHardPlanConstraint(Rcpp::Function fn):
        fn(Rcpp::clone(fn))
            {};

        std::pair<bool, double> compute_plan_constraint_score(const Plan &plan) const override;
        std::pair<bool, double> compute_merged_plan_constraint_score(
            const Plan &plan, int const region1_id, int const region2_id
        ) const override;
};

// Hard constraint which checks the plan has the correct number of valid sized districts
class ValidDistrictsConstraint : public HardPlanConstraint {
    private:
        MapParams const &map_params;

    public:
        ValidDistrictsConstraint(MapParams const &map_params):
        HardPlanConstraint(),
        map_params(map_params){};

        std::pair<bool, double> compute_plan_constraint_score(const Plan &plan) const override;
        std::pair<bool, double> compute_merged_plan_constraint_score(
            const Plan &plan, int const region1_id, int const region2_id
        ) const override {throw Rcpp::exception("ValidDistrictsConstraint Merged version Not implemented yet!\n");};
};

// scoring function 
class ScoringFunction {
    private:
        std::vector<std::unique_ptr<RegionConstraint>> region_constraint_ptrs; // These are constraints called on every split
        std::vector<std::unique_ptr<RegionConstraint>> non_final_region_constraint_ptrs; // these are constraints that are not called on the final round
        std::vector<std::unique_ptr<SoftPlanConstraint>> plan_constraint_ptrs;
        std::vector<std::unique_ptr<SoftPlanConstraint>> non_plan_constraint_ptrs;
        std::vector<std::unique_ptr<HardPlanConstraint>> hard_plan_constraint_ptrs; // These are hard constraints applied to an entire plan. Current not used in SMC code, only diagnostics  
    
    public:
        // rho and district_rho_only help determine computing compactness
        // we only compute log spanning tree if rho != 1 and if 
        // district_rho_only is true 
        // the smc is a legacy flag needed for splits. 
        // Ideally update functions and remove in the future 
        ScoringFunction(
            const MapParams &map_params, Rcpp::List const &constraints, 
            // double const rho, bool const district_rho_only,
            double const pop_temper, bool const smc);

        const MapParams &map_params;
        // 
        // double const excess_rho;
        // bool const any_excess_rho;
        // bool const district_rho_only;

        // counts region constraints
        int num_non_final_soft_region_constraints; // applied to all but final split
        int num_final_soft_region_constraints; // applied to only final split
        int all_rounds_soft_region_constraints; // applied to all splits
        int total_soft_region_constraints;
        // counts plan constraints 
        int num_non_final_soft_plan_constraints; // applied to all but final split
        int num_final_soft_plan_constraints; // applied to only final split
        int all_rounds_soft_plan_constraints; // applied to all splits
        int total_soft_plan_constraints;
        int num_hard_plan_constraints; //

        int total_soft_constraints; // total constraints 
        bool any_soft_region_constraints;
        bool any_soft_plan_constraints;
        bool any_hard_plan_constraints;

        // Rcpp::Function objects are not thread safe so special care is needed
        // For those. We make booleans to keep track of that 
        bool any_soft_custom_constraints;
        bool any_hard_custom_constraints;
        

        // scores individual regions
        double compute_region_score(const Plan &plan, int const region_id, bool const is_final) const;
        double compute_merged_region_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const;
        // scores plans 
        // soft score - ie always finite
        double compute_plan_score(const Plan &plan, bool const is_final) const;
        double compute_merged_plan_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const;
        // hard score - ie possible to be infinity which means probability zero
        std::pair<bool, double> compute_hard_plan_constraints_score(const Plan &plan) const; // false if hard constraint failed 
        std::pair<bool, double> compute_hard_merged_plan_constraints_score(
            const Plan &plan, int const region1_id, int const region2_id
        ) const; // false if hard constraint failed 

        // check if the two new regions or the plan trigger any hard constraints 
        bool new_split_ok(Plan const &plan, RegionID const region1_id, RegionID const region2_id);
};


#endif
