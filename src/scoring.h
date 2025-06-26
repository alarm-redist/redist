#pragma once
#ifndef SCORING_H
#define SCORING_H

#include "base_plan_type.h"
#include "redist_types.h"

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
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const;
};



class GroupHingeConstraint : public RegionConstraint {
    private:
        double const strength;
        arma::vec const tgts_group;
        arma::uvec const group_pop;
        arma::uvec const total_pop;

    public:
        GroupHingeConstraint(double const strength,  arma::vec const &tgts_group, arma::uvec const &group_pop,
            arma::uvec const &total_pop, bool const score_districts_only) :
            RegionConstraint(score_districts_only),
            strength(strength),
            tgts_group(tgts_group),
            group_pop(group_pop), 
            total_pop(total_pop) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const;
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
};

// scoring function 
class ScoringFunction {
    private:
        std::vector<std::unique_ptr<RegionConstraint>> region_constraint_ptrs; // These are constraints called on every split
        std::vector<std::unique_ptr<RegionConstraint>> non_final_region_constraint_ptrs; // these are constraints that are not called on the final round
        std::vector<std::unique_ptr<HardPlanConstraint>> hard_plan_constraint_ptrs; // These are hard constraints applied to an entire plan. Current not used in SMC code, only diagnostics  
    
    public:
        ScoringFunction(
            const MapParams &map_params, Rcpp::List const &constraints, 
            double const pop_temper);

        const MapParams &map_params;
        int num_non_final_soft_constraints; // applied to all but final split
        int num_final_soft_constraints; // applied to only final split
        int all_rounds_soft_constraints; // applied to all splits
        int total_soft_constraints; // total constraints 
        bool any_soft_constraints;
        int num_hard_constraints;

        double compute_region_score(const Plan &plan, int const region_id, bool const is_final) const;
        double compute_merged_region_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const;
        std::pair<bool, double> compute_hard_plan_constraints_score(const Plan &plan) const; // false if hard constraint failed 
};


#endif
