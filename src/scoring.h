#pragma once
#ifndef SCORING_H
#define SCORING_H

#include "base_plan_type.h"
#include "gredist_types.h"

class Constraint {
    public:
        virtual ~Constraint() = default;  
        // log constraint for one region
        virtual double compute_region_constraint_score(const Plan &plan, int region_id) const = 0;
        // log constraint for region made by merging region 1 and 2
        virtual double compute_merged_region_constraint_score(const Plan &plan, int region1_id, int region2_id) const = 0;
};
    

// population tempering constraint 
class PopTemperConstraint : public Constraint {
    private:
        double const target;
        int const ndists;
        double const pop_temper;
        MapParams const &map_params;

    public:
        PopTemperConstraint(double const target, int const ndists, const MapParams &map_params, const double pop_temper) :
            target(target),
            ndists(ndists),
            pop_temper(pop_temper), 
            map_params(map_params) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const;
};



class GroupHingeConstraint : public Constraint {
    private:
        double const strength;
        arma::vec const tgts_group;
        arma::uvec const group_pop;
        arma::uvec const total_pop;

    public:
        GroupHingeConstraint(double const strength,  arma::vec const &tgts_group, arma::uvec const &group_pop,
            arma::uvec const &total_pop) :
            strength(strength),
            tgts_group(tgts_group),
            group_pop(group_pop), 
            total_pop(total_pop) {}
    
        double compute_region_constraint_score(const Plan &plan, int const region_id) const;
        // log constraint for region made by merging region 1 and 2
        double compute_merged_region_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const;
};


// scoring function 
class ScoringFunction {
    private:
        std::vector<std::unique_ptr<Constraint>> constraint_ptrs; // These are constraints called on every split
        std::vector<std::unique_ptr<Constraint>> non_final_plan_constraint_ptrs; // these are constraints that are not called on the final round
    
    public:
        ScoringFunction(
            const MapParams &map_params, Rcpp::List const &constraints, 
            double const pop_temper,
            bool const score_districts_only = false);

        int num_non_final_constraints; // applied to all but final split
        int num_final_constraints; // applied to only final split
        int all_rounds_constraints; // applied to all splits
        bool any_constraints;
        bool score_districts_only; // If the score should only apply to districts (for backwards compatibility)

        double compute_region_score(const Plan &plan, int const region_id, bool const is_final) const;
        double compute_merged_region_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const;
};


#endif