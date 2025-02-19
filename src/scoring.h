#pragma once
#ifndef SCORING_H
#define SCORING_H

#include "base_plan_type.h"
#include "gredist_types.h"

class Constraint {
    public:
        virtual ~Constraint() = default;  
        // log constraint for one region
        virtual double compute_log_region_constraint(const Plan &plan, int region_id) const = 0;
        // log constraint for region made by merging region 1 and 2
        virtual double compute_log_merged_region_constraint(const Plan &plan, int region1_id, int region2_id) const = 0;
};
    


class PopTemperConstraint : public Constraint {
    private:
        double const pop_temper;
        double const target;
        int const ndists;
        MapParams const &map_params;

    public:
        PopTemperConstraint(double const target, int const ndists, const MapParams &map_params, const double pop_temper) :
            target(target),
            ndists(ndists),
            pop_temper(pop_temper), 
            map_params(map_params) {}
    
        double compute_log_region_constraint(const Plan &plan, int const region_id) const;
        // log constraint for region made by merging region 1 and 2
        double compute_log_merged_region_constraint(const Plan &plan, int const region1_id, int const region2_id) const;
};


// scoring function 
class ScoringFunction {
    private:
        std::vector<std::unique_ptr<Constraint>> constraint_ptrs;
        std::vector<std::unique_ptr<Constraint>> non_final_plan_constraint_ptrs; // these are constraints that are not called on the final round
    
    public:
        ScoringFunction(
            const MapParams &map_params, Rcpp::List const &constraints, 
            bool const do_pop_temper, double const pop_temper);

        double compute_log_region_score(const Plan &plan, int const region_id, bool const is_final) const;
        double compute_log_merged_region_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const;
};


#endif