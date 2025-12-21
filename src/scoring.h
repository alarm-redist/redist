#pragma once
#ifndef SCORING_H
#define SCORING_H

#include "base_plan_type.h"
#include "redist_types.h"
#include "map_calc.h"


class Plan;

/****************
 * Simple Hard Constraint Functions
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
    private:
        double const strength; // constraint strength 

    public:
        RegionConstraint(
            bool const score_districts_only, double const strength, 
            bool const hard_constraint, double const hard_threshold): 
            score_districts_only(score_districts_only),
            strength(strength),
            hard_constraint(hard_constraint),
            hard_threshold(hard_threshold) {}; 

        RegionConstraint(
            Rcpp::List const &a_constraint
        ): 
            score_districts_only(a_constraint.containsElementNamed("only_districts") ? as<bool>(a_constraint["only_districts"]) : false),
            strength(a_constraint.containsElementNamed("strength") ? as<double>(a_constraint["strength"]) : 1.0),
            hard_constraint(a_constraint.containsElementNamed("hard_constraint") ? as<bool>(a_constraint["hard_constraint"]) : false),
            hard_threshold(a_constraint.containsElementNamed("hard_threshold") ? as<double>(a_constraint["hard_threshold"]) : 0.0) 
            {};

        virtual ~RegionConstraint() = default;
        
        bool const score_districts_only; // whether or not to score districts only 
        bool const hard_constraint; // whether or not this is a hard constraint
        double const hard_threshold; // If hard constraint then the threshold for becoming zero

        // raw log constraint for one region
        virtual double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const = 0;
        // raw log constraint for region made by merging region 1 and 2
        virtual double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const = 0;

        // Returns whether region is valid and score
        // - value of true means satisfies hard constraint, false means probability zero
        // - double is the constraint score 
        std::pair<bool, double> compute_region_score(const Plan &plan, int region_id) const;
        std::pair<bool, double> compute_merged_region_score(const Plan &plan, int region1_id, int region2_id) const;

        // returns 0 if hard constraint otherwise returns soft score
        double compute_soft_region_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region_id) const;

        // just true or false for if the region is valid 
        bool region_constraint_ok(const Plan &plan, int region_id) const;

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
            RegionConstraint(score_districts_only, 1.0, false, 0.0),
            target(target),
            ndists(ndists),
            pop_temper(pop_temper) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};


class PopDevConstraint : public RegionConstraint {
    private:
        double const parity;
        arma::uvec const total_pop;


    public:
        PopDevConstraint(
            Rcpp::List const &a_constraint,
            double const parity, arma::uvec const &total_pop
        ) :
            RegionConstraint(a_constraint),
            parity(parity),
            total_pop(total_pop) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};




class StatusQuoConstraint : public RegionConstraint {
    private:
        arma::uvec const current;
        arma::uvec const pop;
        int const ndists;
        int const n_current;
        int const V;

    public:
        StatusQuoConstraint(
            Rcpp::List const &a_constraint, 
            arma::uvec const &current, arma::uvec const &pop,
            int const ndists, int const n_current, int const V
        ) :
            RegionConstraint(a_constraint),
            current(current),
            pop(pop),
            ndists(ndists),
            n_current(n_current),
            V(V) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};



class SegregationConstraint : public RegionConstraint {
    private:
        const arma::uvec grp_pop;
        const arma::uvec total_pop;
        int const V;

    public:
        SegregationConstraint(
            Rcpp::List const &a_constraint,
            arma::uvec const &grp_pop, arma::uvec const &total_pop,
            int const V
        ) :
            RegionConstraint(a_constraint),
            grp_pop(grp_pop),
            total_pop(total_pop),
            V(V) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};

class GroupPowerConstraint : public RegionConstraint {
    private:
        int const V;
        arma::uvec const grp_pop;
        arma::uvec const total_pop;
        double const tgt_grp;
        double const tgt_other;
        double const pow;


    public:
        GroupPowerConstraint(
            Rcpp::List const &a_constraint, 
            int const V,
            arma::uvec const &grp_pop, arma::uvec const &total_pop,
            double const tgt_grp, double const tgt_other, double const pow
        ) :
            RegionConstraint(a_constraint),
            V(V),
            grp_pop(grp_pop),
            total_pop(total_pop),
            tgt_grp(tgt_grp), tgt_other(tgt_other), pow(pow)
            {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};

class GroupHingeConstraint : public RegionConstraint {
    private:
        int const V;
        arma::vec const tgts_group;
        arma::uvec const group_pop;
        arma::uvec const total_pop;

    public:
        GroupHingeConstraint(
            Rcpp::List const &a_constraint, 
            int const V, arma::vec const &tgts_group, arma::uvec const &group_pop,
            arma::uvec const &total_pop
        ):
            RegionConstraint(a_constraint),
            V(V), 
            tgts_group(tgts_group),
            group_pop(group_pop), 
            total_pop(total_pop) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};



class IncumbentConstraint : public RegionConstraint {
    private:
        arma::uvec const incumbents; // NOTE: incumbents is 1-indexed

    public:
        IncumbentConstraint(
            Rcpp::List const &a_constraint, 
            const arma::uvec &incumbents
            ) :
            RegionConstraint(a_constraint),
            incumbents(incumbents) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};



class SplitsConstraint : public RegionConstraint {
    private:
        arma::uvec const admin_units;
        int const n_admin_units;
        bool const smc;

    public:
        SplitsConstraint(
            Rcpp::List const &a_constraint, 
            arma::uvec const &admin_units, int const n_admin_units, bool const smc
        ) :
            RegionConstraint(a_constraint),
            admin_units(admin_units),
            n_admin_units(n_admin_units),
            smc(smc) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};


class MultisplitsConstraint : public RegionConstraint {
    private:
        arma::uvec const admin_units;
        int const n_admin_units;
        bool const smc;

    public:
        MultisplitsConstraint(
            Rcpp::List const &a_constraint, 
            arma::uvec const &admin_units, int const n_admin_units, bool const smc
        ) :
            RegionConstraint(a_constraint),
            admin_units(admin_units),
            n_admin_units(n_admin_units),
            smc(smc) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};


class TotalSplitsConstraint : public RegionConstraint {
    private:
        arma::uvec const admin_units;
        int const n_admin_units;
        bool const smc;

    public:
        TotalSplitsConstraint(
            Rcpp::List const &a_constraint, 
            arma::uvec const &admin_units, int const n_admin_units, bool const smc
        ) :
            RegionConstraint(a_constraint),
            admin_units(admin_units),
            n_admin_units(n_admin_units),
            smc(smc) {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};


class PolsbyConstraint : public RegionConstraint {
    private:
        int const V;
        arma::ivec const from;
        arma::ivec const to;
        arma::vec const area; 
        arma::vec const perimeter;
        

    public:
        PolsbyConstraint(
            Rcpp::List const &a_constraint, 
            int const V,
            arma::ivec const &from, arma::ivec const &to, arma::vec const &area,
            arma::vec const &perimeter
        ) :
            RegionConstraint(a_constraint),
            V(V),
            from(from), to(to), area(area), perimeter(perimeter)
            {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};


class CustomRegionConstraint : public RegionConstraint {
    private:
        Rcpp::Function const fn;
        mutable Rcpp::IntegerVector rcpp_plan_wrap;    

    public:
        CustomRegionConstraint(
            Rcpp::List const &a_constraint, 
            int const V,
            Rcpp::Function fn
        ) :
            RegionConstraint(a_constraint),
            fn(Rcpp::clone(fn)),
            rcpp_plan_wrap(V)
            {}
    
        double compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const override;
        // log constraint for region made by merging region 1 and 2
        double compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const override;
};


// Soft plan constraint 
// These constraints should take an entire plan and return a finite score  
class PlanConstraint {
    private:
        double const strength; // constraint strength 

    public:
        PlanConstraint(
            double const strength, std::vector<bool> const &num_regions_to_score,
            bool const hard_constraint, double const hard_threshold):
            strength(strength),
            num_regions_to_score(num_regions_to_score),
            hard_constraint(hard_constraint),
            hard_threshold(hard_threshold) {};


        PlanConstraint(
            Rcpp::List const &a_constraint, int const ndists
            ):
            strength(a_constraint.containsElementNamed("strength") ? as<double>(a_constraint["strength"]) : 1.0),
            num_regions_to_score([ndists, &a_constraint]() {
                    // vector where index i is true iff i seats is a district 
                    std::vector<bool> num_regions_to_score(ndists + 1, false);
                    if (a_constraint.containsElementNamed("nregions_to_score")){
                        // The vector in R is one indexed but c++ is 0 indexed so need to pad an 
                        // extra element
                        num_regions_to_score = Rcpp::as<std::vector<bool>>(a_constraint["nregions_to_score"]);
                        num_regions_to_score.insert(num_regions_to_score.begin(), false);
                    }
                    return num_regions_to_score;
                }()),
            hard_constraint(a_constraint.containsElementNamed("hard_constraint") ? as<bool>(a_constraint["hard_constraint"]) : false),
            hard_threshold(a_constraint.containsElementNamed("hard_threshold") ? as<double>(a_constraint["hard_threshold"]) : 0.0) 
          {};

        virtual ~PlanConstraint() = default;
        // attributes
        
        std::vector<bool> const num_regions_to_score; // Whether or not to score plans with that many regions
        bool const hard_constraint; // whether or not this is a hard constraint
        double const hard_threshold; // If hard constraint then the threshold for becoming zero
        

        // print. Just for debugging
        virtual void print() const;
        // computes score for a plan 
        virtual double compute_raw_plan_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const = 0;
        virtual double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const = 0;
        // 
        
        std::pair<bool, double> compute_plan_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const;
        std::pair<bool, double> compute_merged_plan_score(const Plan &plan, int const region1_id, int const region2_id) const;

        // just checks if plan is ok, doesn't save score
        bool plan_constraint_ok(const Plan &plan) const;

};


// Splits but on the entire plan
// assume admin is 1 indexed and only has values 1:num_admin_units
class PlanSplitsConstraint : public PlanConstraint {
    private:
        arma::uvec const admin_units;
        int const num_admin_units;
        std::vector<std::vector<int>> const admin_vertex_lists;
        mutable std::vector<int> region_reindex_vec;

    public:
        PlanSplitsConstraint(
            double const strength, int const ndists,
            arma::uvec const admin_units, 
            std::vector<std::vector<int>> const &admin_vertex_lists,
            std::vector<bool> const &num_regions_to_score,
            bool const hard_constraint, double const hard_threshold):
            PlanConstraint(strength, num_regions_to_score, hard_constraint, hard_threshold),
            admin_units(admin_units), 
            num_admin_units(arma::max(admin_units)),
            admin_vertex_lists(admin_vertex_lists),
            region_reindex_vec(ndists)
            {};
        // computes score for a plan 
        double compute_raw_plan_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const override;
        double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;

};


// Total Splits but of the entire plan
// assume admin is 1 indexed and only has values 1:num_admin_units
class TotalPlanSplitsConstraint : public PlanConstraint {
    private:
        arma::uvec const admin_units;
        int const num_admin_units;
        std::vector<std::vector<int>> const admin_vertex_lists;
        mutable std::vector<int> region_reindex_vec;
        mutable std::vector<std::set<int>> admin_unit_regions;

    public:
        TotalPlanSplitsConstraint(
            double const strength, int const ndists,
            arma::uvec const admin_units, 
            std::vector<std::vector<int>> const &admin_vertex_lists,
            std::vector<bool> const &num_regions_to_score,
            bool const hard_constraint, double const hard_threshold):
            PlanConstraint(strength, num_regions_to_score, hard_constraint, hard_threshold),
            admin_units(admin_units), 
            num_admin_units(arma::max(admin_units)),
            admin_vertex_lists(admin_vertex_lists),
            region_reindex_vec(ndists),
            admin_unit_regions(std::vector<std::set<int>>(num_admin_units))
            {};
        // computes score for a plan 
        double compute_raw_plan_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const override;
        double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;

};


// Counts how many districts in a plan have more than 1 incumbent in it 
class PlanIncumbentConstraint : public PlanConstraint {
    private:
        std::vector<bool> const is_district; // of length total_seats that says whether or not that size is a district
        arma::uvec const incumbents;
        mutable std::vector<int> region_reindex_vec;
        mutable std::vector<int> region_incumbent_counts;
        mutable std::vector<bool> region_is_district; // for specific plan checks if region is district 

    public:
        PlanIncumbentConstraint(
            double const strength, int const ndists,
            std::vector<bool> const &is_district,
            arma::uvec const &incumbents,
            std::vector<bool> const &num_regions_to_score,
            bool const hard_constraint, double const hard_threshold):
            PlanConstraint(strength, num_regions_to_score, hard_constraint, hard_threshold),
            is_district(is_district), 
            incumbents(incumbents),
            region_reindex_vec(ndists),
            region_incumbent_counts(ndists),
            region_is_district(ndists)
            {};
        // computes score for a plan 
        double compute_raw_plan_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const override;
        double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;

};


// Counts the number of districts where group_pop/target_pop >= min_frac
class MinGroupFracConstraint : public PlanConstraint {
    private:
        std::vector<bool> const is_district;
        mutable std::vector<std::vector<double>> plan_group_pops;
        mutable std::vector<std::vector<double>> plan_total_pops;
        mutable std::vector<bool> region_ids_to_count;
        mutable std::vector<int> region_reindex_vec;
        std::vector<arma::vec> const group_pops;
        std::vector<arma::vec> const total_pops;
        std::vector<double> const min_fracs;
        int const num_populations;


    public:
        MinGroupFracConstraint(
            double const strength, int const ndists,
            std::vector<bool> const &is_district,
            std::vector<arma::vec> const &group_pops, 
            std::vector<arma::vec> const &total_pops, 
            std::vector<double> const &min_fracs,
            int const num_populations, 
            std::vector<bool> const &num_regions_to_score,
            bool const hard_constraint, double const hard_threshold):
            PlanConstraint(strength, num_regions_to_score, hard_constraint, hard_threshold),
            is_district(is_district),
            plan_group_pops(num_populations, std::vector<double>(ndists, 0.0)), 
            plan_total_pops(num_populations, std::vector<double>(ndists, 0.0)), 
            region_ids_to_count(ndists, false),
            region_reindex_vec(ndists),
            group_pops(group_pops), total_pops(total_pops), min_fracs(min_fracs),
            num_populations(num_populations)
            {};
        // computes score for a plan 
        double compute_raw_plan_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const override;
        double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;

};

// custom plan constraint. Assumes an Rcpp::Function is passed in
// which takes a vector 
class CustomPlanConstraint : public PlanConstraint {
    private:
        Rcpp::Function const fn; 

    public:
        CustomPlanConstraint(double const strength, Rcpp::Function fn, 
            std::vector<bool> const &num_regions_to_score,
            bool const hard_constraint, double const hard_threshold):
            PlanConstraint(strength, num_regions_to_score, hard_constraint, hard_threshold),
            fn(Rcpp::clone(fn))
            {};
        // computes score for a plan 
        double compute_raw_plan_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const override;
        double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const override;

};


// Hard constraint which checks the plan has the correct number of valid sized districts
class ValidDistrictsConstraint : public PlanConstraint {
    private:
        MapParams const &map_params;

    public:
        ValidDistrictsConstraint(MapParams const &map_params):
        PlanConstraint(1, 
    [&map_params]() {
                // we only score full plans 
                std::vector<bool> num_regions_to_score(map_params.ndists + 1, false);
                num_regions_to_score[map_params.ndists] = true;
                return num_regions_to_score;
            }(),
            true, .5),
        map_params(map_params){};

        double compute_raw_plan_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
        ) const override;
        double compute_raw_merged_plan_constraint_score(
            const Plan &plan, int const region1_id, int const region2_id
        ) const override {throw Rcpp::exception("ValidDistrictsConstraint Merged version Not implemented yet!\n");};
};

// scoring function 
class ScoringFunction {
    private:
        std::vector<std::unique_ptr<RegionConstraint>> region_constraint_ptrs; // These are constraints called on every split
        std::vector<std::unique_ptr<RegionConstraint>> non_final_region_constraint_ptrs; // these are constraints that are not called on the final round
        std::vector<std::unique_ptr<PlanConstraint>> plan_constraint_ptrs; // Constraints called on every split



    
    public:
        // rho and district_rho_only help determine computing compactness
        // we only compute log spanning tree if rho != 1 and if 
        // district_rho_only is true 
        // the smc is a legacy flag needed for splits. 
        // Ideally update functions and remove in the future 
        // thread id is for dealing with custom constraints 
        ScoringFunction(
            const MapParams &map_params, Rcpp::List const &constraints, 
            // double const rho, 
            double const pop_temper, bool const smc, int const thread_id = 0);

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
        int num_hard_region_constraints;
        // counts plan constraints 
        int total_soft_plan_constraints;
        int num_hard_plan_constraints; //

        int total_soft_constraints; // total constraints 
        bool any_soft_region_constraints;
        bool any_soft_plan_constraints;
        bool any_hard_plan_constraints;
        bool any_hard_region_constraints;
        bool any_hard_constraints;

        // Rcpp::Function objects are not thread safe so special care is needed
        // For those. We make booleans to keep track of that 
        bool any_soft_custom_constraints;
        bool any_hard_custom_constraints;

        // used for constraint splitter
        double compute_full_split_plan_soft_score(int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const split_region1, int const split_region2
        ) const;
        

        // scores individual regions
        std::pair<bool, double> compute_region_full_score(const Plan &plan, int const region_id, bool const is_final) const;
        double compute_region_soft_score(const Plan &plan, int const region_id, bool const is_final) const;
        std::pair<bool, double> compute_merged_region_full_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const;


        // scores plans 
        // false means probability zero 
        // soft score - ie always finite
        std::pair<bool, double> compute_plan_score(const Plan &plan) const;
        std::pair<bool, double> compute_merged_plan_score(const Plan &plan, int const region1_id, int const region2_id, bool const is_final) const;

        // check if the merged region triggers any hard constraints
        bool merged_region_ok(Plan const &plan, int const region1_id, int const region2_id, bool const is_final_split) const;
        // check any hard entire plan constraints on the merged plan
        bool entire_merged_plan_constraint_only_ok(Plan const &plan, int const region1_id, int const region2_id, bool const is_final_split) const;
        // check if the entire merged plan is ok
        bool merged_plan_ok(Plan const &plan, int const region1_id, int const region2_id, bool const is_final_split) const;
        // check if the two new regions or the plan trigger any hard constraints 
        bool new_split_ok(Plan const &plan, int const region1_id, int const region2_id, bool const is_final_split) const;
};


#endif
