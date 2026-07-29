#pragma once
#ifndef SCORING_H
#define SCORING_H


#include "map_calc.h"
#include "advanced_types.h"


class Plan;

inline std::vector<std::vector<int>> build_admin_vertex_lists(const Graph &g,
                                                       const arma::uvec &admin_units) {
    // we assume admin units is 1 indexed and if `k` units then values are in `1:k`
    int const num_counties = arma::max(admin_units);
    std::vector<std::vector<int>> admin_vertex_lists(num_counties);
    // nothing if only 1 county
    if (num_counties == 1)
        return admin_vertex_lists;
    // else walk through the graph and add each vertex to list for each unit
    int const V = g.size();

    for (int v = 0; v < V; v++) {
        int v_admin_unit = admin_units[v] - 1;
        admin_vertex_lists[v_admin_unit].push_back(v);
    }
    return admin_vertex_lists;
}

/****************
 * Simple Hard Constraint Functions
 ********************/

// Checks if plan is connected and if not returns first disconnected region
template <typename PlanID>
std::pair<bool, int> is_plan_connected(Graph const &g, CircularQueue<int> &vertex_queue,
                                       std::vector<bool> &regions_visited,
                                       const PlanID &region_ids) {

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
    RegionConstraint(Rcpp::List const &a_constraint, 
        int const ndists, int const total_seats)
        : strength(a_constraint.containsElementNamed("strength")
                       ? Rcpp::as<double>(a_constraint["strength"])
                       : 1.0),
    num_regions_to_score([ndists, &a_constraint]() {
                  // vector where index i == true means score plan with i regions
                  std::vector<bool> num_regions_to_score(ndists + 1, true);
                  if (a_constraint.containsElementNamed("nregions_to_score")) {
                      // The vector in R is one indexed but c++ is 0 indexed so need to pad an
                      // extra element
                      num_regions_to_score =
                          Rcpp::as<std::vector<bool>>(a_constraint["nregions_to_score"]);
                      num_regions_to_score.insert(num_regions_to_score.begin(), false);
                  }
                  return num_regions_to_score;
              }()),
    region_sizes_to_score([total_seats, &a_constraint]() {
                  // vector where index i == true means score regions with size i
                  std::vector<bool> region_sizes_to_score(total_seats + 1, true);
                  if (a_constraint.containsElementNamed("nseats_to_score")) {
                      // The vector in R is one indexed but c++ is 0 indexed so need to pad an
                      // extra element
                      region_sizes_to_score =
                          Rcpp::as<std::vector<bool>>(a_constraint["nseats_to_score"]);
                      region_sizes_to_score.insert(region_sizes_to_score.begin(), false);
                  }
                  return region_sizes_to_score;
              }()),
          hard_constraint(a_constraint.containsElementNamed("hard_constraint")
                              ? Rcpp::as<bool>(a_constraint["hard_constraint"])
                              : false),
          hard_threshold(a_constraint.containsElementNamed("hard_threshold")
                             ? Rcpp::as<double>(a_constraint["hard_threshold"])
                             : 0.0) {};

    virtual ~RegionConstraint() = default;

    std::vector<bool> const
        num_regions_to_score;    // Whether or not to score plans with that many regions
    std::vector<bool> const
        region_sizes_to_score;    // Whether or not to score regions with that many seats
    bool const hard_constraint;      // whether or not this is a hard constraint
    double const hard_threshold;     // If hard constraint then the threshold for becoming zero

    // raw log constraint for one region
    virtual double compute_raw_region_constraint_score(int const num_regions,
                                                       PlanVector const &region_ids,
                                                       RegionSizes const &region_sizes,
                                                       IntPlanAttribute const &region_pops,
                                                       int region_id) const = 0;
    // raw log constraint for region made by merging region 1 and 2
    virtual double compute_raw_merged_region_constraint_score(
        int const num_regions, PlanVector const &region_ids, RegionSizes const &region_sizes,
        IntPlanAttribute const &region_pops, int const region1_id,
        int const region2_id) const = 0;

    // Returns whether region is valid and score
    // - value of true means satisfies hard constraint, false means probability zero
    // - double is the constraint score
    std::pair<bool, double> compute_region_score(int const num_regions,
                                                       PlanVector const &region_ids,
                                                       RegionSizes const &region_sizes,
                                                       IntPlanAttribute const &region_pops,
                                                       int region_id) const;
    std::pair<bool, double> compute_region_score(const Plan &plan, int region_id) const;
    std::pair<bool, double> compute_merged_region_score(const Plan &plan, int region1_id,
                                                        int region2_id) const;

    // returns 0 if hard constraint otherwise returns soft score
    double compute_soft_region_score(int const num_regions, PlanVector const &region_ids,
                                     RegionSizes const &region_sizes,
                                     IntPlanAttribute const &region_pops,
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
    PopTemperConstraint(double const target, int const ndists, const MapParams &map_params,
                        const double pop_temper)
    : RegionConstraint(
          Rcpp::List::create(
              Rcpp::_["strength"] = 1.0,
              Rcpp::_["hard_constraint"] = false,
              Rcpp::_["hard_threshold"] = 0.0,
              Rcpp::_["nregions_to_score"] = [ndists]() {
                  // make every region but the last one true
                  std::vector<bool> out(ndists, true);
                  out[ndists - 1] = false;

                  return out;
              }()
          ),
          ndists,
          map_params.total_seats
      ),
      target(target),
      ndists(ndists),
      pop_temper(pop_temper) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class PopDevConstraint : public RegionConstraint {
  private:
    double const parity;
    std::vector<unsigned int> const total_pop;

  public:
    PopDevConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats), 
        parity(map_params.target), total_pop(map_params.pop) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class StatusQuoConstraint : public RegionConstraint {
  private:
    std::vector<unsigned int> const current;
    std::vector<unsigned int> const pop;
    int const ndists;
    int const n_current;
    int const V;

  public:
    StatusQuoConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats),
        current(Rcpp::as<std::vector<unsigned int>>(constr_inst["current"])), 
        pop(map_params.pop), ndists(map_params.ndists),
          n_current(Rcpp::as<int>(constr_inst["n_current"])), V(map_params.V) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class SegregationConstraint : public RegionConstraint {
  private:
    const arma::uvec grp_pop;
    const arma::uvec total_pop;
    int const V;

  public:
    SegregationConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats),
        grp_pop(Rcpp::as<arma::uvec>(constr_inst["group_pop"])), 
        total_pop(Rcpp::as<arma::uvec>(constr_inst["total_pop"])), 
        V(map_params.V) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class GroupPowerConstraint : public RegionConstraint {
  private:
    int const V;
    arma::uvec const grp_pop;
    arma::uvec const total_pop;
    double const tgt_grp;
    double const tgt_other;
    double const pow;

    GroupPowerConstraint(
        Rcpp::List const &constr_inst,
        MapParams const &map_params,
        arma::uvec grp_pop_,
        arma::uvec total_pop_,
        double const tgt_grp_,
        double const tgt_other_,
        double const pow_
    )
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats),
          V(map_params.V),
          grp_pop(std::move(grp_pop_)),
          total_pop(std::move(total_pop_)),
          tgt_grp(tgt_grp_),
          tgt_other(tgt_other_),
          pow(pow_) {}

  public:
    static std::unique_ptr<GroupPowerConstraint> from_group_power(
        Rcpp::List const &constr_inst,
        MapParams const &map_params
    ) {
        return std::unique_ptr<GroupPowerConstraint>(
            new GroupPowerConstraint(
                constr_inst,
                map_params,
                Rcpp::as<arma::uvec>(constr_inst["group_pop"]),
                Rcpp::as<arma::uvec>(constr_inst["total_pop"]),
                Rcpp::as<double>(constr_inst["tgt_group"]),
                Rcpp::as<double>(constr_inst["tgt_other"]),
                Rcpp::as<double>(constr_inst["pow"])
            )
        );
    }

    static std::unique_ptr<GroupPowerConstraint> from_competition(
        Rcpp::List const &constr_inst,
        MapParams const &map_params
    ) {
        arma::uvec dvote = Rcpp::as<arma::uvec>(constr_inst["dvote"]);
        arma::uvec rvote = Rcpp::as<arma::uvec>(constr_inst["rvote"]);
        arma::uvec total = dvote + rvote;

        return std::unique_ptr<GroupPowerConstraint>(
            new GroupPowerConstraint(
                constr_inst,
                map_params,
                std::move(dvote),
                std::move(total),
                0.5,
                0.5,
                Rcpp::as<double>(constr_inst["pow"])
            )
        );
    }

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class GroupHingeConstraint : public RegionConstraint {
  private:
    int const V;
    arma::vec const tgts_group;
    arma::uvec const group_pop;
    arma::uvec const total_pop;

  public:
    // This works for both group hinge and group inverse hinge
    GroupHingeConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats), 
        V(map_params.V), 
        tgts_group(Rcpp::as<arma::vec>(constr_inst["tgts_group"])), 
        group_pop(Rcpp::as<arma::uvec>(constr_inst["group_pop"])),
          total_pop(Rcpp::as<arma::uvec>(constr_inst["total_pop"])) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class IncumbentConstraint : public RegionConstraint {
  private:
    arma::uvec const incumbents; // NOTE: incumbents is 1-indexed

  public:
    IncumbentConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats),
        incumbents(Rcpp::as<arma::uvec>(constr_inst["incumbents"])) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class SplitsConstraint : public RegionConstraint {
  private:
    arma::uvec const admin_units;
    int const n_admin_units;
    bool const smc;

  public:
    SplitsConstraint(Rcpp::List const &constr_inst, MapParams const &map_params, bool const smc)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats),
        admin_units(Rcpp::as<arma::uvec>(constr_inst["admin"])),
          n_admin_units(Rcpp::as<int>(constr_inst["n"])), smc(smc) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class MultisplitsConstraint : public RegionConstraint {
  private:
    arma::uvec const admin_units;
    int const n_admin_units;
    bool const smc;

  public:
    MultisplitsConstraint(Rcpp::List const &constr_inst, MapParams const &map_params, bool const smc)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats),
        admin_units(Rcpp::as<arma::uvec>(constr_inst["admin"])),
          n_admin_units(Rcpp::as<int>(constr_inst["n"])), smc(smc) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class TotalSplitsConstraint : public RegionConstraint {
  private:
    arma::uvec const admin_units;
    int const n_admin_units;
    bool const smc;

  public:
    TotalSplitsConstraint(Rcpp::List const &constr_inst, MapParams const &map_params, bool const smc)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats),
        admin_units(Rcpp::as<arma::uvec>(constr_inst["admin"])),
          n_admin_units(Rcpp::as<int>(constr_inst["n"])), smc(smc) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class PolsbyConstraint : public RegionConstraint {
  private:
    int const V;
    arma::ivec const from;
    arma::ivec const to;
    arma::vec const area;
    arma::vec const perimeter;

  public:
    PolsbyConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : RegionConstraint(constr_inst, map_params.ndists, map_params.total_seats), 
        V(map_params.V), 
        from(Rcpp::as<arma::ivec>(constr_inst["from"])), 
        to(Rcpp::as<arma::ivec>(constr_inst["to"])), 
        area(Rcpp::as<arma::vec>(constr_inst["area"])),
          perimeter(Rcpp::as<arma::vec>(constr_inst["perimeter"])) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

class CustomRegionConstraint : public RegionConstraint {
  private:
    Rcpp::Function const fn;
    mutable Rcpp::IntegerVector rcpp_plan_wrap;

  public:
    CustomRegionConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : RegionConstraint(constr_inst,  map_params.ndists, map_params.total_seats),
         fn(Rcpp::clone(Rcpp::as<Rcpp::Function>(constr_inst["fn"]))), 
         rcpp_plan_wrap(map_params.V) {}

    double compute_raw_region_constraint_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops,
                                               int region_id) const override;
    // log constraint for region made by merging region 1 and 2
    double compute_raw_merged_region_constraint_score(int const num_regions,
                                                      PlanVector const &region_ids,
                                                      RegionSizes const &region_sizes,
                                                      IntPlanAttribute const &region_pops,
                                                      int const region1_id,
                                                      int const region2_id) const override;
};

// Soft plan constraint
// These constraints should take an entire plan and return a finite score
class PlanConstraint {
  private:
    double const strength; // constraint strength

  public:
    PlanConstraint(Rcpp::List const &constr_inst, int const ndists)
        : strength(constr_inst.containsElementNamed("strength")
                       ? Rcpp::as<double>(constr_inst["strength"])
                       : 1.0),
          num_regions_to_score([ndists, &constr_inst]() {
              // vector where index i is true iff i seats is a district
              std::vector<bool> num_regions_to_score(ndists + 1, true);
              if (constr_inst.containsElementNamed("nregions_to_score")) {
                  // The vector in R is one indexed but c++ is 0 indexed so need to pad an
                  // extra element
                  num_regions_to_score =
                      Rcpp::as<std::vector<bool>>(constr_inst["nregions_to_score"]);
                  num_regions_to_score.insert(num_regions_to_score.begin(), false);
              }
              return num_regions_to_score;
          }()),
          hard_constraint(constr_inst.containsElementNamed("hard_constraint")
                              ? Rcpp::as<bool>(constr_inst["hard_constraint"])
                              : false),
          hard_threshold(constr_inst.containsElementNamed("hard_threshold")
                             ? Rcpp::as<double>(constr_inst["hard_threshold"])
                             : 0.0) {};

    virtual ~PlanConstraint() = default;
    // attributes

    std::vector<bool> const
        num_regions_to_score;    // Whether or not to score plans with that many regions
    bool const hard_constraint;  // whether or not this is a hard constraint
    double const hard_threshold; // If hard constraint then the threshold for becoming zero

    // print. Just for debugging
    virtual void print() const;
    // computes score for a plan
    virtual double
    compute_raw_plan_constraint_score(int const num_regions, PlanVector const &region_ids,
                                      RegionSizes const &region_sizes,
                                      IntPlanAttribute const &region_pops) const = 0;
    virtual double compute_raw_merged_plan_constraint_score(const Plan &plan,
                                                            int const region1_id,
                                                            int const region2_id) const = 0;
    //

    std::pair<bool, double> compute_plan_score(int const num_regions,
                                               PlanVector const &region_ids,
                                               RegionSizes const &region_sizes,
                                               IntPlanAttribute const &region_pops) const;
    std::pair<bool, double> compute_merged_plan_score(const Plan &plan, int const region1_id,
                                                      int const region2_id) const;

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
    PlanSplitsConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : PlanConstraint(constr_inst, map_params.ndists),
          admin_units(Rcpp::as<arma::uvec>(constr_inst["admin"])), 
          num_admin_units(*std::max_element(admin_units.begin(), admin_units.end())),
          admin_vertex_lists(build_admin_vertex_lists(map_params.g, admin_units)), 
          region_reindex_vec(map_params.ndists) {};
    // computes score for a plan
    double
    compute_raw_plan_constraint_score(int const num_regions, PlanVector const &region_ids,
                                      RegionSizes const &region_sizes,
                                      IntPlanAttribute const &region_pops) const override;
    double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id,
                                                    int const region2_id) const override;
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
    TotalPlanSplitsConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : PlanConstraint(constr_inst,map_params.ndists),
          admin_units(Rcpp::as<arma::uvec>(constr_inst["admin"])), 
          num_admin_units(*std::max_element(admin_units.begin(), admin_units.end())),
          admin_vertex_lists(build_admin_vertex_lists(map_params.g, admin_units)),
          region_reindex_vec(map_params.ndists),
          admin_unit_regions(std::vector<std::set<int>>(num_admin_units)) {};
    // computes score for a plan
    double
    compute_raw_plan_constraint_score(int const num_regions, PlanVector const &region_ids,
                                      RegionSizes const &region_sizes,
                                      IntPlanAttribute const &region_pops) const override;
    double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id,
                                                    int const region2_id) const override;
};

// Counts how many districts in a plan have more than 1 incumbent in it
class PlanIncumbentConstraint : public PlanConstraint {
  private:
    std::vector<bool> const
        is_district; // of length total_seats that says whether or not that size is a district
    arma::uvec const incumbents;
    mutable std::vector<int> region_reindex_vec;
    mutable std::vector<int> region_incumbent_counts;
    mutable std::vector<bool>
        region_is_district; // for specific plan checks if region is district

  public:
    PlanIncumbentConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : PlanConstraint(constr_inst, map_params.ndists),
          is_district(map_params.is_district), 
          incumbents(Rcpp::as<arma::uvec>(constr_inst["incumbents"])), 
          region_reindex_vec(map_params.ndists),
          region_incumbent_counts(map_params.ndists), 
          region_is_district(map_params.ndists) {};
    // computes score for a plan
    double
    compute_raw_plan_constraint_score(int const num_regions, PlanVector const &region_ids,
                                      RegionSizes const &region_sizes,
                                      IntPlanAttribute const &region_pops) const override;
    double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id,
                                                    int const region2_id) const override;
};

// Counts the number of districts where group_pop/target_pop >= min_frac
class MinGroupFracConstraint : public PlanConstraint {
  private:
    std::vector<bool> const is_district;
    int const num_populations;
    mutable std::vector<std::vector<double>> plan_group_pops;
    mutable std::vector<std::vector<double>> plan_total_pops;
    mutable std::vector<bool> region_ids_to_count;
    mutable std::vector<int> region_reindex_vec;
    std::vector<arma::vec> const group_pops;
    std::vector<arma::vec> const total_pops;
    std::vector<double> const min_fracs;
    

  public:
    MinGroupFracConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : PlanConstraint(constr_inst, map_params.ndists),
          is_district(map_params.is_district),
          num_populations(Rcpp::as<double>(constr_inst["num_populations"])),
          plan_group_pops(num_populations, std::vector<double>(map_params.ndists, 0.0)),
          plan_total_pops(num_populations, std::vector<double>(map_params.ndists, 0.0)),
          region_ids_to_count(map_params.ndists, false), 
          region_reindex_vec(map_params.ndists),
          group_pops(Rcpp::as<std::vector<arma::vec>>(constr_inst["group_pops"])), 
          total_pops(Rcpp::as<std::vector<arma::vec>>(constr_inst["total_pops"])), 
          min_fracs(Rcpp::as<std::vector<double>>(constr_inst["min_fracs"])) {};
    // computes score for a plan
    double
    compute_raw_plan_constraint_score(int const num_regions, PlanVector const &region_ids,
                                      RegionSizes const &region_sizes,
                                      IntPlanAttribute const &region_pops) const override;
    double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id,
                                                    int const region2_id) const override;
};

// custom plan constraint. Assumes an Rcpp::Function is passed in
// which takes a vector
class CustomPlanConstraint : public PlanConstraint {
  private:
    Rcpp::Function const fn;

  public:
    CustomPlanConstraint(Rcpp::List const &constr_inst, MapParams const &map_params)
        : PlanConstraint(constr_inst, map_params.ndists),
          fn(Rcpp::clone(Rcpp::as<Rcpp::Function>(constr_inst["fn"]))) {};
    // computes score for a plan
    double
    compute_raw_plan_constraint_score(int const num_regions, PlanVector const &region_ids,
                                      RegionSizes const &region_sizes,
                                      IntPlanAttribute const &region_pops) const override;
    double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id,
                                                    int const region2_id) const override;
};

// Hard constraint which checks the plan has the correct number of valid sized districts
class ValidDistrictsConstraint : public PlanConstraint {
  private:
    MapParams const &map_params;

  public:
    ValidDistrictsConstraint(MapParams const &map_params)
        : PlanConstraint(
                    Rcpp::List::create(
              Rcpp::_["strength"] = 1.0,
              Rcpp::_["hard_constraint"] = true,
              Rcpp::_["hard_threshold"] = 0.5,
              Rcpp::_["nregions_to_score"] = [&map_params]() {
                  // we only score full plans
                  std::vector<bool> num_regions_to_score(map_params.ndists + 1, false);
                  num_regions_to_score[map_params.ndists] = true;
                  return num_regions_to_score;
              }()), 
              map_params.ndists
          ),
          map_params(map_params) {};

    double
    compute_raw_plan_constraint_score(int const num_regions, PlanVector const &region_ids,
                                      RegionSizes const &region_sizes,
                                      IntPlanAttribute const &region_pops) const override;
    double compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id,
                                                    int const region2_id) const override {
        throw Rcpp::exception("ValidDistrictsConstraint Merged version Not implemented yet!\n");
    };
};

// scoring function
class ScoringFunction {
  private:
    std::vector<std::unique_ptr<RegionConstraint>>
        region_constraint_ptrs; // 
    std::vector<std::unique_ptr<PlanConstraint>>
        plan_constraint_ptrs; // 

  public:
    // rho and district_rho_only help determine computing compactness
    // we only compute log spanning tree if rho != 1 and if
    // district_rho_only is true
    // the smc is a legacy flag needed for splits.
    // Ideally update functions and remove in the future
    // thread id is for dealing with custom constraints
    ScoringFunction(const MapParams &map_params, Rcpp::List const &constraints,
                    // double const rho,
                    double const pop_temper, bool const smc, int const thread_id = 0);

    const MapParams &map_params;

    // counts region constraints
    int total_soft_region_constraints;
    int total_soft_plan_constraints;
    int total_soft_constraints; // total constraints
    int num_hard_region_constraints;
    // counts plan constraints
    
    int num_hard_plan_constraints; //

    
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
    std::pair<bool, double> compute_full_split_plan_score(int const num_regions,
                                              PlanVector const &region_ids,
                                              RegionSizes const &region_sizes,
                                              IntPlanAttribute const &region_pops,
                                              int const split_region1,
                                              int const split_region2) const;

    // scores individual regions
    std::pair<bool, double> compute_region_full_score(const Plan &plan, int const region_id
                                                      ) const;
    double compute_region_soft_score(const Plan &plan, int const region_id) const;
    std::pair<bool, double> compute_merged_region_full_score(const Plan &plan,
                                                             int const region1_id,
                                                             int const region2_id
                                                             ) const;

    // scores plans
    // false means probability zero
    // soft score - ie always finite
    std::pair<bool, double> compute_plan_score(const Plan &plan) const;
    std::pair<bool, double> compute_merged_plan_score(const Plan &plan, int const region1_id,
                                                      int const region2_id) const;

    // check if the merged region triggers any hard constraints
    bool merged_region_ok(Plan const &plan, int const region1_id, int const region2_id) const;
    // check any hard entire plan constraints on the merged plan
    bool entire_merged_plan_constraint_only_ok(Plan const &plan, int const region1_id,
                                               int const region2_id) const;
    // check if the entire merged plan is ok
    bool merged_plan_ok(Plan const &plan, int const region1_id, int const region2_id) const;
    // check if the plan satisfies all hard constraints 
    bool satisfies_hard_constraints(Plan const &plan, int const region1_id, int const region2_id) const;
};

#endif
