#ifndef CW_PROPOSAL_H
#define CW_PROPOSAL_H

#include "lct_graph_plan_type.h"
#include "random.h"
#include "redist_types.h"
#include "scoring.h"

#include <utility>
#include <vector>

struct CycleWalkUpdate {
    DistrictPair changed_districts;
    std::vector<std::pair<int, int>> links;
    std::vector<std::pair<int, int>> cuts;
    bool swap_link11;
    bool valid;

    CycleWalkUpdate() : swap_link11(false), valid(false) {}
};

struct CycleWalkDiagnostics {
    int status;         // 1=accept, 0=reject, <0=failure code
    double accept_prob; // MH acceptance probability in [0, 1]
    int cycle_length;
    int n_valid_cuts;

    CycleWalkDiagnostics() : status(0), accept_prob(0.0), cycle_length(0), n_valid_cuts(0) {}
};

/*
 * One cycle-walk step on `plan`.
 *
 * Returns 1 if accepted, 0 if rejected, <0 on a sampling failure:
 *   -1 no adjacent district pair
 *   -2 fewer than 2 boundary edges between picked pair
 *   -3 couldn't construct cycle paths
 *   -4 no valid cut pairs satisfying population bounds
 */
int cycle_walk(LCTGraphPlan &plan, MapParams const &map_params,
               ScoringFunction const &scoring_function, double const compactness,
               RNGState &rng_state, CycleWalkDiagnostics &diagnostics);

bool get_random_adjacent_districts(LCTGraphPlan const &plan, RNGState &rng_state, int &d1,
                                   int &d2);

bool get_random_edge_pair(LCTGraphPlan const &plan, int d1, int d2, RNGState &rng_state,
                          CWEdge &e1, CWEdge &e2);

bool get_cycle_paths(LCTGraphPlan &plan, CWEdge const &e1, CWEdge const &e2,
                     std::vector<int> &path1, std::vector<int> &path2);

std::vector<int> get_collapsed_cycle_weights(LCTGraphPlan &plan, std::vector<int> const &path1,
                                             std::vector<int> const &path2,
                                             arma::uvec const &pop);

// Finds valid (cut1, cut2) pairs. For MMD, the two new regions inherit the
// seat counts of d1 and d2 in some order, so a pair is accepted if EITHER
// assignment of (pop1, pop2) to (d1, d2) falls inside per-region bounds.
// For SMD (s1 == s2) the two bounds collapse and this matches old behavior.
std::vector<std::pair<int, int>> find_valid_cut_pairs(std::vector<int> const &cycle_pops,
                                                      int initial_cut, int total_pop,
                                                      double lower1, double upper1,
                                                      double lower2, double upper2);

void apply_update(LCTGraphPlan &plan, MapParams const &map_params,
                  CycleWalkUpdate const &update);

#endif
