#ifndef CW_FOREST_WALK_H
#define CW_FOREST_WALK_H

#include "lct_graph_plan_type.h"
#include "random.h"
#include "advanced_types.h"

/*
 * One internal forest-walk step. Reshuffles the spanning tree within a
 * single region without changing region assignments. Returns 0 on success,
 * 1 if no internal edge could be found.
 */
int internal_forest_walk(LCTGraphPlan &plan, MapParams const &map_params, RNGState &rng_state,
                         int max_attempts = 100);

bool get_random_internal_edge(LCTGraphPlan &plan, MapParams const &map_params,
                              RNGState &rng_state, int &u, int &v, int max_attempts = 100);

#endif
