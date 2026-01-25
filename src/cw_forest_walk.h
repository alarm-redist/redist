/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Internal Forest Walk Proposal
 *
 * A simpler proposal that shuffles spanning trees within districts
 * without changing district assignments. Helps with mixing.
 ********************************************************/

#ifndef CW_FOREST_WALK_H
#define CW_FOREST_WALK_H

#include "cw_partition.h"
#include "random.h"

/*
 * Perform one internal forest walk step.
 *
 * Algorithm:
 * 1. Pick a random internal edge (both endpoints in same district)
 * 2. Find the path between endpoints in the spanning tree
 * 3. Sample an edge from the path (weighted by 1/edge_weight)
 * 4. If we sample the new edge, do nothing
 * 5. Otherwise, cut the sampled edge and link the new edge
 *
 * This changes the spanning tree but not district assignments.
 * It's a reversible move that helps explore the space of spanning trees.
 *
 * Returns: 0 on success, 1 on failure (no internal edges found)
 */
int internal_forest_walk(LCTPartition& partition, int max_attempts = 100);

/*
 * Get a random internal edge (endpoints in same district).
 * Returns true if found, false if no internal edge exists.
 */
bool get_random_internal_edge(const LCTPartition& partition,
                               int& u, int& v,
                               int max_attempts = 100);

#endif // CW_FOREST_WALK_H
