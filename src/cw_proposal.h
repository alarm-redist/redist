/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Cycle Walk Proposal
 *
 * The main proposal mechanism that can change district assignments
 * by finding and resampling cycles between adjacent districts.
 ********************************************************/

#ifndef CW_PROPOSAL_H
#define CW_PROPOSAL_H

#include "cw_partition.h"
#include "random.h"
#include <utility>

/*
 * CycleWalkUpdate: Stores information about a proposed update.
 */
struct CycleWalkUpdate {
    DistrictPair changed_districts;
    std::vector<std::pair<int, int>> links;  // Edges to add
    std::vector<std::pair<int, int>> cuts;   // Edges to remove
    bool valid;                               // Whether this is a valid proposal

    CycleWalkUpdate() : valid(false) {}
};

/*
 * Perform one cycle walk step.
 *
 * Algorithm:
 * 1. Pick a random pair of adjacent districts
 * 2. Pick two random boundary edges between them
 * 3. Find the cycle formed by the paths + boundary edges
 * 4. Find valid cut pairs that satisfy population constraints
 * 5. Sample a cut pair and compute MH ratio
 * 6. Accept/reject and update partition
 *
 * Returns: 1 if accepted, 0 if rejected, -1 if no valid proposal
 */
int cycle_walk(LCTPartition& partition,
               double lower, double upper,
               double& accept_ratio);

/*
 * Get a random pair of adjacent districts.
 * Returns true if found, false if no adjacent pairs exist.
 */
bool get_random_adjacent_districts(const LCTPartition& partition,
                                    int& d1, int& d2);

/*
 * Get two random boundary edges between districts d1 and d2.
 * Returns true if found (need at least 2 edges), false otherwise.
 */
bool get_random_edge_pair(const LCTPartition& partition,
                          int d1, int d2,
                          CWEdge& e1, CWEdge& e2);

/*
 * Compute the cycle formed by the two boundary edges.
 * Returns the paths in each district that form the cycle.
 *
 * The cycle is: u1 -> (path in d1) -> u2 -> (edge e2) -> v2 -> (path in d2) -> v1 -> (edge e1) -> u1
 */
bool get_cycle_paths(LCTPartition& partition,
                     const CWEdge& e1, const CWEdge& e2,
                     std::vector<int>& path1, std::vector<int>& path2);

/*
 * Compute cumulative population weights along the cycle.
 * Used to find valid cut positions.
 */
std::vector<int> get_cycle_populations(const LCTPartition& partition,
                                        const std::vector<int>& path1,
                                        const std::vector<int>& path2);

/*
 * Find all valid cut pairs that satisfy population constraints.
 * Returns pairs of indices (i, j) where cutting at positions i and j
 * results in valid district populations.
 */
std::vector<std::pair<int, int>> find_valid_cut_pairs(
    const std::vector<int>& cycle_pops,
    int initial_cut,
    int total_pop,
    double lower, double upper);

/*
 * Apply an update to the partition.
 */
void apply_update(LCTPartition& partition,
                  const CycleWalkUpdate& update);

#endif // CW_PROPOSAL_H
