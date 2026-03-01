/********************************************************
 * BUD MCMC Redistricting Sampler
 * BUD Proposal Mechanism
 *
 * Core proposal that operates on multi-district cycles
 * using the Balanced Up-Down Walk algorithm.
 * Closely follows the Julia implementation in BUD.jl.
 ********************************************************/

#ifndef BUD_PROPOSAL_H
#define BUD_PROPOSAL_H

#include "bud_partition.h"
#include "mcmc_gibbs.h"
#include "random.h"
#include <kirchhoff_inline.h>

/*
 * Update struct matching Julia's Update{T}.
 * Stores information needed to accept a proposal.
 */
struct BUDUpdate {
    std::vector<DistrictPair> district_path; // old district path edges
    int link_u, link_v;                      // the non-tree edge to link
    int cut_ind;                             // which cut is the actual cut (0-based)
    std::vector<std::pair<int, int>> cuts;   // all cuts in cycle order
    std::vector<int> cycle_path;             // cycle path vertex indices
    std::map<int, int> fragment_map;         // vertex → cycle_path vertex (for all affected)
    bool valid;

    BUDUpdate() : link_u(-1), link_v(-1), cut_ind(-1), valid(false) {}
};

/*
 * Result from getBalancedCuts.
 */
struct BalancedCutsResult {
    std::vector<int> edge_inds;      // balanced cut positions (0-based into cycle_path)
    std::vector<double> pathWeights; // cumulative inverse edge weights
    double cumWeight;
};

/*
 * Perform one full BUD step: proposal + MH accept/reject with energy.
 * Returns: 1 if accepted, 0 if rejected, <0 if no valid proposal
 */
int bud_step(BUDPartition& partition,
             double lower, double upper,
             double target,
             double compactness,
             const arma::uvec& counties,
             Rcpp::List constraints,
             double& accept_prob_out);

/*
 * Apply an accepted BUD update to the partition.
 * Follows Julia's update_partition!.
 */
void apply_bud_update(BUDPartition& partition, const BUDUpdate& update);

#endif // BUD_PROPOSAL_H
