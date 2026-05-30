#pragma once
#ifndef LCT_GRAPH_PLAN_TYPE_H
#define LCT_GRAPH_PLAN_TYPE_H

#include "cw_lct.h"
#include "graph_plan_type.h"

#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>

/*
 * Cyclewalk-flavored GraphPlan: same region_ids/region_sizes/region_pops
 * machinery, with a Link-Cut Tree representation of a spanning forest
 * (one tree per region) layered on top.
 *
 * Reuses GraphPlan's implementations of Plan's pure virtuals; only
 * update_vertex_and_plan_specific_info_from_cut is overridden so the LCT
 * stays in sync when the plan is mutated by a split (i.e. when used as an
 * SMC step in the future).
 */

struct CWEdge {
    int u, v;
    double weight;

    CWEdge(int a, int b, double w = 1.0);
    bool operator<(const CWEdge &other) const;
    bool operator==(const CWEdge &other) const;
};

using DistrictPair = std::pair<int, int>;
using EdgeSet = std::set<CWEdge>;
using CrossEdgeMap = std::map<DistrictPair, EdgeSet>;

struct EdgePairHash {
    std::size_t operator()(const std::pair<int, int> &p) const {
        int u = std::min(p.first, p.second);
        int v = std::max(p.first, p.second);
        return std::hash<int>()(u) ^ (std::hash<int>()(v) << 1);
    }
};

struct EdgePairEqual {
    bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const {
        return (a.first == b.first && a.second == b.second) ||
               (a.first == b.second && a.second == b.first);
    }
};

using EdgeWeightMap =
    std::unordered_map<std::pair<int, int>, double, EdgePairHash, EdgePairEqual>;

class USTSampler;

class LCTGraphPlan : public GraphPlan {
  public:
    using GraphPlan::GraphPlan;

    // Spanning forest state — populated by init_lct_from_regions, not by the
    // Plan/PlanEnsemble constructors.
    std::unique_ptr<LinkCutTree> lct;
    std::vector<int> district_roots;
    CrossEdgeMap cross_edges;
    EdgeWeightMap edge_weights;
    double default_edge_weight = 1.0;

    // Build a spanning tree per region via USTSampler and populate LCT, roots,
    // and cross-region edge index. region_ids/region_sizes/region_pops must
    // already be set (this is what PlanEnsemble's loading ctor does).
    // Returns true on success.
    bool init_lct_from_regions(MapParams const &map_params, USTSampler &ust_sampler,
                               RNGState &rng_state, int const max_attempts_per_region = 25);

    // R-side custom edge weights: list of list(edge=c(u,v), weight=w), 1-indexed.
    void set_edge_weights(MapParams const &map_params, Rcpp::List const &edge_weights_list);
    double get_edge_weight(int u, int v) const;

    // Cross-edge queries.
    const EdgeSet &get_cross_edges(int d1, int d2) const;
    bool districts_adjacent(int d1, int d2) const;
    std::vector<DistrictPair> get_adjacent_district_pairs() const;

    // Recompute the cross_edges map from current region_ids. Cheap (O(E))
    // and idempotent; used after a proposal accept/reject.
    void rebuild_cross_edges(Graph const &g);

    // For inspection / debug.
    void print_state(int verbosity = 1) const;

    // SMC-step hook: after a split applied to this plan, update LCT to
    // reflect the new region assignments. Falls through to GraphPlan for
    // region_ids bookkeeping. Not exercised by the MCMC inner loop.
    void update_vertex_and_plan_specific_info_from_cut(
        TreeSplitter const &tree_splitter, USTSampler &ust_sampler,
        EdgeCut const cut_edge, int const split_region1_id, int const split_region2_id,
        bool const add_region) override;

  private:
    // Helper: load `tree` (rooted at `root`) into the LCT under `district`.
    void load_tree_into_lct(Tree const &tree, int root, int district);

    // Empty set returned by get_cross_edges when no adjacency exists.
    static const EdgeSet empty_edge_set;
};

#endif
