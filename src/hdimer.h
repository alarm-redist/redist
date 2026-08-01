#ifndef HDIMER_H
#define HDIMER_H

#include "redist_types.h"
#include <cstdint>
#include <vector>

struct HDimerForestEdge {
    int u;
    int v;
};

struct HDimerStepDiagnostics {
    bool changed = false;
    bool planar = true;
    bool numerical_failure = false;
    bool no_matching = false;
    bool constraint_rejected = false;
    bool used_subset_dp = false;
    int matching_edges_changed = 0;
    int vertices_changed = 0;
    double old_matching_probability = NA_REAL;
    double cut_seconds = 0.0;
    double fragment_seconds = 0.0;
    double candidate_seconds = 0.0;
    double matching_seconds = 0.0;
    double relabel_seconds = 0.0;
};

bool sample_hierarchical_forest(const Graph &graph, const arma::uvec &plan, int n_districts,
                                const arma::umat &hierarchy,
                                std::vector<HDimerForestEdge> &forest);

#endif
