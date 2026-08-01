#pragma once

#include <limits>
#include <vector>

enum class HDimerMatchingStatus {
    success,
    invalid_input,
    no_perfect_matching,
    nonplanar,
    numerical_failure
};

// Endpoints are zero-based fragment indices. log_weight must be finite. payload
// is opaque to the matching backend and is returned for a selected edge.
struct HDimerMatchingEdge {
    int u;
    int v;
    double log_weight;
    int payload;
};

struct HDimerMatchingResult {
    HDimerMatchingStatus status = HDimerMatchingStatus::invalid_input;
    // `planar` is meaningful only when planarity_checked is true. The subset
    // dynamic program does not need to test planarity.
    bool planarity_checked = false;
    bool planar = false;
    bool used_subset_dp = false;
    std::vector<int> selected_edge_indices;
    double log_partition = -std::numeric_limits<double>::infinity();
    double max_probability_error = std::numeric_limits<double>::infinity();
    double min_probability = std::numeric_limits<double>::infinity();
};

// Samples a perfect matching with probability proportional to the product of
// edge weights. Parallel edges are treated as distinct payload choices. When
// has_perfect_matching_witness is true, the caller guarantees that the supplied
// edges contain a perfect matching; this permits skipping a separate witness
// search. Every failure, including a numerical one, is returned in `status`.
HDimerMatchingResult
sample_weighted_perfect_matching(int n_vertices, const std::vector<HDimerMatchingEdge> &edges,
                                 bool has_perfect_matching_witness = false);
