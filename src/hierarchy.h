#ifndef HIERARCHY_H
#define HIERARCHY_H

#include "smc_base.h"

enum class HierarchyMode {
    none = 0,
    speedup = 1,
    strict = 2
};

struct HierarchySpec {
    HierarchyMode mode = HierarchyMode::none;
    std::vector<uvec> levels;
    std::vector<uvec> sampler_levels;

    bool enabled() const {
        return mode != HierarchyMode::none && !levels.empty();
    }
};

HierarchySpec hierarchy_spec_from_control(const List &control,
                                          const uvec &fallback_counties);

bool hierarchy_region_connected(const Graph &g,
                                const std::vector<int> &region_vertices,
                                const HierarchySpec &hierarchy);

bool hierarchy_plan_connected(const Graph &g,
                              const uvec &plan,
                              const HierarchySpec &hierarchy);

bool hierarchy_plan_valid(const Graph &g,
                          const uvec &plan,
                          const HierarchySpec &hierarchy);

int admissible_boundary_count(const Graph &g,
                              const uvec &plan,
                              int first_label,
                              const std::vector<int> &other_labels,
                              const std::vector<int> &region_vertices,
                              const HierarchySpec &hierarchy);

double log_hierarchical_tree_count(const Graph &g,
                                   const uvec &plan,
                                   int district,
                                   const HierarchySpec &hierarchy);

#endif
