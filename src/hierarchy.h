#ifndef HIERARCHY_H
#define HIERARCHY_H

#include "smc_base.h"

#include <cstdint>

enum class HierarchyMode {
    none = 0,
    connected = 1,
    strict = 2,
    heuristic = 3
};

struct HierarchyUnitCache {
    std::vector<int> vertices;
    std::vector<std::uint64_t> connected_masks;
};

struct HierarchyLevelCache {
    std::vector<int> labels;
    std::vector<int> local_positions;
    std::vector<HierarchyUnitCache> units;
};

struct HierarchySpec {
    HierarchyMode mode = HierarchyMode::none;
    std::vector<uvec> levels;
    std::vector<HierarchyLevelCache> level_caches;
    std::vector<uvec> sampler_levels;
    std::vector<std::vector<int>> sampler_labels;
    std::vector<int> sampler_group_counts;
    std::vector<std::vector<int>> sampler_parents;
    std::vector<int> sampler_finest_adj;
    std::vector<int> sampler_finest_off;

    bool enabled() const {
        return (mode == HierarchyMode::connected ||
                mode == HierarchyMode::strict) &&
               !levels.empty();
    }
};

HierarchySpec hierarchy_spec_from_control(const List &control,
                                          const uvec &fallback_counties,
                                          const Graph &g);

bool hierarchy_region_connected(const Graph &g,
                                const std::vector<int> &region_vertices,
                                const std::vector<char> &region_mark,
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
                              const std::vector<char> &region_mark,
                              const HierarchySpec &hierarchy);

double log_hierarchical_tree_count(const Graph &g,
                                   const uvec &plan,
                                   int district,
                                   const HierarchySpec &hierarchy);

#endif
