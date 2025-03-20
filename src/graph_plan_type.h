
#include "base_plan_type.h"
#include "tree_op.h"
#include "weights.h"

class GraphPlan : public Plan {
    // private member variable
    

public:
    using Plan::Plan;
    // implementation of the pure virtual function
    

    void update_vertex_and_plan_specific_info_from_cut(
        TreeSplitter const &tree_splitter,
        USTSampler &ust_sampler, EdgeCut const cut_edge, 
        const int split_region1_id, const int split_region2_id,
        bool const add_region
    ) override;


    std::vector<std::tuple<int, int, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter,
        std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map = {}
    ) const override;


    double get_log_eff_boundary_len(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter, 
        const int region1_id, int const region2_id
    ) const override;

    // Clone method to create a copy of the ForestPlan object
    std::unique_ptr<Plan> deep_clone() const override {
        return std::make_unique<GraphPlan>(*this); // Copy the entire object
    }
};