
#include "base_plan_type.h"
#include "tree_op.h"
#include "weights.h"

class GraphPlan : public Plan {
    // private member variable
    

public:
    using Plan::Plan;
    // implementation of the pure virtual function
    

    void update_vertex_info_from_cut(
            Tree &ust, EdgeCut cut_edge, 
            const int split_region1_id, const int split_region2_id,
            bool split_district_only
    ) override;


    std::vector<std::tuple<int, int, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter,
        std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map = {}
    ) const override;

    // Clone method to create a copy of the ForestPlan object
    std::unique_ptr<Plan> deep_clone() const override {
        return std::make_unique<GraphPlan>(*this); // Copy the entire object
    }
};