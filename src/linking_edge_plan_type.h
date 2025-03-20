#include "base_plan_type.h"
#include "weights.h"

class LinkingEdgePlan : public Plan {

    public:
        // constructor
        LinkingEdgePlan(arma::subview_col<arma::uword> region_ids_col, 
                   arma::subview_col<arma::uword> region_sizes_col, 
                   int ndists, int num_regions, const arma::uvec &pop, 
                   bool split_district_only,
                  const std::vector<std::tuple<int, int, double>> &initial_linking_edge_list = {}
                );

        void update_vertex_and_plan_specific_info_from_cut(
            TreeSplitter const &tree_splitter,
            USTSampler &ust_sampler, EdgeCut const cut_edge, 
            const int split_region1_id, const int split_region2_id
        ) override;

        double get_log_eff_boundary_len(
            const MapParams &map_params, const SplittingSchedule &splitting_schedule,
            TreeSplitter const &tree_splitter, 
            const int region1_id, int const region2_id
        ) const override;

};