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


        VertexGraph forest_graph; 
        std::vector<std::tuple<int, int, double>> linking_edges;

        std::unique_ptr<Plan> deep_clone() const override {
            return std::make_unique<LinkingEdgePlan>(*this); // Copy the entire object
        }

        VertexGraph get_forest_adj() override{return forest_graph;};

        void update_vertex_and_plan_specific_info_from_cut(
            TreeSplitter const &tree_splitter,
            USTSampler &ust_sampler, EdgeCut const cut_edge, 
            const int split_region1_id, const int split_region2_id,
            bool const add_region
        ) override;

        double get_log_eff_boundary_len(
            const MapParams &map_params, const SplittingSchedule &splitting_schedule,
            TreeSplitter const &tree_splitter, 
            const int region1_id, int const region2_id
        ) const override;

        std::vector<std::tuple<int, int, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
            const MapParams &map_params, const SplittingSchedule &splitting_schedule,
            TreeSplitter const &tree_splitter,
            std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map = {}
        ) const override;

        // Count the number of valid adj regions in a map
        int count_valid_adj_regions(
            MapParams const &map_params, SplittingSchedule const &splitting_schedule
        ) const override;

        // For a plan and edge across two regions it returns the log probability 
        // that edge would be split in the merged tree across the two regions
        double get_regions_log_splitting_prob(
            TreeSplitter const &tree_splitter, USTSampler &ust_sampler,
            const int region1_root, const int region2_root
        ) const;

        // Get a vector of all valid adj region pairs
        std::vector<std::pair<int,int>> get_valid_adj_regions(
            MapParams const &map_params, SplittingSchedule const &splitting_schedule,
            bool const check_split_constraint = true
        ) const override;

        std::vector<std::array<double, 3>> get_linking_edges() override{
            std::vector<std::array<double, 3>> output;
            output.reserve(linking_edges.size());
            for(auto const& an_edge: linking_edges){
                output.push_back({
                    static_cast<double>(std::get<0>(an_edge)),
                    static_cast<double>(std::get<1>(an_edge)),
                    std::get<2>(an_edge)
                });
            }
            return output;
        };

};