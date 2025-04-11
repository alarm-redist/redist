#include "graph_plan_type.h"


constexpr bool DEBUG_GRAPH_PLANS_VERBOSE = false; // Compile-time constant

void GraphPlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter,
    USTSampler &ust_sampler, EdgeCut const cut_edge, 
    const int split_region1_id, const int split_region2_id,
    bool const add_region
){
    // Get the root of the tree associated with region 1 and 2
    int split_region1_tree_root, split_region2_tree_root;
    int split_region1_size, split_region2_size;
    int split_region1_pop, split_region2_pop;

    cut_edge.get_split_regions_info(
        split_region1_tree_root, split_region1_size, split_region1_pop,
        split_region2_tree_root, split_region2_size, split_region2_pop
    );

    // update the vertex labels
    assign_region_id_from_tree(ust_sampler.ust, region_ids,
        split_region1_tree_root, split_region1_id);

    assign_region_id_from_tree(ust_sampler.ust, region_ids,
        split_region2_tree_root, split_region2_id);

    return;
}





/*  
 * @title Get Valid Adjacent Regions to Boundary Lengths Map
 * 
 * Returns a unordered_map mapping valid pairs of adjacent regions to 
 * the length of the boundary between them in the graph
 * 
 * Finds all valid pairs of adjacent regions (meaning their sizes are valid 
 * to merge) and returns a hash map mapping the pairs of valid adjacent 
 * regions to the length of the boundary between them in `g`
 * 
 * 
 * @param g A graph (adjacency list) passed by reference
 * @param plan A plan object
 * @param check_adj_to_regions A vector tracking whether or not we should 
 * check for edges adjacent to vertices in a region of a particular size. For
 * example, `check_adj_to_regions[i] == true` means we will attempt to find 
 * edges adjacent to any vertex in a region of size i. This vector is 1-indexed
 * meaning we don't need to subtract region size by 1 when accessing.
 * @param valid_merge_pairs A 2D `ndists+1` by `ndists+1` boolean matrix that
 * uses region size to check whether or not two regions are considered a valid
 * merge that can be counted in the map. For example `valid_merge_pairs[i][j]`
 * being true means that any regions where the sizes are (i,j) are considered
 * valid to merge. 2D matrix is 1 indexed (don't need to subtract region size)
 * @param existing_pair_map A hash map mapping pairs of regions to a double. 
 * This is an optional parameter and if its empty then its ignored. If its not
 * empty then pairs already in the hash won't be counted added to the output.
 * 
 * @details No modifications to inputs made
 * @return A hash map mapping (std::pair<int,int> to int) that maps a pair of
 * region ids in the form (smaller id, bigger id) to the length of the boundary
 * between them in `g`.
 */
std::unordered_map<std::pair<int, int>, int, bounded_hash> NEW_get_valid_pairs_and_graph_boundary_len_map(
    Graph const &g, Plan const &plan,
    std::vector<bool> const &check_adj_to_regions,
    std::vector<std::vector<bool>> const &valid_merge_pairs,
    std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map = {}
){
    // In FUTURE PROBABLY REMOVE THE EXISTING PAIR, NOT WORTH THE HASSLE WHEN 
    // JUST INCREMENTING, ITS REALLY ONLY VALUABLE FOR THE FOREST STUFF TBH

    // NOTE: In the case where its one district split you maybe don't even need
    // a hash map since you can just index by the not remainder region but nbd
    // for now
    
    // Initialize unordered_map with num_region * 2.5 buckets
    // Hueristic. Bc we know planar graph has at most 3|V| - 6 edges
    int init_bucket_size = std::ceil(2.5*plan.num_regions);

    // create the hash map
    std::unordered_map<std::pair<int, int>, int, bounded_hash> region_pair_map(
        init_bucket_size, bounded_hash(plan.num_regions-1)
        );

    int const V = g.size();
    // boolean for whether or not we should 
    bool const check_existing_map = existing_pair_map.size() != 0;


    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = plan.region_ids(v);
        auto v_region_size = plan.region_sizes(v_region_num);

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!check_adj_to_regions.at(v_region_size)){
            continue;
        }

        // get neighbors
        std::vector<int> nbors = g.at(v);

        // now iterate over its neighbors
        for (int v_nbor : nbors) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = plan.region_ids(v_nbor);

            // ignore if they are in the same region
            if(v_region_num == v_nbor_region_num) continue;
            // ignore if pair can't be merged
            auto v_nbor_region_size = plan.region_sizes(v_nbor_region_num);
            if(!valid_merge_pairs[v_region_size][v_nbor_region_size]) continue;

            // else they are different so regions are adj
            // check if we need to be aware of double counting. IE if both regions
            // are ones where check adjacent is true then its possible to double count edges
            bool double_counting_not_possible = !check_adj_to_regions[v_nbor_region_size];
            // Rprintf("Neighbor size is %d and Double counting is %d\n", 
            //     (int) v_nbor_region_size, (int) double_counting_not_possible);

            // Now add this edge if either we can't double count it 
            // or `v_region_num` is the smaller of the pair 
            if(double_counting_not_possible || v_region_num < v_nbor_region_num){
                // Rprintf("Counting (%d,%d)\n", v, v_nbor);
                int smaller_region_id  = std::min(v_nbor_region_num,v_region_num);
                int bigger_region_id  = std::max(v_nbor_region_num,v_region_num);

                // if we need to see if the pair is already in the hash and its not then ignore
                if(check_existing_map && 
                    existing_pair_map.find({smaller_region_id, bigger_region_id}) == existing_pair_map.end()){
                        Rprintf("Skipping!!\n");
                        continue;
                    }
                region_pair_map[{smaller_region_id, bigger_region_id}]++;
            }
        }
    }

    // adj_pairs_vec
    return region_pair_map;
}


// TEMP FIX
std::vector<std::unordered_set<std::pair<int, int>, bounded_hash>> TEMP_inner_county_boundaries_crossed(
    MapParams const &map_params, Plan const &plan
){
    std::vector<std::unordered_set<std::pair<int, int>, bounded_hash>> county_crossings_sets;
    county_crossings_sets.reserve(plan.num_regions);

    for (size_t i = 0; i < plan.num_regions; i++)
    {
        // Initialize unordered_map with num_counties/region size * 2.5 buckets
        // Hueristic. Bc we know planar graph has at most 3|V| - 6 edges
        int init_bucket_size = std::ceil((2 * map_params.num_counties) / plan.region_sizes(i));
        county_crossings_sets.emplace_back(init_bucket_size, bounded_hash(map_params.num_counties));
    }
    
    

    for (int v = 0; v < map_params.V; v++) {
        // Find out which region this vertex corresponds to
        int v_region = plan.region_ids(v);
        int v_county = map_params.counties(v);

        // now iterate over its neighbors
        for (int u : map_params.g[v]) {
            // find which region neighbor corresponds to
            int u_region = plan.region_ids(u);

            // ignore if they are not in the same region 
            if(v_region != u_region) continue;

            int u_county = map_params.counties(u);

            // if the edge crosses a county boundary then hash it
            if(v_county != u_county){
                // pair is always (smaller, bigger)
                if(v_county < u_county){
                    county_crossings_sets[v_region].insert({v_county, u_county});
                }else{
                    county_crossings_sets[v_region].insert({u_county, v_county});
                }
            }
        }
    }

    return county_crossings_sets;
}


// temp fix 
// given a pair of adjacent regions it county all illegal county boundaries
// meaning boundaries between counties where either region already has an internal 
// edge which crosses that county boundary
int count_illegal_county_crossing_edges(
    MapParams const &map_params, Plan const &plan, 
    int const region1_id, int const region2_id,
    std::vector<std::unordered_set<std::pair<int, int>, bounded_hash>> county_crossings_sets
){
    if(map_params.num_counties <= 1) return 0;
    // REprintf("\nHey! {");
    // for(const auto &el: split_counties){
    //     REprintf("%d, ", el);
    // }
    // REprintf("}\n");
    int remove_count = 0;
    // Now count boundaries between two regions that 
    for (int v = 0; v < map_params.V; v++)
    {
        int v_region = plan.region_ids(v);
        int v_county = map_params.counties(v);
        // ignore if not pairs 
        if(v_region != region1_id) continue;

        for(auto const &u: map_params.g[v]){
            // ignore if not regions we want
            int u_region = plan.region_ids(u);
            if(u_region != region2_id) continue;
            int u_county = map_params.counties(u);
            // ignore if same county 
            if(v_county == u_county) continue;
            bool v_county_bigger = v_county > u_county;
            std::pair<int, int> county_pair = v_county_bigger ? std::make_pair(u_county, v_county) : std::make_pair(v_county, u_county);
            // if either region already has that pair then count 
            if(county_crossings_sets[v_region].count(county_pair) > 0 || 
               county_crossings_sets[u_region].count(county_pair) > 0){
                    remove_count++;
               }
        }
    }

    return remove_count;
}

std::vector<std::tuple<int, int, double>> GraphPlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter,
    std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map
) const{
    // get pairs with graph theoretic boundary 
    auto valid_adj_region_pairs_to_boundary_map = NEW_get_valid_pairs_and_graph_boundary_len_map(
        map_params.g, *this,
        splitting_schedule.check_adj_to_regions,
        splitting_schedule.valid_merge_pair_sizes
    );

    // declare the vector of tuples
    std::vector<std::tuple<int, int, double>> region_pairs_tuple_vec;
    region_pairs_tuple_vec.reserve(
        valid_adj_region_pairs_to_boundary_map.size()
    );

    if(DEBUG_GRAPH_PLANS_VERBOSE){
        REprintf("Size is %d with %d counties\n", 
            valid_adj_region_pairs_to_boundary_map.size(),
            map_params.num_counties);
    }

    std::vector<bool> visited(map_params.V);


    auto county_crossings_sets = TEMP_inner_county_boundaries_crossed(
        map_params, *this
    );

    for(auto const &a_pair: valid_adj_region_pairs_to_boundary_map){
        // see if illegal merge 
        int merged_plan_county_splits = count_merged_county_splits(
            map_params, visited,
            a_pair.first.first, a_pair.first.second);

        if(merged_plan_county_splits > num_regions - 2){
            if(DEBUG_GRAPH_PLANS_VERBOSE){
                REprintf("(%d,%d):%d regions when merged but merged has %d splits!\n", 
                    a_pair.first.first, a_pair.first.second,
                num_regions-1, merged_plan_county_splits);
            }
            continue;
        }

        // remove county correction 
        int correction_term = count_illegal_county_crossing_edges(
            map_params, *this, 
            a_pair.first.first, a_pair.first.second,
            county_crossings_sets
        );
        if(DEBUG_GRAPH_PLANS_VERBOSE){
            REprintf("(%d,%d): %d-%d\n",
                a_pair.first.first, a_pair.first.second, a_pair.second,
                correction_term
            );
        }

        // double testy = log_graph_boundary(map_params.g, region_ids,
        //     a_pair.first.first, a_pair.first.second,
        //      map_params.num_counties, map_params.counties
        // );
        // double diff = std::fabs(
        //     std::exp(testy) - (a_pair.second - correction_term)
        // );

        // if(diff > 1e10){
        //     REprintf("Uh oh %f vs %d!\n",
        //         std::exp(testy), (a_pair.second - correction_term)
        //     );
        // }

        // skip if not edges
        if(correction_term == a_pair.second){
            continue;
        } 
        
        // add the region ids and the log of the boundary 
        region_pairs_tuple_vec.emplace_back(
            a_pair.first.first, 
            a_pair.first.second,
            std::log(static_cast<double>(a_pair.second-correction_term))
        );
    }

    return region_pairs_tuple_vec;

}


double GraphPlan::get_log_eff_boundary_len(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, 
    const int region1_id, int const region2_id
) const{
    // Return the log of the graph theoretic boundary 
    return log_graph_boundary(map_params.g, region_ids,
        region1_id, region2_id, map_params.num_counties, map_params.counties);
}