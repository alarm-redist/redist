#include "wilson.h"

namespace{

/*
 * Builds a deterministic spanning tree on a county using depth first search
 */
static void add_county_to_tree_dfs(
    Graph const &county_restricted_g,
    int const county_id,
    std::vector<int> const &county_vertices, // vector of vertices in the county 
    std::vector<bool> &visited,
    DummyTreeQueue &queue,
    Tree &ust
) {
    int const V = static_cast<int>(county_restricted_g.size());
    int const n_vtx = static_cast<int>(county_vertices.size());

    // Optional checks 
    if constexpr(perf_config::supposedly_safe_input_checks){
        if (n_vtx <= 0) {
            return;
        }
    }

    // set the root as the first county vertex 
    int root = -1;

    // Walk through the vertices in the county 
    for (int const v : county_vertices) {
        // optional bounds checking 
        if constexpr (perf_config::bounds_checking){
            if (v < 0 || v >= V) {
                std::ostringstream oss;
                oss << "add_county_dfs_tree_edges got invalid county vertex.\n";
                oss << "county_id=" << county_id << "\n";
                oss << "v=" << v << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        // Check if we have already visited this vertex in the Wilson call (sample_sub_ust)
        // If yes then that means this is the vertex to use as the root for this county's portion
        // of the tree. 
        if (visited[v]) {
            // Checks we haven't already found a root 
            if constexpr(perf_config::supposedly_safe_input_checks){
                // This shouldn't be neccesary if the tree is truly a directed tree
                if (root != -1) {
                    std::ostringstream oss;
                    oss << "add_county_dfs_tree_edges found multiple already-visited "
                        << "vertices in county.\n";
                    oss << "county_id=" << county_id << "\n";
                    oss << "first_root=" << root << "\n";
                    oss << "second_root=" << v << "\n";
                    throw std::runtime_error(oss.str());
                }
            }
            // we've found the vertex to use as the root so we can exit if not error checking
            root = v;
            // if not error checking then break out now since the tree is supposed to be
            // a directed tree 
            if constexpr(!perf_config::supposedly_safe_input_checks){
                break;
            }
        }

    }

    if constexpr(perf_config::bounds_checking){
        if (root < 0) {
            std::ostringstream oss;
            oss << "add_county_dfs_tree_edges could not find already-visited county root.\n";
            oss << "county_id=" << county_id << "\n";
            oss << "county_vertices.size()=" << county_vertices.size() << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    // clear the stack and start from the root
    queue.clear();
    queue.push(root);

    int seen_count = 1;

    // perform DFS
    while (!queue.empty()) {
        int const v = queue.pop();

        if constexpr(perf_config::bounds_checking){
            if (v < 0 || v >= V) {
                std::ostringstream oss;
                oss << "add_county_dfs_tree_edges popped invalid vertex.\n";
                oss << "county_id=" << county_id << "\n";
                oss << "v=" << v << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }


        for (auto const u : county_restricted_g[v]) {
            if constexpr(perf_config::bounds_checking){
                if (u < 0 || u >= V) {
                    std::ostringstream oss;
                    oss << "add_county_dfs_tree_edges saw invalid graph neighbor.\n";
                    oss << "county_id=" << county_id << "\n";
                    oss << "v=" << v << "\n";
                    oss << "u=" << u << "\n";
                    oss << "V=" << V << "\n";
                    throw std::runtime_error(oss.str());
                }
            }

            // because this graph only contains within county edges we don't need to check
            // county status

            // ignore if we've already visited this vertex since it was added to the stack
            if (visited[u]) {
                continue;
            }

            // This is a real graph edge v--u, oriented away from the county root.
            ust[v].push_back(u);
            // mark as visited and add this to stack
            visited[u] = true;
            queue.push(u);

            if constexpr(perf_config::redundancy_checks) ++seen_count;
        }
    }

    // This just checks we visited the number of vertices we expected to visit 
    if constexpr(perf_config::redundancy_checks){
        if (seen_count != n_vtx) {
            std::ostringstream oss;
            oss << "add_county_dfs_tree_edges could not span county-induced subgraph.\n";
            oss << "county_id=" << county_id << "\n";
            oss << "root=" << root << "\n";
            oss << "seen_count=" << seen_count << "\n";
            oss << "n_vtx=" << n_vtx << "\n";

            oss << "Unvisited county vertices: ";
            for (int const v : county_vertices) {
                if (!visited[v]) {
                    oss << v << " ";
                }
            }
            oss << "\n";

            throw std::runtime_error(oss.str());
        }
    }

}


/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
int walk_until(const Graph &g, int root, std::vector<int> &path, int MAX,
               const std::vector<bool> &visited, const std::vector<bool> &ignore,
               const uvec &counties, RNGState &rng_state);

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
void loop_erase(std::vector<int> &path, int proposal);

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
void walk_until_cty(const Multigraph &mg, int root, std::vector<std::array<int, 3>> &path,
                    const std::vector<bool> &visited, const std::vector<bool> &ignore,
                    RNGState &rng_state);

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase_cty(std::vector<std::array<int, 3>> &path, int proposal, int root);

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
// TESTED
int sample_sub_ust(MapParams const &map_params, Tree &tree, int &root, double const lower,
                   double const upper, std::vector<bool> &visited,
                   const std::vector<bool> &ignore, Tree &cty_tree, TreePopStack &county_stack,
                   DummyTreeQueue &dummy_county_tree_queue,
                   arma::uvec &county_pop, std::vector<std::vector<int>> &county_members,
                   std::vector<bool> &c_visited, std::vector<int> &cty_pop_below,
                   std::vector<std::array<int, 3>> &county_path, std::vector<int> &path,
                   RNGState &rng_state, bool const draw_fake_dummy_trees) {
    // auto t1_start = std::chrono::steady_clock::now();
    int const n_county = map_params.num_counties;
    std::fill(c_visited.begin(), c_visited.end(), true);
    // reset county pops to zero
    county_pop.zeros();
    int tot_pop = 0;
    // reset the county members inner vectors
    for (size_t i = 0; i < map_params.num_counties; i++) {
        county_members[i].clear();
    }

    int remaining = 0;
    int const V = map_params.V;
    for (int i = 0; i < V; i++) {
        if (ignore[i]) {
            visited[i] = true;
        } else {
            visited[i] = false;
            remaining++;
            int county = map_params.counties[i] - 1;
            tot_pop += map_params.pop[i];
            county_pop(county) += map_params.pop[i];
            if (c_visited[county]) {
                c_visited[county] = false;
            }
            county_members[county].push_back(i);
        }
    }

    int c_remaining = 0;
    for (int i = 0; i < n_county; i++) {
        c_remaining += 1 - c_visited[i];
    }

    // pick root
    int lower_i = 0;
    int lower_c = 0;
    root = rvtx(visited, V, remaining, lower_i, rng_state);
    visited[root] = true;
    remaining--;
    c_visited[map_params.counties[root] - 1] = true;
    c_remaining--;

    // Connect counties
    // clear the tree
    clear_tree(cty_tree);
    county_path.clear();
    while (c_remaining > 0) {
        int add = rvtx(c_visited, n_county, c_remaining, lower_c, rng_state);
        // random walk from `add` until we hit the path
        walk_until_cty(map_params.cg, add, county_path, c_visited, ignore, rng_state);
        // update visited list and constructed tree
        int added = county_path.size();
        if (added == 0) { // bail
            return 1;
        }
        c_remaining -= added;
        c_visited[add] = true;
        for (int i = 0; i < added; i++) {
            c_visited[county_path[i][0]] = true;
            // reverse path so that arrows point away from root
            tree[county_path[i][2]].push_back(county_path[i][1]);
            cty_tree[county_path[i][0]]
                .push_back(map_params.counties(county_path[i][1]) - 1);

            visited[county_path[i][1]] = true; // root for next district
            remaining--;
        }
    }

    // figure out which counties will not need to be split
    if (n_county > 1) {
        // don't need to fill pop below since it gets reset
        get_tree_pops_below(cty_tree, map_params.counties[root] - 1, county_stack, county_pop,
                            cty_pop_below);
        for (int i = 0; i < n_county; i++) {
            int n_vtx = county_members[i].size();
            if (n_vtx <= 1)
                continue;
            // check child counties
            int children = cty_tree[i].size();
            int split_ub = cty_pop_below[i];
            int split_lb = split_ub - county_pop[i];
            if (lower - 1 < county_pop[i])
                split_lb = (int)lower;
            for (int j = 0; j < children; j++) {
                int pop_child = cty_pop_below[cty_tree[i][j]];
                if (pop_child >= 0 && pop_child < split_lb) {
                    split_lb = pop_child;
                }
            }

            // split_lb < split_ub so the smallest possible population is
            // min(split_lb, tot_pop - split_ub)
            // its impossible to split if smallest possible size is bigger than largest ub
            // bool miss_first = std::min(split_lb, tot_pop - split_ub) > upper;
            // // biggest possible population is max(split_ub, total_pop - split_lb)
            // // its impossible to split if largest possible size is smaller than smallest lb
            // bool miss_second = std::max(split_ub, tot_pop - split_lb) < lower;

            bool miss_first = split_ub < lower || split_lb > upper;
            bool miss_second = (tot_pop - split_lb) < lower || (tot_pop - split_ub) > upper;

            // impossible for this county to need to be split
            // so we fill in a dummy tree 
            if (cty_pop_below[i] >= 0 && (miss_first && miss_second)) {
                // check if we need the dummy tree to actually be a subset of g
                if(!draw_fake_dummy_trees){
                    add_county_to_tree_dfs(
                        map_params.county_restricted_graph,
                        map_params.counties[county_members[i][0]], // pass in the county CAREFUL COUNTIES ARE 1-indexed
                        county_members[i],
                        visited,
                        dummy_county_tree_queue,
                        tree
                    );
                    remaining -= n_vtx - 1; // already visited county root
                }else{
                    // fill in with a fake dummy tree
                    remaining -= n_vtx - 1; // already visited county root
                    int cty_root = -1;
                    for (int j = 0; j < n_vtx; j++) {
                        int vtx_idx = county_members[i][j];
                        if (visited[vtx_idx]) { // county root
                            cty_root = j;
                        }
                        if (j > 0 && j != cty_root + 1) {
                            tree[vtx_idx].push_back(county_members[i][j - 1]);
                        }
                        visited[vtx_idx] = true;
                    }

                    if (cty_root < n_vtx - 1) {
                        tree[county_members[i][cty_root]]
                            .push_back(county_members[i][n_vtx - 1]);
                    }
                }
            }
        }
    }

    // Generate tree within each county
    if (remaining > 0) {
        path.clear();
        path.resize(remaining + 2);
        int max_try = 50 * remaining * (static_cast<int>(std::log(remaining)) + 1);
        while (remaining > 0) {
            int add = rvtx(visited, V, remaining, lower_i, rng_state);
            // random walk from `add` until we hit the path
            int added = walk_until(map_params.g, add, path, max_try, visited, ignore,
                                   map_params.counties, rng_state);
            // update visited list and constructed tree
            if (added == 0) { // bail
                return 1;
            }
            remaining -= added - 1; // minus 1 because ending vertex already in tree
            for (int i = 0; i < added - 1; i++) {
                visited[path[i]] = true;
                // reverse path so that arrows point away from root
                tree[path[i + 1]].push_back(path[i]);
            }
        }
    }

    return 0;
}

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
int walk_until(const Graph &g, int root, std::vector<int> &path, int MAX,
               const std::vector<bool> &visited, const std::vector<bool> &ignore,
               const uvec &counties, RNGState &rng_state) {
    path[0] = root;
    // walk until we hit something in `visited`
    int curr = root;
    int county = counties[root];
    int added = 1; // cursor
    int i;
    for (i = 0; i < MAX; i++) {
        int proposal = rnbor(g, curr, rng_state);
        if (ignore[proposal] || counties[proposal] != county) {
            continue;
        } else if (!visited[proposal]) {
            for (int j = added - 1; j >= 0; j--) {
                if (path[j] == proposal) { // if yes, restart from there
                    added = j;
                    break;
                }
            }
            path[added++] = proposal;
        } else { // reached something in `visited`
            path[added++] = proposal;
            break;
        }
        curr = proposal;
    }
    if (i == MAX) {
        added = 0;
    }

    return added;
}

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
void walk_until_cty(const Multigraph &mg, int root, std::vector<std::array<int, 3>> &path,
                    const std::vector<bool> &visited, const std::vector<bool> &ignore,
                    RNGState &rng_state) {
    path.clear();

    // walk until we hit something in `visited`
    int curr = root;
    // while (true) {
    int i;
    int max = visited.size() * 500;
    for (i = 0; i < max; i++) {
        int prop_idx = rng_state.r_int((int)mg[curr].size());
        int proposal = mg[curr][prop_idx][0];
        if (ignore[mg[curr][prop_idx][2]] || ignore[mg[curr][prop_idx][1]]) {
            continue;
        } else if (!visited[proposal]) {
            path.push_back(mg[curr][prop_idx]);
            loop_erase_cty(path, proposal, root);
        } else {
            path.push_back(mg[curr][prop_idx]);
            break;
        }
        curr = proposal;
    }
    if (i == max) {
        path.clear();
    }
}

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase_cty(std::vector<std::array<int, 3>> &path, int proposal, int root) {
    int length = path.size();
    if (proposal == root) {
        path.erase(path.begin(), path.begin() + length);
        return;
    }

    int idx;
    for (idx = 0; idx < length - 1; idx++) {
        if (path[idx][0] == proposal)
            break;
    }

    if (idx != length - 1) { // a loop
        path.erase(path.begin() + idx + 1, path.begin() + length);
    }
}

}

// [[Rcpp::export]]
Tree sample_ust(List l, const arma::uvec &pop, double lower, double upper,
                const arma::uvec &counties, const std::vector<bool> ignore) {
    RNGState rng_state((int)Rcpp::sample(INT_MAX, 1)[0]);
    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();

    int FAKE_NDISTS = 6;
    double FAKE_TARGET = 6.6;

    MapParams map_params(l, counties, pop, FAKE_NDISTS, FAKE_NDISTS, std::vector<int>{}, lower,
                         FAKE_TARGET, upper, SamplingSpace::GraphSpace);

    Tree tree = init_tree(V);
    int root;
    std::vector<bool> visited(V);
    Tree county_tree = init_tree(map_params.num_counties);
    TreePopStack county_stack(map_params.num_counties + 1);
    DummyTreeQueue dummy_county_tree_queue(map_params.V);
    arma::uvec county_pop(map_params.num_counties, arma::fill::zeros);
    std::vector<std::vector<int>> county_members(map_params.num_counties, std::vector<int>{});
    std::vector<bool> c_visited(map_params.num_counties, true);
    std::vector<int> cty_pop_below(map_params.num_counties, 0);
    std::vector<std::array<int, 3>> county_path;
    std::vector<int> path;

    sample_sub_ust(map_params, tree, root, lower, upper, visited, ignore, county_tree,
                   county_stack, dummy_county_tree_queue, 
                   county_pop, county_members, c_visited, cty_pop_below,
                   county_path, path, rng_state, false);
    return tree;
}





/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2025/3
 * Purpose: Encapsulation of uniform spanning tree sampler functions
 ********************************************************/


bool USTSampler::draw_ust(
      int &root, double const lower, double const upper,
        RNGState &rng_state) {
    // We assume that ignore has already been properly set 
    // Now get a uniform spanning tree drawn on the subgraph denoted by the ignore
    // vertices 
    int result = sample_sub_ust(map_params, ust, root, lower, upper, visited, ignore,
                                county_tree, county_stack, dummy_county_tree_queue,
                                county_pop, county_members,
                                c_visited, cty_pop_below, county_path, path, rng_state,
                                fake_dummy_trees_ok);
    // result == 0 means it was successful
    return (result == 0);
}


std::pair<bool, int> USTSampler::attempt_to_draw_tree_on_region(RNGState &rng_state, Plan const &plan,
                                                const int region_to_draw_tree_on) {
    int V = map_params.V;
    // optional for checking 
    int num_region_vertices = 0;

    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++) {
        ignore[i] = plan.region_ids[i] != region_to_draw_tree_on;
        // count vertices we're not ignoring 
        num_region_vertices += !ignore[i];
    }

    // get upper and lower bounds on region pops
    auto min_max_pair = splitting_schedule.all_regions_min_and_max_possible_cut_sizes
                            [plan.region_sizes[region_to_draw_tree_on]];

    // clear the tree
    clear_tree(ust);

    // Now sample a uniform spanning tree drawn on that region
    bool const valid_tree = draw_ust(root, 
        min_max_pair.first * map_params.lower, // lower, 
        min_max_pair.second * map_params.upper, // upper
        rng_state
        );

    if constexpr(perf_config::object_integrity_checking){
        if (valid_tree){
            check_tree_integrity(
                ust,
                "Just called `sample_sub_ust` in attempt_to_draw_tree_on_region\n",
                root,
                num_region_vertices,
                true
            );
        }
    }

    return std::make_pair(valid_tree, num_region_vertices);
}

std::pair<bool, int> USTSampler::attempt_to_draw_tree_on_merged_region(RNGState &rng_state, Plan const &plan,
                                                       const int region1_to_draw_tree_on,
                                                       const int region2_to_draw_tree_on) {
    int V = map_params.V;
    // optional for checking 
    int num_merged_region_vertices = 0;

    // Mark it as ignore if its not in either of the two regions
    for (int i = 0; i < V; i++) {
        ignore[i] = plan.region_ids[i] != region1_to_draw_tree_on &&
                    plan.region_ids[i] != region2_to_draw_tree_on;
        num_merged_region_vertices += !ignore[i];
    }

    int merged_region_size =
        plan.region_sizes[region1_to_draw_tree_on] + plan.region_sizes[region2_to_draw_tree_on];
    // get upper and lower bounds on region pops
    auto min_max_pair =
        splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];

    // clear the tree
    clear_tree(ust);

    // Now sample a uniform spanning tree drawn on that region
    bool const valid_tree = draw_ust(root, 
        min_max_pair.first * map_params.lower, // lower, 
        min_max_pair.second * map_params.upper, // upper
        rng_state
        );

    if constexpr(perf_config::object_integrity_checking){
        if (valid_tree){
            check_tree_integrity(
                ust,
                "Just called `sample_sub_ust` in attempt_to_draw_tree_on_merged_region\n",
                root,
                num_merged_region_vertices,
                true
            );
        }
    }

    return std::make_pair(valid_tree, num_merged_region_vertices);
}

std::pair<bool, EdgeCut> USTSampler::try_to_sample_splittable_tree(
    Plan const &plan, int const split_region1, int const split_region2,
    ScoringFunction const &scoring_function, RNGState &rng_state, TreeSplitter &tree_splitter,
    int const region_populations, int const region_size, bool const save_selection_prob) {
    // We assume a tree has already been successfully drawn so
    // try to find a valid cut
    auto cut_size_bounds =
        splitting_schedule.all_regions_min_and_max_possible_cut_sizes[region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    // REprintf("Remainder Size: %d Cut sizes:", region_size);
    // for(auto const &v: splitting_schedule.all_regions_smaller_cut_sizes_to_try[region_size]){
    //     REprintf("%d, ", v);
    // }
    // REprintf("\n");

    std::pair<bool, EdgeCut> edge_search_result = tree_splitter.attempt_to_find_edge_to_cut(
        map_params, scoring_function, rng_state, plan, split_region1, split_region2, ust, root,
        stack, pops_below_vertex, ignore, region_populations, region_size,
        min_possible_cut_size, max_possible_cut_size,
        splitting_schedule.all_regions_smaller_cut_sizes_to_try[region_size],
        save_selection_prob);

    bool search_successful = std::get<0>(edge_search_result);
    // return false if unsuccessful
    if (!search_successful)
        return std::make_pair(false, EdgeCut());

    // If successful extract the edge cut info
    EdgeCut cut_edge = std::get<1>(edge_search_result);
    // Now erase the cut edge in the tree
    erase_tree_edge(ust, cut_edge);

    return edge_search_result;
}

std::pair<bool, EdgeCut> USTSampler::attempt_to_find_valid_tree_split(
    RNGState &rng_state, ScoringFunction const &scoring_function, TreeSplitter &tree_splitter,
    Plan const &plan, int const region_to_split, int const new_region_id,
    bool const save_selection_prob) {
    // Try to draw a tree
    bool tree_drawn = attempt_to_draw_tree_on_region(rng_state, plan, region_to_split).first;
    // return false if unsuccessful
    if (!tree_drawn)
        return std::make_pair(false, EdgeCut());

    // Else try to find a valid cut
    int region_to_split_size = plan.region_sizes[region_to_split];
    int region_to_split_population = plan.region_pops[region_to_split];

    return try_to_sample_splittable_tree(plan, region_to_split, new_region_id, scoring_function,
                                         rng_state, tree_splitter, region_to_split_population,
                                         region_to_split_size, save_selection_prob);
}

std::pair<bool, EdgeCut> USTSampler::attempt_to_find_valid_tree_mergesplit(
    RNGState &rng_state, ScoringFunction const &scoring_function, TreeSplitter &tree_splitter,
    Plan const &plan, int const merge_region1, int const merge_region2,
    bool const save_selection_prob) {
    // Try to draw a tree
    bool tree_drawn =
        attempt_to_draw_tree_on_merged_region(rng_state, plan, merge_region1, merge_region2).first;
    // return false if unsuccessful
    if (!tree_drawn)
        return std::make_pair(false, EdgeCut());

    // Else try to find a valid cut
    int region_to_split_size =
        plan.region_sizes[merge_region1] + plan.region_sizes[merge_region2];
    int region_to_split_population =
        plan.region_pops[merge_region1] + plan.region_pops[merge_region2];

    return try_to_sample_splittable_tree(plan, merge_region1, merge_region2, scoring_function,
                                         rng_state, tree_splitter, region_to_split_population,
                                         region_to_split_size, save_selection_prob);
}


void USTSampler::check_tree_integrity(
      Tree const &a_ust,
      std::string_view where,
      int root,
      int expected_tree_vertices,
      bool check_vertex_count
    ){
    // check no garbage vertices in the tree 
    check_directed_tree_edges_are_graph_edges(
        a_ust, map_params.graph_edge_index,
        where
    );
    // now check the tree returned is actually a directed tree
    check_is_directed_tree(
        a_ust,
        where,
        root,
        expected_tree_vertices,
        check_vertex_count,
        visited,
        stack
    );
}