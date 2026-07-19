#include "lct_graph_plan_type.h"

#include "wilson.h"

#include <algorithm>
#include <queue>

const EdgeSet LCTGraphPlan::empty_edge_set;

CWEdge::CWEdge(int a, int b, double w) : u(std::min(a, b)), v(std::max(a, b)), weight(w) {}

bool CWEdge::operator<(const CWEdge &other) const {
    if (u != other.u)
        return u < other.u;
    if (v != other.v)
        return v < other.v;
    return weight < other.weight;
}

bool CWEdge::operator==(const CWEdge &other) const { return u == other.u && v == other.v; }

bool LCTGraphPlan::init_lct_from_regions(MapParams const &map_params, USTSampler &ust_sampler,
                                         RNGState &rng_state,
                                         int const max_attempts_per_region) {
    int const V = map_params.V;
    lct = std::make_unique<LinkCutTree>(V);
    district_roots.assign(num_regions, -1);

    for (int d = 0; d < num_regions; ++d) {
        bool ok = false;
        USTDrawResult ust_draw_result;
        for (int attempt = 0; attempt < max_attempts_per_region && !ok; ++attempt) {
            ust_draw_result = ust_sampler.attempt_to_draw_tree_on_region(rng_state, *this, d);
            ok = ust_draw_result.successful;
        }
        if (!ok) {
            return false;
        }

        // ust_sampler.ust is rooted at ust_sampler.root after the draw.
        district_roots[d] = ust_draw_result.root;
        load_tree_into_lct(ust_sampler.ust, ust_draw_result.root, d);
    }

    rebuild_cross_edges(map_params.g);
    return true;
}

void LCTGraphPlan::load_tree_into_lct(Tree const &tree, int root, int /*district*/) {
    // BFS from root; for each child, evert + link to its parent.
    std::queue<int> q;
    q.push(root);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (int child : tree[v]) {
            lct->evert(child);
            lct->link(child, v);
            q.push(child);
        }
    }
}

void LCTGraphPlan::rebuild_cross_edges(Graph const &g) {
    cross_edges.clear();
    int const V = static_cast<int>(g.size());
    for (int u = 0; u < V; ++u) {
        int du = region_ids[u];
        for (int v : g[u]) {
            if (v <= u)
                continue;
            int dv = region_ids[v];
            if (du != dv) {
                DistrictPair key(std::min(du, dv), std::max(du, dv));
                cross_edges[key].insert(CWEdge(u, v));
            }
        }
    }
}

const EdgeSet &LCTGraphPlan::get_cross_edges(int d1, int d2) const {
    DistrictPair key(std::min(d1, d2), std::max(d1, d2));
    auto it = cross_edges.find(key);
    if (it != cross_edges.end())
        return it->second;
    return empty_edge_set;
}

bool LCTGraphPlan::districts_adjacent(int d1, int d2) const {
    DistrictPair key(std::min(d1, d2), std::max(d1, d2));
    return cross_edges.find(key) != cross_edges.end();
}

std::vector<DistrictPair> LCTGraphPlan::get_adjacent_district_pairs() const {
    std::vector<DistrictPair> pairs;
    pairs.reserve(cross_edges.size());
    for (auto const &kv : cross_edges) {
        pairs.push_back(kv.first);
    }
    return pairs;
}

void LCTGraphPlan::set_edge_weights(MapParams const &map_params,
                                    Rcpp::List const &edge_weights_list) {
    edge_weights.clear();
    if (edge_weights_list.size() == 0)
        return;

    int const V = map_params.V;
    for (int i = 0; i < edge_weights_list.size(); ++i) {
        Rcpp::List entry = edge_weights_list[i];
        if (!entry.containsElementNamed("edge"))
            Rcpp::stop("Edge weight entry %d missing 'edge' field", i + 1);
        Rcpp::IntegerVector edge = entry["edge"];
        if (edge.size() != 2)
            Rcpp::stop("Edge weight entry %d: 'edge' must be length 2", i + 1);
        if (!entry.containsElementNamed("weight"))
            Rcpp::stop("Edge weight entry %d missing 'weight' field", i + 1);
        double weight = Rcpp::as<double>(entry["weight"]);
        if (weight <= 0)
            Rcpp::stop("Edge weight entry %d: weight must be positive", i + 1);

        // R is 1-indexed.
        int u = edge[0] - 1;
        int v = edge[1] - 1;
        if (u < 0 || u >= V || v < 0 || v >= V)
            Rcpp::stop("Edge weight entry %d: vertices out of range [1, %d]", i + 1, V);

        auto const &neighbors = map_params.g[u];
        if (std::find(neighbors.begin(), neighbors.end(), v) == neighbors.end()) {
            Rcpp::stop("Edge weight entry %d: edge (%d, %d) not in adjacency graph", i + 1,
                       edge[0], edge[1]);
        }
        edge_weights[{std::min(u, v), std::max(u, v)}] = weight;
    }
}

double LCTGraphPlan::get_edge_weight(int u, int v) const {
    auto key = std::make_pair(std::min(u, v), std::max(u, v));
    auto it = edge_weights.find(key);
    if (it != edge_weights.end())
        return it->second;
    return default_edge_weight;
}

void LCTGraphPlan::print_state(int verbosity) const {
    if (verbosity < 3)
        return;
    Rcpp::Rcout << "[LCTGraphPlan] " << num_regions << " regions\n";
    Rcpp::Rcout << "[LCTGraphPlan] District roots: ";
    for (int i = 0; i < num_regions; ++i) {
        Rcpp::Rcout << district_roots[i];
        if (i < num_regions - 1)
            Rcpp::Rcout << ", ";
    }
    Rcpp::Rcout << "\n[LCTGraphPlan] Region populations: ";
    for (int i = 0; i < num_regions; ++i) {
        Rcpp::Rcout << region_pops[i];
        if (i < num_regions - 1)
            Rcpp::Rcout << ", ";
    }
    Rcpp::Rcout << "\n[LCTGraphPlan] Adjacent pairs: " << cross_edges.size() << "\n";
}

void LCTGraphPlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter, USTSampler &ust_sampler, EdgeCut const cut_edge,
    int const split_region1_id, int const split_region2_id, bool const add_region) {
    // Update region_ids via GraphPlan's behavior.
    GraphPlan::update_vertex_and_plan_specific_info_from_cut(
        tree_splitter, ust_sampler, cut_edge, split_region1_id, split_region2_id, add_region);
    // LCT consistency after a split is only needed when LCTGraphPlan is used
    // as an SMC move. Cyclewalk's MCMC inner loop manages LCT itself in
    // attempt_cycle_walk / attempt_forest_walk and never goes through this
    // hook. When the SMC integration lands, rebuild spanning trees on the
    // two affected regions via USTSampler and rebuild_cross_edges here.
    throw Rcpp::exception(
        "LCTGraphPlan::update_vertex_and_plan_specific_info_from_cut not implemented; "
        "LCTGraphPlan is currently only usable as an MCMC state, not as an SMC step.");
}
