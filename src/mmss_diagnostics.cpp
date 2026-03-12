#include "mmss.h"

#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>

namespace {

using Count = long double;
using StateMap = std::vector<std::unordered_map<int, Count>>;

inline bool pop_ok(double pop_sum, double lower, double upper) {
    return lower <= pop_sum && pop_sum <= upper;
}

int find_root(const Tree &tree) {
    int V = tree.size();
    std::vector<int> indeg(V, 0);
    for (int u = 0; u < V; u++) {
        for (int v : tree[u]) indeg[v]++;
    }

    for (int v = 0; v < V; v++) {
        if (indeg[v] == 0) return v;
    }
    return -1;
}

std::vector<std::vector<int>> tree_to_undirected(const Tree &tree) {
    int V = tree.size();
    std::vector<std::vector<int>> undirected(V);
    for (int u = 0; u < V; u++) {
        for (int v : tree[u]) {
            undirected[u].push_back(v);
            undirected[v].push_back(u);
        }
    }
    return undirected;
}

bool edge_is_cut(int u, int v, const std::vector<int> &parent,
                 const std::vector<bool> &cut_child) {
    return (parent[u] == v && cut_child[u]) || (parent[v] == u && cut_child[v]);
}

bool valid_cut_set(const std::vector<std::vector<int>> &undirected, const arma::uvec &pop,
                   const std::vector<int> &parent, const std::vector<bool> &cut_child,
                   int l_split, double lower, double upper) {
    int V = undirected.size();
    std::vector<bool> seen(V, false);
    std::vector<int> stack;
    stack.reserve(V);
    int n_comp = 0;

    for (int start = 0; start < V; start++) {
        if (seen[start]) continue;
        n_comp++;

        double comp_pop = 0.0;
        stack.clear();
        stack.push_back(start);
        seen[start] = true;

        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            comp_pop += pop[u];

            for (int v : undirected[u]) {
                if (seen[v] || edge_is_cut(u, v, parent, cut_child)) continue;
                seen[v] = true;
                stack.push_back(v);
            }
        }

        if (!pop_ok(comp_pop, lower, upper)) return false;
    }

    return n_comp == l_split;
}

Count count_partitions_enum(const Tree &tree, const arma::uvec &pop, int root,
                            const std::vector<int> &parent, int l_split,
                            double lower, double upper) {
    int V = tree.size();
    int n_cuts = l_split - 1;
    if (n_cuts == 0) {
        double total_pop = arma::accu(pop);
        return pop_ok(total_pop, lower, upper) ? 1.0L : 0.0L;
    }

    std::vector<int> cut_candidates;
    cut_candidates.reserve(V - 1);
    for (int v = 0; v < V; v++) {
        if (v != root) cut_candidates.push_back(v);
    }

    std::vector<std::vector<int>> undirected = tree_to_undirected(tree);
    std::vector<bool> cut_child(V, false);
    Count total = 0.0L;

    std::function<void(int, int)> recurse = [&](int start, int left) {
        if (left == 0) {
            if (valid_cut_set(undirected, pop, parent, cut_child, l_split, lower, upper)) {
                total += 1.0L;
            }
            return;
        }

        int max_start = (int) cut_candidates.size() - left;
        for (int i = start; i <= max_start; i++) {
            int v = cut_candidates[i];
            cut_child[v] = true;
            recurse(i + 1, left - 1);
            cut_child[v] = false;
        }
    };

    recurse(0, n_cuts);
    return total;
}

StateMap count_partition_states(const Tree &tree, const arma::uvec &pop, int node,
                                int max_closed, double lower, double upper) {
    StateMap states(max_closed + 1);
    states[0][pop[node]] = 1.0L;

    for (int child : tree[node]) {
        StateMap child_states = count_partition_states(tree, pop, child, max_closed, lower, upper);
        StateMap next(max_closed + 1);

        for (int k0 = 0; k0 <= max_closed; k0++) {
            for (const auto &state0 : states[k0]) {
                int open0 = state0.first;
                Count count0 = state0.second;

                for (int k1 = 0; k0 + k1 <= max_closed; k1++) {
                    for (const auto &state1 : child_states[k1]) {
                        int open1 = state1.first;
                        Count count1 = state1.second;
                        Count ways = count0 * count1;

                        next[k0 + k1][open0 + open1] += ways;
                        if (k0 + k1 + 1 <= max_closed && pop_ok(open1, lower, upper)) {
                            next[k0 + k1 + 1][open0] += ways;
                        }
                    }
                }
            }
        }

        states = std::move(next);
    }

    return states;
}

Count count_partitions_dp(const Tree &tree, const arma::uvec &pop, int root,
                          int l_split, double lower, double upper) {
    int n_cuts = l_split - 1;
    StateMap states = count_partition_states(tree, pop, root, n_cuts, lower, upper);
    Count total = 0.0L;
    for (const auto &state : states[n_cuts]) {
        if (pop_ok(state.first, lower, upper)) total += state.second;
    }
    return total;
}

int resolve_method(const std::string &method, int l_split, int region_size) {
    if (method == "enumeration") return 1;
    if (method == "tree_dp") return 2;
    if (method == "auto") return (l_split <= 3 && region_size <= 30) ? 1 : 2;
    Rcpp::stop("`method` must be one of 'auto', 'enumeration', or 'tree_dp'.");
}

Count count_partitions(const Tree &tree, const arma::uvec &pop, int l_split,
                       double lower, double upper, int method_code) {
    int root = find_root(tree);
    if (root < 0) Rcpp::stop("Invalid tree: could not determine root.");

    std::vector<int> pop_below(tree.size(), 0);
    std::vector<int> parent(tree.size(), -1);
    parent[root] = -1;
    Tree tree_copy = tree;
    tree_pop(tree_copy, root, pop, pop_below, parent);

    if (method_code == 1) {
        return count_partitions_enum(tree, pop, root, parent, l_split, lower, upper);
    }
    return count_partitions_dp(tree, pop, root, l_split, lower, upper);
}

bool reference_partition_in_tree(const Tree &tree, const arma::uvec &pop,
                                 const std::vector<int> &ref_label,
                                 int l_split, double lower, double upper) {
    std::vector<double> group_pop(l_split, 0.0);
    std::vector<bool> seen_group(l_split, false);

    for (int v = 0; v < (int) ref_label.size(); v++) {
        int label = ref_label[v];
        if (label < 0 || label >= l_split) return false;
        group_pop[label] += pop[v];
        seen_group[label] = true;
    }

    for (int j = 0; j < l_split; j++) {
        if (!seen_group[j] || !pop_ok(group_pop[j], lower, upper)) return false;
    }

    int crossing_edges = 0;
    for (int u = 0; u < (int) tree.size(); u++) {
        for (int v : tree[u]) {
            if (ref_label[u] != ref_label[v]) crossing_edges++;
        }
    }

    return crossing_edges == l_split - 1;
}

void extract_region_tree(const Tree &full_tree, const std::vector<int> &region_vertices,
                         arma::uvec &region_pop, std::vector<int> &ref_label,
                         Tree &local_tree) {
    int V = full_tree.size();
    std::vector<int> local_idx(V, -1);
    int region_size = region_vertices.size();
    local_tree = init_tree(region_size);

    for (int i = 0; i < region_size; i++) {
        local_idx[region_vertices[i]] = i;
    }

    for (int i = 0; i < region_size; i++) {
        int old_u = region_vertices[i];
        for (int old_v : full_tree[old_u]) {
            int new_v = local_idx[old_v];
            if (new_v >= 0) local_tree[i].push_back(new_v);
        }
    }
}

} // namespace

// [[Rcpp::export]]
Rcpp::List count_single_tree_partitions_impl(Tree tree, const arma::uvec &pop,
                                             int l_split, double lower, double upper,
                                             std::string method = "auto") {
    if ((int) pop.n_elem != (int) tree.size()) {
        Rcpp::stop("Population vector must match the tree size.");
    }
    if (l_split < 1 || l_split > (int) tree.size()) {
        Rcpp::stop("`l_split` must be between 1 and the number of tree vertices.");
    }

    int method_code = resolve_method(method, l_split, (int) tree.size());
    Count count = count_partitions(tree, pop, l_split, lower, upper, method_code);

    return Rcpp::List::create(
        Rcpp::_["count"] = (double) count,
        Rcpp::_["method"] = method_code == 1 ? "enumeration" : "tree_dp"
    );
}

// [[Rcpp::export]]
Rcpp::List diag_single_tree_partitions_impl(List l, const arma::uvec &pop,
                                            const Rcpp::IntegerMatrix &plans,
                                            int n_distr, double pop_tol,
                                            int l_split, int n_trees_per_region) {
    Graph g = list_to_graph(l);
    int V = g.size();
    int n_plans = plans.ncol();

    if ((int) pop.n_elem != V || plans.nrow() != V) {
        Rcpp::stop("Plan matrix and population vector must match the map size.");
    }
    if (l_split < 1 || l_split > n_distr) {
        Rcpp::stop("`l` must be between 1 and the number of districts in the map.");
    }
    if (n_plans < 1 || n_trees_per_region < 1) {
        Rcpp::stop("`n_plans` and `n_trees_per_region` must be positive.");
    }

    arma::uvec counties(V, arma::fill::ones);
    Multigraph cg = county_graph(g, counties);
    Tree full_tree = init_tree(V);
    std::vector<bool> visited(V, false);

    std::vector<double> m_values;
    std::vector<double> z_values;
    m_values.reserve(n_plans * n_trees_per_region);
    z_values.reserve(n_plans * n_trees_per_region);
    std::map<long long, int> dist_m;

    for (int plan_idx = 0; plan_idx < n_plans; plan_idx++) {
        arma::uvec plan(V);
        for (int v = 0; v < V; v++) plan[v] = plans(v, plan_idx);

        Graph dist_g = district_graph(g, plan, n_distr);
        double log_prob = 0.0;
        std::vector<int> selected = select_l_districts(n_distr, dist_g, l_split, log_prob);
        if ((int) selected.size() != l_split || !std::isfinite(log_prob)) {
            Rcpp::stop("Failed to select a connected set of districts.");
        }

        std::vector<int> district_to_ref(n_distr + 1, -1);
        for (int j = 0; j < l_split; j++) district_to_ref[selected[j]] = j;

        std::vector<bool> ignore(V, true);
        std::vector<int> region_vertices;
        region_vertices.reserve(V);
        double region_pop_total = 0.0;
        for (int v = 0; v < V; v++) {
            int ref_idx = district_to_ref[plan[v]];
            if (ref_idx >= 0) {
                ignore[v] = false;
                region_vertices.push_back(v);
                region_pop_total += pop[v];
            }
        }

        double target = region_pop_total / l_split;
        double lower = (1.0 - pop_tol) * target;
        double upper = (1.0 + pop_tol) * target;

        arma::uvec region_pop(region_vertices.size());
        std::vector<int> ref_label(region_vertices.size());
        for (int i = 0; i < (int) region_vertices.size(); i++) {
            int v = region_vertices[i];
            region_pop[i] = pop[v];
            ref_label[i] = district_to_ref[plan[v]];
        }

        int method_code = resolve_method("auto", l_split, (int) region_vertices.size());

        for (int tree_idx = 0; tree_idx < n_trees_per_region; tree_idx++) {
            int attempts = 0;
            int result = 1;
            int root = -1;
            while (attempts < 100 && result != 0) {
                clear_tree(full_tree);
                result = sample_sub_ust(g, full_tree, V, root, visited, ignore, pop,
                                        0.0, region_pop_total, counties, cg);
                attempts++;
            }
            if (result != 0) {
                Rcpp::stop("Failed to draw a spanning tree on a sampled merged region.");
            }

            Tree region_tree;
            extract_region_tree(full_tree, region_vertices, region_pop, ref_label, region_tree);

            Count m_count = count_partitions(region_tree, region_pop, l_split, lower, upper, method_code);
            double m_value = (double) m_count;
            bool ref_in_tree = reference_partition_in_tree(region_tree, region_pop, ref_label,
                                                           l_split, lower, upper);
            double z_value = (ref_in_tree && m_value > 0.0) ? 1.0 / m_value : 0.0;

            m_values.push_back(m_value);
            z_values.push_back(z_value);
            dist_m[(long long) std::llround(m_value)]++;
        }
    }

    int n_total = m_values.size();
    double hit_rate = 0.0;
    double mean_m = 0.0;
    double max_m = 0.0;
    for (double value : m_values) {
        if (value > 0.0) hit_rate += 1.0;
        mean_m += value;
        if (value > max_m) max_m = value;
    }
    hit_rate /= n_total;
    mean_m /= n_total;

    double mean_z = 0.0;
    for (double z : z_values) mean_z += z;
    mean_z /= n_total;

    double var_z = 0.0;
    for (double z : z_values) {
        double diff = z - mean_z;
        var_z += diff * diff;
    }
    if (n_total > 1) var_z /= (n_total - 1);

    Rcpp::IntegerVector dist_vals(dist_m.size());
    Rcpp::CharacterVector dist_names(dist_m.size());
    int idx = 0;
    for (const auto &entry : dist_m) {
        dist_vals[idx] = entry.second;
        dist_names[idx] = std::to_string(entry.first);
        idx++;
    }
    dist_vals.attr("names") = dist_names;

    return Rcpp::List::create(
        Rcpp::_["hit_rate"] = hit_rate,
        Rcpp::_["mean_m"] = mean_m,
        Rcpp::_["max_m"] = max_m,
        Rcpp::_["dist_m"] = dist_vals,
        Rcpp::_["cv2_Z"] = mean_z > 0.0 ? var_z / (mean_z * mean_z) : R_PosInf
    );
}
