#include <set>
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;

// forward declaration
// [[Rcpp::export]]
std::vector<std::set<int>> get_plan_graph(List l, int V,
                                          IntegerVector plan, int n_distr);

// helper
bool deg_sort(const std::pair<int, int> &a, const std::pair<int, int> &b) {
    return a.second > b.second;
}

// [[Rcpp::export]]
IntegerVector color_graph(List l, IntegerVector plan) {
    int n_distr = max(plan);
    int V = l.size();
    std::vector<std::set<int>> dist_gr = get_plan_graph(l, V, plan, n_distr);
    std::vector<std::pair<int, int>> degs(n_distr);

    // iterate from largest to smallest degree
    for (int i = 0; i < n_distr; i++) {
        degs[i] = std::make_pair(i, dist_gr[i].size());
    }
    std::sort(degs.begin(), degs.end(), deg_sort);

    std::vector<int> color(n_distr, 0);
    int colors = 4;
    color[degs[0].first] = 1; // first color
    for (int i = 1; i < n_distr; i++) {
        int curr = degs[i].first;
        std::vector<bool> seen(colors);
        std::set<int> nbors = dist_gr[curr];
        for (int nbor : nbors) {
            int nbor_color = color[nbor] - 1;
            if (nbor_color >= 0)
                seen[nbor_color] = true;
        }

        auto idx = std::find(seen.begin(), seen.end(), false);
        if (idx == seen.end()) colors++;
        color[curr] = idx - seen.begin() + 1;
    }


    IntegerVector out(V);
    for (int i = 0; i < V; i++) {
        out[i] = color[plan[i] - 1];
    }

    return out;
}

std::vector<std::set<int>> get_plan_graph(List l, int V,
                                          IntegerVector plan, int n_distr) {
    std::vector<std::set<int>> dist_gr;
    for (int i = 0; i < n_distr; i++) {
        dist_gr.push_back(std::set<int>());
    }

    for (int i = 0; i < V; i++) {
        IntegerVector nbors = (IntegerVector) l[i];
        int length = nbors.size();
        int distr = plan.at(i) - 1;
        for (int j = 0; j < length; j++) {
            int nbor_distr = plan.at(nbors[j]) - 1;
            if (distr == nbor_distr) continue;
            dist_gr.at(distr).insert(nbor_distr);
        }
    }

    return dist_gr;
}
