#include "labeling.h"

double counter_helper(std::vector<bool> &A_in, int n_in, int add, const Graph &g,
                      std::map<std::vector<bool>, double> &memos);

double log_labelings_exact(const Graph &g) {
    int n = g.size();
    std::map<std::vector<bool>, double> memos;
    std::vector<bool> A_in(n, false);

    // start recursion from every vertex
    std::vector<double> xchild(n);
    double max_x = 0.0;
    for (int i = 0; i < n; i++) {
        xchild[i] = counter_helper(A_in, 0, i, g, memos);
        if (xchild[i] > max_x) {
            max_x = xchild[i];
        }
    }

    double accuml = 0.0;
    for (int i = 0; i < n; i++) {
        accuml += std::exp(xchild[i] - max_x);
    }

    return std::log(accuml);
}

// recursion helper
// Add `add` to the labelled set `A_in`, generate the new frontier, and recurse across it
// add to `memos` as we go to save computation
double counter_helper(std::vector<bool> &A_in, int n_in, int add, const Graph &g,
                      std::map<std::vector<bool>, double> &memos) {
    int n = A_in.size();
    if (n_in >= n - 2) return 0.0;
    A_in[add] = true;
    n_in++;
    auto search = memos.find(A_in);
    if (search != memos.end()) {
        A_in[add] = false;
        return search->second;
    } else {
        std::vector<double> xchild(n);
        std::vector<bool> skip(n, true);
        double max_x = 0.0;
        for (int i = 0; i < n; i++) {
            if (A_in[i]) continue;
            bool touches = false;
            std::vector<int> nbors = g[i];
            for (int nbor : nbors) {
                if (A_in[nbor]) {
                    touches = true;
                    break;
                }
            }
            if (!touches) continue;
            skip[i] = false;

            xchild[i] = counter_helper(A_in, n_in, i, g, memos);
            if (xchild[i] > max_x)
                max_x = xchild[i];
        }

        double accuml = 0.0;
        for (int i = 0; i < n; i++) {
            if (!skip[i]) {
                accuml += std::exp(xchild[i] - max_x);
            }
        }

        std::vector<bool> A_key(A_in);
        double out = std::log(accuml) + max_x;
        memos.emplace(A_key, out);
        A_in[add] = false;
        return out;
    }
}


double log_labelings_IS(const Graph &g, int n) {
    int V = g.size();
    vec weights(V);
    double tot_wgt = 0.0;
    for (int i = 0; i < V; i++) {
        weights[i] = std::sqrt(g[i].size());
        tot_wgt += weights[i];
    }

    vec lp(n);
    double min_lp = 1e6;
    double log_tot_wgt = std::log(tot_wgt);
    for (int i = 0; i < n; i++) {
        std::vector<bool> candidate(V, false);
        std::vector<bool> visited(V, false);

        double idx = tot_wgt * unif(generator);
        double accuml = 0;
        int vtx;
        for (vtx = 0; vtx < V - 1; vtx++) {
            accuml += weights[vtx];
            if (accuml >= idx) break;
        }
        lp[i] = std::log(weights[vtx]) - log_tot_wgt;

        visited[vtx] = true;
        std::vector<int> nbors = g[vtx];
        int n_nbors = nbors.size();
        double n_cands = 0;
        for (int k = 0; k < n_nbors; k++) {
            candidate.at(nbors[k]) = true;
            n_cands += weights.at(nbors[k]);
        }

        for (int j = 1; j < V; j++) {
            double idx = n_cands * unif(generator);
            double accuml = 0;
            int vtx;
            for (int k = 0; k < V; k++) {
                if (candidate[k]) {
                    vtx = k;
                    accuml += weights[vtx];
                    if (accuml >= idx) break;
                }
            }
            lp[i] += std::log(weights.at(vtx)) - std::log(n_cands);

            candidate[vtx] = false;
            visited[vtx] = true;
            std::vector<int> nbors = g[vtx];
            n_cands -= weights[vtx];
            int n_nbors = nbors.size();
            for (int k = 0; k < n_nbors; k++) {
                if (!visited.at(nbors[k]) && !candidate.at(nbors[k])) {
                    n_cands += weights[nbors[k]];
                    candidate[nbors[k]] = true;
                }
            }
        }

        if (lp[i] < min_lp) min_lp = lp[i];
    }

    return std::log(sum(exp(min_lp - lp)));
}
