#include <RcppArmadillo.h>
#include "redist_types.h"

// anonymous namespace used here so we no longer need to link to Tree_ops
namespace {
    Tree init_tree(int V) {
    Tree tree;
    for (int i = 0; i < V; i++) {
        tree.push_back(std::vector<int>());
    }
    return tree;
}
}

// [[Rcpp::export]]
Rcpp::List reduce_adj(Rcpp::List adj_list, Rcpp::IntegerVector prec_map, int n_keep) {
    Rcpp::List adj(n_keep);
    Rcpp::IntegerVector nbors, sub;

    int V = adj_list.size();
    for (int i = 0; i < V; i++) {
        int new_i = prec_map[i];
        if (new_i == -1)
            continue;

        sub = Rcpp::IntegerVector::create();

        nbors = adj_list[i];
        int n_nbors = nbors.size();
        for (int j = 0; j < n_nbors; j++) {
            int new_j = prec_map[nbors[j]];
            if (new_j >= 0)
                sub.push_back(new_j);
        }

        adj[new_i] = sub.sort();
    }

    return adj;
}

// Create the quotient graph of `graph` by the equiv relation induced by `idxs`,
// with multiple edges and self-loops removed.
// idxs should run continuously from 0 to n_groups-1
// [[Rcpp::export]]
Graph collapse_adj(Rcpp::List graph, const arma::uvec &idxs) {
    int V = graph.size();
    int V_new = max(idxs) + 1;
    Graph collapsed = init_tree(V_new);

    for (int i = 0; i < V; i++) {
        int from = idxs(i);
        std::vector<int> *nbors = &collapsed[from];
        int length = ((Rcpp::IntegerVector)graph[i]).size();
        for (int j = 0; j < length; j++) {
            int to = idxs(((Rcpp::IntegerVector)graph[i])[j]);

            if (from != to && std::find(nbors->begin(), nbors->end(), to) == nbors->end()) {
                nbors->push_back(to);
            }
        }
    }

    return collapsed;
}
