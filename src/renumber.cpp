#include <Rcpp.h>
using namespace Rcpp;

#include "hungarian.h"

// row of output is m1 districts
// col of output is m2 districts
// [[Rcpp::export]]
NumericMatrix plan_joint(IntegerVector m1, IntegerVector m2, NumericVector pop) {
    int k = max(m1);
    int V = m1.size();
    NumericMatrix joint(k);
    NumericVector p1(k);
    NumericVector p2(k);

    for (int i = 0; i < V; i++) {
        joint(m1[i]-1, m2[i]-1) += pop[i];
        p1[m1[i]-1] += pop[i];
        p2[m2[i]-1] += pop[i];
    }

    return joint;
}

// [[Rcpp::export]]
IntegerMatrix renumber_matrix(IntegerMatrix plans, IntegerVector renumb) {
    int V = plans.nrow();
    int N = plans.ncol();
    int n_distr = max(plans(_, 0));
    IntegerMatrix out(V, N);

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < V; i++) {
            out(i, j) = renumb(n_distr*j + plans(i, j) - 1);
        }
    }

    return out;
}

// Hungarian Algorithm Solver
//
// Wrapper code from Justin Silverman, 2019
//
// Solves weighted bipartite matching problems (e.g., optimal matching of people to cars
// or optimal matching of students to colleges, etc...)
//
// @param costMatrix matrix giving cost of each possible pairing - can be rectangular
// @return List with cost and parings, pairings are given as an Nx2 matrix
// giving edges that are matched (1-indexed rather than 0-indexed as it will be returned to R)
// @details this is a copy/wrapper for the code developed by Cong Ma and made available
// as a github repository (mcximing/hungarian-algorithm-cpp). Code was
// changed to a header only file for use in other Rcpp packages.
//
// [[Rcpp::export]]
IntegerMatrix solve_hungarian(NumericMatrix costMatrix) {
    int nr = costMatrix.nrow();
    int nc = costMatrix.ncol();

    vector<double> c(nc);
    vector<vector<double>> cm(nr, c);
    for (int i=0; i < nr; i++){
        for (int j=0; j < nc; j++){
            c[j] = costMatrix(i,j);
        }
        cm[i] = c;
    }

    HungarianAlgorithm HungAlgo;
    vector<int> assignment;
    HungAlgo.Solve(cm, assignment);
    IntegerMatrix assign(nr, 2);
    for (int i=0; i < nr; i++){
        assign(i,0) = i+1;
        assign(i,1) = assignment[i]+1;
    }

    return assign;
}
