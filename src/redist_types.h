#pragma once
#ifndef TYPES_H
#define TYPES_H

#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

typedef std::vector<std::vector<int>> Tree;
typedef std::vector<std::vector<int>> Graph;
typedef std::vector<std::vector<std::vector<int>>> Multigraph;

#endif
