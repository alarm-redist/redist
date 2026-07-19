#pragma once
#ifndef GREDIST_TYPES_H
#define GREDIST_TYPES_H

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#define PRINT_LN Rcout << __func__ << "(), " << __FILE__ << ":" << __LINE__ << "\n";

#include <vector>
#include <unordered_map>


// old types
typedef std::vector<std::vector<int>> Tree;
typedef std::vector<std::vector<int>> Graph;
typedef std::vector<std::vector<std::array<int, 3>>> Multigraph;

typedef std::vector<std::unordered_map<int, int>> RegionMultigraphCount;



#endif