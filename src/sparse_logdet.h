// sparse_logdet.h
#pragma once

#include <Eigen/Sparse>
#include <vector>

double compute_log_det_from_triplets(
    std::vector<Eigen::Triplet<double, int>> const &trips,
    int const num_rows
);