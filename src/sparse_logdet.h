// sparse_logdet.h
#pragma once

#include <vector>

// Minimal stand-in for Eigen::Triplet: one instance per nonzero, same layout.
//
// Eigen is deliberately kept out of this header. <Eigen/Sparse> is expensive to
// parse and this header is reached by map_calc.h, which nearly every
// translation unit includes. Callers build entries with this type and Eigen
// stays confined to sparse_logdet.cpp.
struct SparseEntry {
    int row;
    int col;
    double value;

    SparseEntry(int const row, int const col, double const value)
        : row(row), col(col), value(value) {}
};

// Log determinant of a sparse symmetric positive definite matrix, given only
// its upper triangle (ie every entry (i, j, value) has i <= j). Repeated
// coordinates are summed, so degrees may be appended as separate diagonal
// entries.
double compute_log_det_from_entries(std::vector<SparseEntry> const &entries,
                                    int const num_rows);
