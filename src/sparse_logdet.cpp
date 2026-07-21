
#include "sparse_logdet.h"


// alias for SparseMatrix
using SparseMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

// This assumes the matrix is stored as upper triangular one via
// triplets (ie every triplet (i, j, value) we have i <= j )
double compute_log_det_from_triplets(std::vector<Eigen::Triplet<double, int>> const &trips,
                                     int const num_rows) {
    // now make the sparse matrix
    SparseMat sparse_adj_mat(num_rows, num_rows);
    sparse_adj_mat.setFromTriplets(trips.begin(), trips.end());
    sparse_adj_mat.makeCompressed();
    // now factor it

    // 1) Create a sparse Cholesky (LLᵀ) factorization object.
    Eigen::SimplicialLLT<SparseMat> chol;

    // 2) Analyze + factorize the matrix.
    //    - analyzePattern() inspects the sparsity structure (symbolic factorization / ordering)
    //    - factorize() does the numeric factorization
    //    compute() does both in one call (analyzePattern + factorize).
    // Note we assume it is an upper triangular matrix
    chol.compute(sparse_adj_mat.selfadjointView<Eigen::Upper>());

    // 3) Check whether factorization succeeded.
    //    Failure usually means: matrix is not SPD (singular/indefinite) or numerical breakdown.
    if (chol.info() != Eigen::Success) {
        std::ostringstream oss;
        oss << "Sparse Matrix in compute_log_det_from_triplets ";
        oss << "was not semi-positive definite" << "\n";
        throw std::runtime_error(oss.str());

        return -INFINITY; // or throw, depending on how you want to handle failure
    }

    // Avoid materializing L: iterate the diagonal directly in the stored sparse factor
    const auto Lview = chol.matrixL();           // TriangularView (in your build)
    const auto &Lmat = Lview.nestedExpression(); // underlying SparseMatrix

    double sumlog = 0.0;
    for (int j = 0; j < num_rows; ++j) {
        bool found = false;

        // Column j scan
        for (SparseMat::InnerIterator it(Lmat, j); it; ++it) {
            if (it.row() == j) { // diagonal entry
                const double d = it.value();
                if (!(d > 0.0)) {
                    std::ostringstream oss;
                    oss << "In compute_log_det_from_triplets d was not positive";
                    oss << "\n";
                    throw std::runtime_error(oss.str());
                    return -INFINITY;
                }
                sumlog += std::log(d);
                found = true;
                break;
            }
            // For lower-triangular L, rows are >= j; if we passed j, diagonal is missing
            if (it.row() > j)
                break;
        }

        if (!found) {
            std::ostringstream oss;
            oss << "In compute_log_det_from_triplets found was false";
            oss << "\n";
            throw std::runtime_error(oss.str());
            return -INFINITY; // should not happen for successful LLT
        }
    }

    return 2.0 * sumlog;
}