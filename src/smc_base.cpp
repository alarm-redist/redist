#include "smc_base.h"

/*
 * Partition `x` and its indices `idxs` between `right` and `left` by `pivot`
 */
// TESTED
void partition_vec(std::vector<double> &x, std::vector<int> &idxs, int left,
                   int right, int &pivot) {
    double pivot_value = x[pivot];
    std::swap(x[pivot], x[right]);
    std::swap(idxs[pivot], idxs[right]);
    pivot = left;
    for (int i = left; i < right; i++) {
        if (x[i] < pivot_value) {
            std::swap(x[pivot], x[i]);
            std::swap(idxs[pivot], idxs[i]);
            pivot++;
        }
    }
    std::swap(x[right], x[pivot]);
    std::swap(idxs[right], idxs[pivot]);
}

/*
 * Get the index of the k-th smallest element of x
 */
// TESTED
int global_rng_select_k(std::vector<double> x, int k) {
    if(k > x.size()){
        REprintf("k=%d bigger than number of edges=%d!\n", k, (int) x.size());
        throw Rcpp::exception("k bigger than number of edges!\n");
    }
    int right = x.size() - 1;
    int left = 0;
    std::vector<int> idxs(right + 1);
    for (int i = 0; i <= right; i++) idxs[i] = i;

    k--;
    while (true) {
        if (left == right)
            return idxs[left];
        int pivot = left + GLOBAL_RNG.r_int(right - left + 1);
        partition_vec(x, idxs, left, right, pivot);
        if (k == pivot) {
            return idxs[k];
        } else if (k < pivot) {
            right = pivot - 1;
        } else {
            left = pivot + 1;
        }
    }
}


/*
 * Get the index of the k-th smallest element of x
 */
// TESTED
int select_k(std::vector<double> x, int k, RNGState &rng_state) {
    if(k > x.size()){
        REprintf("k=%d bigger than number of edges=%d!\n", k, (int) x.size());
        throw Rcpp::exception("k bigger than number of edges!\n");
    }
    int right = x.size() - 1;
    int left = 0;
    std::vector<int> idxs(right + 1);
    for (int i = 0; i <= right; i++) idxs[i] = i;

    k--;
    while (true) {
        if (left == right)
            return idxs[left];
        int pivot = left + rng_state.r_int(right - left + 1);
        partition_vec(x, idxs, left, right, pivot);
        if (k == pivot) {
            return idxs[k];
        } else if (k < pivot) {
            right = pivot - 1;
        } else {
            left = pivot + 1;
        }
    }
}

List cli_config(bool clear, const char * fmt) {
    return List::create(_["clear"]=clear, _["show_after"]=0.25,
                        _["format"]=fmt);
}
