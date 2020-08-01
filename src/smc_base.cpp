#include "smc_base.h"

std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<double> unif(0.0, 1.0);

/*
 * Generate a uniform random integer in [0, max).
 */
int rint(int max) {
    return std::floor(max * unif(generator));
}


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
int select_k(std::vector<double> x, int k) {
    int right = x.size() - 1;
    int left = 0;
    std::vector<int> idxs(right + 1);
    for (int i = 0; i <= right; i++) idxs[i] = i;

    k--;
    while (true) {
        if (left == right)
            return idxs[left];
        int pivot = left + rint(right - left + 1);
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
