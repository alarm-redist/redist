#ifndef UTILS_H
#define UTILS_H

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#include "redist_types.h"
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]


#include "random.h"


/*
 * Creates a Rcpp Threadpool
 *
 *
 */
RcppThread::ThreadPool get_thread_pool(int const num_threads);

// Returns the number of threads being used. 
// When single threading it returns as one thread
int get_num_threads(const RcppThread::ThreadPool &pool);

/*
 * Make a progress bar configuration with format string `fmt`
 */
Rcpp::List cli_config(bool clear = false,
                const char *fmt = "{cli::pb_bar} {cli::pb_percent} | ETA:{cli::pb_eta}");


// Only relevant for debugging 
struct ActiveUserGuard {
    std::atomic<int> &counter;

    explicit ActiveUserGuard(std::atomic<int> &counter_) : counter(counter_) {
        int const old = counter.fetch_add(1, std::memory_order_acq_rel);

        if (old != 0) {
            throw std::runtime_error(
                "Two workers are using the same per-thread scratch slot."
            );
        }
    }

    ~ActiveUserGuard() {
        counter.fetch_sub(1, std::memory_order_acq_rel);
    }
};




// Prints a vector as `name` = c(<VECTOR CONTENT>)
template <typename T>
std::string vec_to_string(std::vector<T> const &x, std::string_view name) {
    std::ostringstream oss;

    oss << name << " = c(";
    for (std::size_t i = 0; i < x.size(); ++i) {
        oss << x[i];
        if (i + 1 < x.size()) {
            oss << ", ";
        }
    }
    oss << ")\n";

    return oss.str();
};



// Simple helper to make sure no tree size indexing 
void tree_size_check(
    MapParams const &map_params, 
    int const proposed_tree_size, 
    std::vector<int> const &tree_sizes,
    std::string_view const where
);


#endif
