#ifndef UTILS_H
#define UTILS_H

#include <atomic>
#include <limits>
#include <vector>

#include "advanced_types.h"
#include <Rcpp.h>

/*
 * Make a progress bar configuration with format string `fmt`
 */
Rcpp::List cli_config(bool clear = false,
                const char *fmt = "{cli::pb_bar} {cli::pb_percent} | ETA:{cli::pb_eta}");

// Only relevant for debugging 



// gets current time if tracking granular times
inline Clock::time_point maybe_now() {
    if constexpr (perf_config::track_granular_times) {
        return Clock::now();
    } else {
        return {};
    }
}

// Adds the time elapsed to slot 
inline void add_elapsed(
    double &slot,
    Clock::time_point const start
) {
    if constexpr (perf_config::track_granular_times) {
        slot += std::chrono::duration<double, std::ratio<1>>(
            Clock::now() - start
        ).count();
    }else{
        throw Rcpp::exception("Time elapsed is being called when TRACK_GRANULAR_PERFORMANCE_TIMES = false");
    }
}


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
