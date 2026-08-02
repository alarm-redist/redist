#ifndef RANDOM_H
#define RANDOM_H

#include <RcppArmadillo.h>
#include <cstdint>
#include <random>
#include <vector>
#include "redist_constants.h"

// random number generator class
class RNGState {
  private:
    uint64_t state_sr;
    uint32_t state_xo[4];

    // rotate-left helper for xoshiro128++
    static inline uint32_t rotl32(uint32_t const x, int const k) {
        return (x << k) | (x >> (32 - k));
    }

    void long_jump(); // Does Long Jump https://prng.di.unimi.it/

  public:
    RNGState(int seed, int num_jumps = 0); // seed
    RNGState(void);                        // seed

    uint64_t next_sr();                    // moves to next state

    // xoshiro128++ generator
    inline uint32_t generator() {
        uint32_t const result = rotl32(state_xo[0] + state_xo[3], 7) + state_xo[0];

        uint32_t const t = state_xo[1] << 9;

        state_xo[2] ^= state_xo[0];
        state_xo[3] ^= state_xo[1];
        state_xo[1] ^= state_xo[2];
        state_xo[0] ^= state_xo[3];

        state_xo[2] ^= t;

        state_xo[3] = rotl32(state_xo[3], 11);

        return result;
    }

    void seed_rng(int seed, int num_jumps = 0); // seed it
    void do_long_jumps(int num_jumps);          // Does num_jumps number of long jumps

    // Generate a uniform random integer in [0, max). Slightly biased.
    inline int r_int(uint32_t const max) {
        if constexpr (perf_config::supposedly_safe_input_checks) {
            if (max == 0) {
                throw std::runtime_error("RNGState::r_int called with max=0.");
            }
        }

        uint32_t const x = generator();
        uint64_t const m = static_cast<uint64_t>(x) * static_cast<uint64_t>(max);
        return static_cast<int>(m >> 32);
    }

    // Generate a uniform random double in [0, 1). Slightly biased.
    inline double r_unif() {
        return 0x1.0p-32 * generator();
    }

    // Generate a random integer in [0, cum_wgts.size()) according
    // to cumulative normalized weights.
    int r_int_wgt(const arma::vec &cum_wgts);

    // Generate random integer with probability proporitional
    // to weights
    int r_int_unnormalized_wgt(const arma::vec &unnormalized_wgts);

    // Delete copy operator
    RNGState(const RNGState &) = delete;
    RNGState &operator=(const RNGState &) = delete;

    // Allow move semantics
    RNGState(RNGState &&) = default;
    RNGState &operator=(RNGState &&) = default;
};



/*
 * Set global RNGState seed
 */
void global_seed_rng(int seed, int num_jumps = 1);

/*
 * Generate a random integer within a stratum with some probability p
 */
int r_int_mixstrat(int max, int stratum, double p, arma::vec cum_wgts);

/*
 * Get the index of the k-th smallest element of x using global state
 * NOT THREAD SAFE
 */
// TESTED
int global_rng_select_k(std::vector<double> x, int k);


// legacy code. Global RNG state
// NOT THREAD SAFE
extern std::random_device GLOBAL_RD;
extern RNGState GLOBAL_RNG;

#endif
