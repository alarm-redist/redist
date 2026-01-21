#ifndef RANDOM_H
#define RANDOM_H

#include <RcppArmadillo.h>
#include <vector>
#include <cstdint>
#include <random>

using namespace arma;

// random number generator class
class RNGState{
    private:
        uint64_t state_sr;
        uint32_t state_xo[4];
        void long_jump(); // Does Long Jump https://prng.di.unimi.it/
        

    public:
        RNGState(int seed, int num_jumps = 0); // seed 
        RNGState(void); // seed 
        uint64_t next_sr(); // moves to next state
        uint32_t generator();
        void seed_rng(int seed, int num_jumps = 0); // seed it 
        void do_long_jumps(int num_jumps); // Does num_jumps number of long jumps
        

        int r_int(uint32_t max); // Generate a uniform random integer in [0, max). Slightly biased.
        double r_unif(); // Generate a uniform random double in [0, 1). Slightly biased.
        int r_int_wgt(vec cum_wgts); // Generate a random integer in [0, cum_wgts.size()) according to cumulative normalized weights.
        int r_int_unnormalized_wgt(const vec &unnormalized_wgts); // Generate random integer with probability proporitional to weights

        // Delete copy operator 
        RNGState(const RNGState&) = delete;
        RNGState& operator=(const RNGState&) = delete;

        // Allow move semantics
        RNGState(RNGState&&) = default;
        RNGState& operator=(RNGState&&) = default;
};

/*
 * Set global RNGState seed
 */
void global_seed_rng(int seed, int num_jumps = 1);



/*
 * Generate a random integer within a stratum with some probability p
 */
int r_int_mixstrat(int max, int stratum, double p, vec cum_wgts);

/*
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
// [[Rcpp::export]]
arma::ivec resample_lowvar(arma::vec wgts);



// legacy code. Global RNG state
// NOT THREAD SAFE
static std::random_device GLOBAL_RD;
static RNGState GLOBAL_RNG(GLOBAL_RD());

#endif
