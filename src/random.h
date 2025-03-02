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

        

    public:
        RNGState(int seed); // seed 
        RNGState(void); // seed 
        uint64_t next_sr(); // moves to next state
        uint32_t generator();
        void seed_rng(int seed); // seed it 

        int r_int(uint32_t max); // Generate a uniform random integer in [0, max). Slightly biased.
        double r_unif(); // Generate a uniform random double in [0, 1). Slightly biased.
        int r_int_wgt(int max, vec cum_wgts); // Generate a random integer in [0, max) according to cumulative normalized weights.
        int r_int_unnormalized_wgt(const vec &unnormalized_wgts); // Generate random integer with probability proporitional to weights



        RNGState(const RNGState&) = delete;
        RNGState& operator=(const RNGState&) = delete;
};

/*
 * Set global RNGState seed
 */
void global_seed_rng(int seed);



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
