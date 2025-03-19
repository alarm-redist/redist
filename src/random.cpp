#include "random.h"



/* This is a fixed-increment version of Java 8's SplittableRandom generator
 See http://dx.doi.org/10.1145/2714064.2660195 and
 http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html
 Written in 2015 by Sebastiano Vigna (vigna@acm.org)
 [Public Domain]
 */

/* This is xoshiro128++ 1.0.
 Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 [Public domain]
 */
static inline uint32_t rotl(const uint32_t x, int k) {
    return (x << k) | (x >> (32 - k));
}



// Set the state for RNG state object
void RNGState::seed_rng(int seed, int num_jumps){
    state_sr = seed;
    // seed xoshiro128++ with SplittableRandom, as recommended by authors
    state_xo[0] = (uint32_t) (next_sr() >> 32);
    state_xo[1] = (uint32_t) (next_sr() >> 32);
    state_xo[2] = (uint32_t) (next_sr() >> 32);
    state_xo[3] = (uint32_t) (next_sr() >> 32);
    // Now do jumps 
    do_long_jumps(num_jumps);
}

// Construct RNG state object 
RNGState::RNGState(int seed, int num_jumps){
    seed_rng(seed, num_jumps);
}

RNGState::RNGState(void){
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);
}

uint64_t RNGState::next_sr(){
    uint64_t z = (state_sr += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

// https://prng.di.unimi.it/xoshiro128plusplus.c
uint32_t RNGState::generator(void) {
    const uint32_t result = rotl(state_xo[0] + state_xo[3], 7) + state_xo[0];

    const uint32_t t = state_xo[1] << 9;

    state_xo[2] ^= state_xo[0];
    state_xo[3] ^= state_xo[1];
    state_xo[1] ^= state_xo[2];
    state_xo[0] ^= state_xo[3];

    state_xo[2] ^= t;

    state_xo[3] = rotl(state_xo[3], 11);

    return result;
}


// https://prng.di.unimi.it/xoshiro128plusplus.c
void RNGState::long_jump(){
    static const uint32_t LONG_JUMP[] = { 0xb523952e, 0x0b6f099f, 0xccf5a0ef, 0x1c580662 };

	uint32_t s0 = 0;
	uint32_t s1 = 0;
	uint32_t s2 = 0;
	uint32_t s3 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 32; b++) {
			if (LONG_JUMP[i] & UINT32_C(1) << b) {
				s0 ^= state_xo[0];
				s1 ^= state_xo[1];
				s2 ^= state_xo[2];
				s3 ^= state_xo[3];
			}
			generator();	
		}
		
	state_xo[0] = s0;
	state_xo[1] = s1;
	state_xo[2] = s2;
	state_xo[3] = s3;
}

void RNGState::do_long_jumps(int num_jumps){
    for (size_t i = 0; i < num_jumps; i++)
    {
        long_jump();
    }
}

// Rest of file is original code --------------------------------

/*
 * Generate a uniform random integer in [0, max). Slightly biased.
 */
int RNGState::r_int(uint32_t max){
    uint32_t x = generator();
    uint64_t m = uint64_t(x) * uint64_t(max);
    return (int) (m >> 32);
}

/*
 * Generate a uniform random double in [0, 1). Slightly biased.
 */
double RNGState::r_unif(){
    return 0x1.0p-32 * generator();
}

/*
 * Generate a uniform random integer in [0, max).
 */
int r_int_exact(uint32_t max) {
    uint32_t x = GLOBAL_RNG.generator();
    uint64_t m = uint64_t(x) * uint64_t(max);
    uint32_t l = (uint32_t) m;
    if (l < max) {
        uint32_t t = -max;
        if (t >= max) {
            t -= max;
            if (t >= max)
                t %= max;
        }
        while (l < t) {
            x = (uint32_t) GLOBAL_RNG.generator();
            m = (uint64_t) x * (uint64_t) max;
            l = (uint32_t) m;
        }
    }
    return (int) (m >> 32);
}



// [[Rcpp::export]]
Rcpp::IntegerVector rint1(int n, int max) {
    Rcpp::IntegerVector out(n);
    for (int i = 0; i < n; i++) {
        out[i] = GLOBAL_RNG.r_int(max);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector runif1(int n, int max) {
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; i++) {
        out[i] = GLOBAL_RNG.r_unif();
    }
    return out;
}


// helper
int find_u(double u, int max, vec cum_wgts) {
    int low = 0, high = max - 1;

    if (cum_wgts[0] > u)
        return 0;

    while (high - low > 1) {
        int midpt = std::ceil((high + low) / 2.0);
        if (cum_wgts[midpt] <= u)
            low = midpt;
        else
            high = midpt;
    }

    return high;
}

/*
 * Generate a random integer in [0, cum_wgts.size()) according to normalized cumulative weights.
 */
int RNGState::r_int_wgt(vec cum_wgts) {
    return find_u(r_unif(), cum_wgts.size(), cum_wgts);
}

/* 
 *  Generate a random index of `unnormalized_wgts` with probability proportional to its weight
 * 
 *  Takes a vector of strictly positive weights and returns an index with probability 
 *  proportional to its weight. In other words, it selects index `i` with probability
 *  proporitional to `unnormalized_wgts[i]` 
 *  (or exactly `unnormalized_wgts[i]/sum(unnormalized_wgts)`). This does not support
 *  inputs where some of the weights are zero. This has positive probability of 
 *  returning indices that have weight zero. 
 * 
 * 
 *  @param unnormalized_wgts An arma vector of positive numbers
 * 
 *  @details no Modifications to inputs made
 * 
 *  @returns An integer in [0, `unnormalized_wgts.size()`)
 * 
 */
int RNGState::r_int_unnormalized_wgt(const vec &unnormalized_wgts) {
    // Get the unnormalized cumulative weights 
    arma::vec cum_wgts = arma::cumsum(unnormalized_wgts); 
    // now normalize them
    cum_wgts = cum_wgts / cum_wgts(cum_wgts.size()-1);
    return r_int_wgt(cum_wgts);
}


/*
 * Generate a random integer within a stratum
 */
int r_int_mixstrat(int max, int stratum, double p, vec cum_wgts) {
    double u;
    if (GLOBAL_RNG.r_unif() > p) {
        u = (stratum + GLOBAL_RNG.r_unif()) / max;
    } else {
        u = GLOBAL_RNG.r_unif();
    }

    return find_u(u, max, cum_wgts);
}

/*
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
ivec resample_lowvar(vec wgts) {
    int N = wgts.n_elem;

    double r = GLOBAL_RNG.r_unif() / N;
    double cuml = wgts[0];
    ivec out(N);

    int i = 0;
    for (int n = 0; n < N; n++) {
        double u = r + n / (double) N;
        while (u > cuml) {
            cuml += wgts[++i]; // increment then access
        }
        out[n] = i + 1;
    }

    return out;
}



/*
 * Set RNG seed globally. Be very careful using this, it is not thread safe. 
 */
void global_seed_rng(int seed, int num_jumps) {
    GLOBAL_RNG.seed_rng(seed, num_jumps);
}

