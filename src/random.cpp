#include "random.h"

std::random_device rd;


/* This is a fixed-increment version of Java 8's SplittableRandom generator
 See http://dx.doi.org/10.1145/2714064.2660195 and
 http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html
 Written in 2015 by Sebastiano Vigna (vigna@acm.org)
 [Public Domain]
 */
static uint64_t state_sr; /* The state can be seeded with any value. */
uint64_t next_sr() {
    uint64_t z = (state_sr += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}


/* This is xoshiro128++ 1.0.
 Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 [Public domain]
 */
static inline uint32_t rotl(const uint32_t x, int k) {
    return (x << k) | (x >> (32 - k));
}

static uint32_t state_xo[4] = {rd(), rd(), rd(), rd()};

uint32_t generator(void) {
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
// Rest of file is original code --------------------------------


/*
 * Set RNG seed
 */
void seed_rng(int seed) {
    state_sr = seed;
    // seed xoshiro128++ with SplittableRandom, as recommended by authors
    state_xo[0] = (uint32_t) (next_sr() >> 32);
    state_xo[1] = (uint32_t) (next_sr() >> 32);
    state_xo[2] = (uint32_t) (next_sr() >> 32);
    state_xo[3] = (uint32_t) (next_sr() >> 32);
}

/*
 * Generate a uniform random integer in [0, max).
 */
int r_int_exact(uint32_t max) {
    uint32_t x = generator();
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
            x = (uint32_t) generator();
            m = (uint64_t) x * (uint64_t) max;
            l = (uint32_t) m;
        }
    }
    return (int) (m >> 32);
}

/*
 * Generate a uniform random integer in [0, max). Slightly biased.
 */
int r_int(uint32_t max) {
    uint32_t x = generator();
    uint64_t m = uint64_t(x) * uint64_t(max);
    return (int) (m >> 32);
}


/*
 * Generate a uniform random double in [0, 1). Slightly biased.
 */
double r_unif() {
    return 0x1.0p-32 * generator();
}

// [[Rcpp::export]]
Rcpp::IntegerVector rint1(int n, int max) {
    Rcpp::IntegerVector out(n);
    for (int i = 0; i < n; i++) {
        out[i] = r_int(max);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector runif1(int n, int max) {
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; i++) {
        out[i] = r_unif();
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
 * Generate a random integer in [0, max) according to weights.
 */
int r_int_wgt(int max, vec cum_wgts) {
    return find_u(r_unif(), max, cum_wgts);
}

/*
 * Generate a random integer within a stratum
 */
int r_int_mixstrat(int max, int stratum, double p, vec cum_wgts) {
    double u;
    if (r_unif() > p) {
        u = (stratum + r_unif()) / max;
    } else {
        u = r_unif();
    }

    return find_u(u, max, cum_wgts);
}

/*
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
ivec resample_lowvar(vec wgts) {
    int N = wgts.n_elem;

    double r = r_unif() / N;
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
