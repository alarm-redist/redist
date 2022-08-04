#include "random.h"

std::random_device rd;
std::mt19937 generator(rd());

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
