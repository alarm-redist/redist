#include "random.h"

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
 * Generate a random integer in [0, max) according to weights.
 */
int rint(int max, vec cum_wgts) {
    int low = 0, high = max - 1;
    double u = unif(generator);

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
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
ivec resample_lowvar(vec wgts) {
    int N = wgts.n_elem;

    double r = unif(generator) / N;
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
