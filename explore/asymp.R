library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

# MCMC: how does # dists -> autocorrelation -> samp size
ess_mcmc = Vectorize(\(rho) 1 / (2/(1 - rho) - 1))

cor_mcmc = function(revisit=1.0, N=100, alpha=10, beta=1) { # pr(visit the district with max order stat)
    ac = list(
        lengths = 1 + rgeom(N, revisit),
        values =  rbeta(N, alpha, beta)
    ) |>
        inverse.rle() |>
        acf(lag.max=1, plot=FALSE)
    ac$acf[2]
}

# cor really is just 1 - rev_pr
rev_pr = seq(0.1, 1.0, 0.1)
rev_ac = map_dbl(rev_pr, ~ cor_mcmc(., 1e4, alpha=2, beta=1))
plot(rev_pr, rev_ac); abline(a=1, b=-1, col='red')


# SMC: how does # dists -> # ancestors (proxy for ess)
esss_smc = function(n) {
    b = numeric(n)
    b[1] = 1
    for (i in seq(2, n)) {
        b[i] = 1 - exp(-b[i-1])
    }
    b
}


n_max = 400
d = tibble(
    n_distr = seq_len(n_max),
    ess_mcmc = ess_mcmc(1 - 1 / seq_len(n_max)),
    ess_smc = esss_smc(n_max)
)

d |>
    pivot_longer(-n_distr, names_prefix="ess_", names_to="algo", values_to="ess") |>
    mutate(min_s = 400 / ess) |>
ggplot(aes(n_distr, min_s, color=algo)) +
    geom_line(linewidth=1.0) +
    geom_hline(yintercept=5000, lty="dashed") +
    geom_hline(yintercept=10000, lty="dashed") +
    geom_hline(yintercept=40000, lty="dashed") +
    geom_vline(xintercept=2, lty="dashed") +
    geom_vline(xintercept=3, lty="dashed") +
    geom_vline(xintercept=7, lty="dashed") +
    geom_vline(xintercept=28, lty="dashed") +
    geom_vline(xintercept=52, lty="dashed") +
    scale_x_log10() +
    scale_y_log10("Minimum sample size needed") +
    theme_bw()

plot(log(d$ess_smc), type='l')
curve(log(1/x), add=T, from=1, to=400, col='red')
