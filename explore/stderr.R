library(dplyr)
library(tidyr)
library(purrr)

data(iowa)

ia = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.02)

plans = redist_smc(ia, 2000, runs=8, seq_alpha=0.5, silent=T) %>%
    mutate(dem = group_frac(ia, dem_08, dem_08 + rep_08),
           comp = distr_compactness(ia),
           splits = county_splits(ia, region)) %>%
    number_by(dem)

plot(plans, dem, sort=F, geom="boxplot")
hist(plans, comp)
summary(plans)

calc_se = function(x1) {
    run_means = tapply(x1, plans$chain[plans$district == 1], mean) %>%
        `names<-`(NULL)
    std_err = sd(run_means)

    x1 = x1[coalesce(plans$chain[plans$district == 1] == 1, F)]
    stopifnot(all.equal(mean(x1), run_means[1]))

    N = length(x1)
    d = tibble(ancestor = factor(attr(plans, "diagnostics")[[1]]$ancestors, levels=1:N),
               resid = x1 - mean(x1))
               # resid = x1)
    cross_sum = expand_grid(rename_with(d, \(x) paste0(x, "_i")),
                rename_with(d, \(x) paste0(x, "_j"))) %>%
        filter(ancestor_i != ancestor_j) %>%
        summarize(sum(resid_i * resid_j)) %>%
        pull()

    # var_est = mean((x1 - mean(x1)))^2 - exp((4 - 3) * log(N) - (4 - 1) * log(N - 1)) * cross_sum
    # var_est = mean(d$resid)^2 - exp((4 - 3) * log(N) - (4 - 1) * log(N - 1)) * cross_sum
    var_est = - (N/(N-1))^2 * (cross_sum / N / (N-1))
    # var_est = mean(d$resid)^2 - (N/(N-1))^2 * (cross_sum / N / (N-1)) #* mean(1/get_sampling_info(plans)$diagnostics[[1]]$accept_rate)

    n_eff = attr(plans, "diagnostics")[[1]]$n_eff
    tibble(true = std_err,
           est1 = sqrt(var_est),
           est2 = sd(x1) / sqrt(n_eff))
}

qs = seq(0.02, 0.98, 0.02)
x1 = plans$comp[plans$district == 1]
refs = quantile(x1, qs)

calc_se(x1)
ses = map_dfr(refs, ~ calc_se(x1 > .))
ggplot(ses, aes(refs)) +
    geom_line(aes(y=true), color="red", lty="dashed") +
    geom_line(aes(y=est1)) +
    geom_line(aes(y=est2), color="grey")


res = map_dfr(cli::cli_progress_along(1:100), function(i) {
    plans = redist_smc(ia, 1000, cores=4, seq_alpha=0.5, silent=T) %>%
        mutate(dem = group_frac(ia, dem_08, dem_08 + rep_08)) %>%
        number_by(dem)

    x1 = plans$dem[plans$district == 1][-1]
    x1 = x1 >= plans$dem[1]

    N = length(x1)
    d = tibble(ancestor = factor(attr(plans, "diagnostics")[[1]]$ancestors, levels=1:N),
               resid = x1 - mean(x1))
    cross_sum = expand_grid(rename_with(d, \(x) paste0(x, "_i")),
                            rename_with(d, \(x) paste0(x, "_j"))) %>%
        filter(ancestor_i != ancestor_j) %>%
        summarize(sum(resid_i * resid_j)) %>%
        pull()

    var_est_ancestor = - (N/(N-1))^2 * (cross_sum / N / (N-1))
    var_est_neff = var(x1) / attr(plans, "diagnostics")[[1]]$n_eff

    tibble(mean = mean(x1),
           se1 = sqrt(var_est_ancestor),
           se2 = sqrt(var_est_neff))
})

sd(res$mean)
mean(res$se1)
mean(res$se2)
mean(res$se1) + 2*sd(res$se1)/sqrt(200)


if (FALSE) {
    library(tictoc)
    library(ggplot2)

    pa = alarmdata::alarm_50state_map("PA")
    {
        tic()
        plans_smc = redist_smc(pa, 1000, counties=county,
                               pop_temper=0.01, seq_alpha=0.65,
                               n_runs=4, cores=2, verbose=TRUE)
        toc()

        tic()
        plans_ms = redist_mergesplit_parallel(pa, 2000, chains=4, warmup=1000,
                                              counties=county, verbose=FALSE)
        is
        toc()
    }

    plans_smc = plans_smc %>%
        mutate(comp = comp_frac_kept(., pa),
               splits = county_splits(pa, county),
               dem = group_frac(pa, ndv, ndv + nrv),
               min = group_frac(pa, vap - vap_white, vap),
               egap = part_egap(., pa, ndv, nrv)) %>%
        number_by(dem)
    plans_ms = plans_ms %>%
        mutate(comp = comp_frac_kept(., pa),
               splits = county_splits(pa, county),
               dem = group_frac(pa, ndv, ndv + nrv),
               min = group_frac(pa, vap - vap_white, vap),
               egap = part_egap(., pa, ndv, nrv)) %>%
        number_by(dem)

    summary(plans_smc)
    summary(plans_ms)

    as_tibble(plans_ms) %>%
        filter(district == 1) %>%
        ggplot(aes((as.integer(draw) - 1) %% 1000, comp, color=chain, group=chain)) +
        geom_line()

    ggplot(NULL) +
        geom_density(aes(plans_diversity(plans_smc)), color='red') +
        geom_density(aes(plans_diversity(plans_ms)))

    sapply(get_sampling_info(plans_smc)$diagnostics, \(x) n_distinct(x$ancestors))
    sapply(get_sampling_info(plans_smc)$diagnostics, \(x) sd(table(x$ancestors)))

    plot(pa, rowMeans(as.matrix(plans_smc)))

    plans_smc %>%
        as_tibble() %>%
        filter(district == 17) %>%
        mutate(x = dem <= dem[1]) %>%
        with(., c(mean(x), diag_rhat(x, chain)))
}
