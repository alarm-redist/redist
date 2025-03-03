map = alarmdata::alarm_50state_map("CT") |>
    set_pop_tol(0.05)

# pl1 = redist_mergesplit_parallel(map, 2200, warmup=200, chains=2)
# pl2 = redist_mergesplit(map, 2200, warmup=200)
perims = prep_perims(map)

pl_ms = redist_mergesplit(map, 8200, warmup=200, thin=2, init_plan="sample", adapt_k_thresh=1) |>
    mutate(polsby = comp_polsby(pl(), map, perim_df=perims, ncores=4),
           dem = group_frac(map, ndv, ndv + nrv)) |>
    subset_sampled() |>
    summarize(polsby = mean(polsby),
              e_dem = sum(dem > 0.5),
              .by=any_of(c("draw", "chain")))
summary(pl_ms)
redist.plot.trace(pl_ms, polsby)
redist.plot.trace(pl_ms, e_dem)
redist.plot.plans(pl_ms, 4000, map)

pl_smc = redist_smc(map, 2000, num_threads_per_process=2, runs=2, pop_temper=0.007, adapt_k_thresh=1) |>
    mutate(polsby = comp_polsby(pl(), map, perim_df=perims, ncores=4),
           dem = group_frac(map, ndv, ndv + nrv)) |>
    subset_sampled() |>
    summarize(polsby = mean(polsby),
              e_dem = sum(dem > 0.5),
              .by=c("draw", "chain"))
summary(pl_smc)
redist.plot.plans(pl_smc, 4000, map)

qqplot(pl_smc$polsby, pl_ms$polsby, cex=0.1); abline(a=0, b=1, col='red');
redist_ci(pl_smc, 8==e_dem) * 4e3
redist_ci(pl_ms, 8==e_dem) * 4e3
redist_ci(pl_smc, polsby)
redist_ci(pl_ms, polsby)

ks.test(pl_smc$polsby, pl_ms$polsby)
t.test(pl_smc$polsby, pl_ms$polsby)
t.test(pl_smc$polsby^2, pl_ms$polsby^2)
t.test(pl_smc$polsby^3, pl_ms$polsby^3)
t.test(with(pl_smc, (polsby - mean(polsby))^2), with(pl_ms, (polsby - mean(polsby))^2))
t.test(pl_smc$polsby < 0.18, pl_ms$polsby < 0.18)
t.test(pl_smc$polsby < 0.16, pl_ms$polsby < 0.16)
t.test(pl_smc$polsby > 0.23, pl_ms$polsby > 0.23)
t.test(pl_smc$polsby > 0.20 & pl_smc$polsby < 0.21,
       pl_ms$polsby > 0.20 & pl_ms$polsby < 0.21)

hist.default(pl_smc$polsby, breaks=seq(0.1, 0.6, 0.01), col="#70000070")
hist.default(rep(pl_ms$polsby, 1), breaks=seq(0.1, 0.6, 0.01), col="#00007070", add=TRUE)

hist.default(pl_smc$e_dem, breaks=0:10, col="#70000070")
hist.default(rep(pl_ms$e_dem, 1), breaks=0:10, col="#00007070", add=TRUE)
