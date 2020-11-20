context("Sequential Monte Carlo")


test_that("SMC runs without errors", {
    N = 20
    expect_silent({
        res = redist.smc(g, pop, N, 3, popcons=0.1, silent=T)
        NULL
    })
    expect_equal(res$nsims, N)
    expect_true(all(res$maxdev <= 0.1))
    expect_equal(res$wgt, rep(1, N)/N)
    expect_equal(range(res$cdvec), c(1, 3))
})

test_that("Subset population bounds are calculated correctly", {
    bounds = redist.subset_bounds(pop, algdat.p10$cdmat[,1] %in% 1:2, 3, 2, 0.1)
    expect_equal(bounds, c(52513, 56666, 64182))
})

test_that("Precise population bounds are enforced", {
    res = redist.smc(g, pop, 20, 3, pop_bounds=c(52e3, 58e3, 60e3), silent=T)
    distr_pop = apply(res$cdvec, 2, function(x) tapply(pop, x, sum))
    expect_true(all(apply(distr_pop, 2, max) <= 60e3))
    expect_true(all(apply(distr_pop, 2, min) >= 52e3))
})

test_that("SMC checks arguments", {
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.0), "positive")
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.1, compactness=-1), "non-negative")
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.1, seq_alpha=1.5), "0, 1")
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.1, adapt_k_thresh=1.5), "0, 1")
    expect_error(redist.smc(g, pop, 0, 3, popcons=0.1), "positive")
    expect_error(redist.smc(g, pop, 10, 3, counties=rep(2, length(g)), popcons=0.1), "must run from")
})
