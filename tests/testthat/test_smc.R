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

test_that("SMC checks arguments", {
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.0), "positive")
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.1, compactness=-1), "non-negative")
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.1, seq_alpha=1.5), "0, 1")
    expect_error(redist.smc(g, pop, 10, 3, popcons=0.1, adapt_k_thresh=1.5), "0, 1")
    expect_error(redist.smc(g, pop, 0, 3, popcons=0.1), "positive")
    expect_error(redist.smc(g, pop, 10, 3, counties=rep(2, length(g)), popcons=0.1), "must run from")
})
