test_that("SMC runs without errors", {
    N = 20
    expect_silent({
        res = redist.smc(adj, pop, N, 3, pop_tol=0.1, silent=T)
        NULL
    })
    expect_equal(res$nsims, N)
    expect_true(all(res$maxdev <= 0.1))
    expect_equal(res$wgt, rep(1, N)/N)
    expect_equal(range(res$plans), c(1, 3))
})

test_that("Precise population bounds are enforced", {
    res = redist.smc(adj, pop, 20, 3, pop_bounds=c(52e3, 58e3, 60e3), silent=T)
    distr_pop = apply(res$plans, 2, function(x) tapply(pop, x, sum))
    expect_true(all(apply(distr_pop, 2, max) <= 60e3))
    expect_true(all(apply(distr_pop, 2, min) >= 52e3))
})

test_that("SMC checks arguments", {
    expect_error(redist.smc(adj, pop, 10, 3, pop_tol=0.0), "positive")
    expect_error(redist.smc(adj, pop, 10, 3, pop_tol=0.1, compactness=-1), "non-negative")
    expect_error(redist.smc(adj, pop, 10, 3, pop_tol=0.1, seq_alpha=1.5), "0, 1")
    expect_error(redist.smc(adj, pop, 10, 3, pop_tol=0.1, adapt_k_thresh=1.5), "0, 1")
    expect_error(redist.smc(adj, pop, 0, 3, pop_tol=0.1), "positive")
    expect_error(redist.smc(adj, pop, 10, 3, counties=rep(2, length(adj)), pop_tol=0.1), "must run from")
})
