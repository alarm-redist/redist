test_that("short bursts run and improve the score over time", {
    fl = suppressWarnings(redist_map(fl25, ndists=3, pop_tol=0.2, adj=adj))
    plans = redist_shortburst(fl, scorer_frackept(fl), max_bursts=20, verbose=F)

    expect_true(inherits(plans, "redist_plans"))
    expect_equal(dim(get_plan_matrix(plans)), c(25, 21))
    expect_true(all(diff(plans$score) >= 0))
})
