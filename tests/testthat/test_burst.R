test_that("short bursts run and improve the score over time", {
    iowa = suppressWarnings(redist_map(iowa, ndists=4, pop_tol=0.2))
    plans = redist_shortburst(iowa, scorer_frac_kept(iowa), max_bursts=20, verbose=F)

    expect_true(inherits(plans, "redist_plans"))
    expect_equal(dim(get_plans_matrix(plans)), c(99, 21))
    expect_true(all(diff(plans$score) >= 0))
})
