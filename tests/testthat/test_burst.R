test_that("short bursts run and improve the score over time", {
    iowa_map <- suppressWarnings(redist_map(iowa, ndists = 4, pop_tol = 0.2))
    plans <- redist_shortburst(iowa_map, scorer_frac_kept(iowa_map),
        max_bursts = 20, verbose = F)

    expect_true(inherits(plans, "redist_plans"))
    expect_equal(dim(get_plans_matrix(plans)), c(99, 21))
    expect_true(all(diff(plans$score) >= 0))

    plans <- redist_shortburst(iowa_map, scorer_frac_kept(iowa_map),
        max_bursts = 20, thin = 2, verbose = F)
    expect_true(ncol(as.matrix(plans)) == 11)
})
