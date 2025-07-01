test_that("pareto domination is calculated correctly", {
    x = matrix(1:5, nrow=1)
    expect_equal(pareto_dominated(x), c(FALSE, rep(TRUE, 4)))
    x = matrix(5:1, nrow=1)
    expect_equal(pareto_dominated(x), c(rep(TRUE, 4), FALSE))

    x = matrix(c(0, 0,
                 1, 0,
                 0, 1), nrow=2)
    expect_equal(pareto_dominated(x), c(FALSE, TRUE, TRUE))

    x = matrix(c(0, 1,
                 0, 0,
                 1, 0), nrow=2)
    expect_equal(pareto_dominated(x), c(TRUE, FALSE, TRUE))

    # remove dupes in right order
    x = matrix(c(0, 0,
                 0, 0,
                 0, 1), nrow=2)
    expect_equal(pareto_dominated(x), c(FALSE, TRUE, TRUE))
})

test_that("short bursts run and improve the score over time", {
    iowa_map <- suppressWarnings(redist_map(iowa, ndists = 4, pop_tol = 0.2))
    plans <- redist_shortburst(iowa_map, scorer_frac_kept(iowa_map),
        max_bursts = 20, verbose = FALSE)

    expect_true(inherits(plans, "redist_plans"))
    expect_equal(dim(get_plans_matrix(plans)), c(99, 21))
    expect_true(all(diff(plans$score) >= 0))

    plans <- redist_shortburst(iowa_map, scorer_frac_kept(iowa_map),
                               max_bursts = 20, maximize = FALSE, verbose = FALSE)
    expect_true(all(diff(plans$score) <= 0))

    plans <- redist_shortburst(iowa_map, scorer_frac_kept(iowa_map),
        max_bursts = 20, thin = 2, verbose = FALSE)
    expect_true(ncol(as.matrix(plans)) == 11)
})


test_that("short bursts work with multiple scorers", {
    iowa_map <- suppressWarnings(redist_map(iowa, ndists = 4, pop_tol = 0.2))
    scorer = cbind(
        comp = scorer_frac_kept(iowa_map),
        dem = scorer_group_pct(iowa_map, dem_08, tot_08, k=2)
    )
    plans <- redist_shortburst(iowa_map, scorer, maximize=c(comp=FALSE, dem=TRUE),
                               max_bursts = 20, verbose = FALSE)

    expect_true(inherits(plans, "redist_plans"))
    expect_equal(dim(get_plans_matrix(plans)), c(99, 21))
    # check frontier is monotone in right direction
    expect_equal(1, cor(attr(plans, "pareto_score"), method="spearman")[1, 2])
})
