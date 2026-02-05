test_that("skinny flips works", {

    # Compute pop_bounds from pop_tol = 0.15
    ndists <- 3
    parity <- sum(fl25$pop) / ndists
    pop_bounds <- c(parity * 0.85, parity, parity * 1.15)

    plans <- skinny_flips(adj = adj, init_plan = plans_10[, 1],
        total_pop = fl25$pop, pop_bounds = pop_bounds,
        nsims = 5, eprob = 0.05, lambda = 1, constraints = redist_constr())


    expect_true("matrix" %in% class(plans))
    expect_true(min(plans) == 1)

})
