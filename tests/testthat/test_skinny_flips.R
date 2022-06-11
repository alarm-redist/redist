test_that("skinny flips works", {

    plans <- skinny_flips(adj = adj, init_plan = plans_10[, 1],
        total_pop = fl25$pop, pop_tol = 0.15,
        nsims = 5, eprob = 0.05, lambda = 1, constraints = redist_constr())


    expect_true("matrix" %in% class(plans))
    expect_true(min(plans) == 1)

})
