test_that("redist_mms works", {
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

    nsims <- 10
    out <- redist_mms(fl_map, nsims, nsims %/% 2, init_plan = plans_10[, 1], silent = TRUE)
    par <- redist.parity(as.matrix(out), total_pop = pop)

    expect_equal(range(as.matrix(out)), c(1, 3))
    expect_true(all(par <= 0.1))

    out <- redist_mms(fl_map, 20, 5, thin = 4, init_plan = plans_10[, 1], silent = TRUE)
    expect_equal(ncol(as.matrix(out)), 5L)
})

test_that("Additional constraints work", {
    skip_on_cran()
    iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.05)

    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) %>%
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) %>%
        add_constr_custom(1e6, function(plan, distr) plan[7] == 2)

    plans <- redist_mms(iowa_map, 100, 20, init_plan = iowa$cd_2010,
                               constraints = constr, silent = TRUE)
    skip = seq(1, which.max(by_plan(plans$mcmc_accept, 4)) - 1) # skip to first acceptance
    expect_false(any(as.matrix(plans)[7, -skip] == 2))
})
