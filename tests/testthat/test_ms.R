test_that("redist_mergesplit works", {
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

    nsims <- 10
    out <- redist_mergesplit(fl_map, nsims, nsims %/% 2, init_plan = plans_10[, 1], silent = TRUE)
    par <- redist.parity(as.matrix(out), total_pop = pop)

    expect_equal(range(as.matrix(out)), c(1, 3))
    expect_true(all(par <= 0.1))

    out <- redist_mergesplit(fl_map, 20, 5, thin = 4, init_plan = plans_10[, 1], silent = TRUE)
    expect_equal(ncol(as.matrix(out)), 5L)
})

test_that("Additional constraints work", {
    skip_on_cran()
    iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.05)

    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) %>%
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) %>%
        add_constr_custom(1e6, function(plan, distr) plan[7] == 2)

    plans <- redist_mergesplit(iowa_map, 100, 20, init_plan = iowa$cd_2010,
                               constraints = constr, silent = TRUE)
    skip = seq(1, which.max(by_plan(plans$mcmc_accept, 4)) - 1) # skip to first acceptance
    expect_false(any(as.matrix(plans)[7, -skip] == 2))
})

test_that("redist_mergesplit_parallel works", {
    skip_on_os("windows")
    data(fl25)
    fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1) %>% suppressMessages()

    N <- 20
    chains <- 2

    pl1 <- redist_mergesplit(fl_map, nsims = N, warmup = N/2, chains = chains,
        ncores = 2, silent = TRUE)
    pl2 <- redist_mergesplit(fl_map, nsims = N, chains = chains,
        ncores = 2, warmup = 0, init_name = F, silent = TRUE)
    pl3 <- redist_mergesplit(fl_map, nsims = N, warmup = N/2, chains = chains,
        ncores = 2, return_all = F, init_name = F, silent = TRUE)

    expect_equal(get_n_ref(pl1), chains)
    expect_equal(get_n_ref(pl2), 0)
    expect_equal(get_n_ref(pl3), 0)

    expect_equal(nrow(pl1), 3*chains*(N/2 + 1))
    expect_equal(nrow(pl2), 3*chains*N)
    expect_equal(nrow(pl3), 3*chains)
})

test_that("Parallel runs are reproducible", {
    set.seed(5118)
    pl1 <- redist_mergesplit(fl_map, 100, warmup = 50, chains = 2, silent = TRUE)
    set.seed(5118)
    pl2 <- redist_mergesplit(fl_map, 100, warmup = 50, chains = 2, silent = TRUE)

    # runtime is the only thing that shouldn't be identical
    for (i in 1:2) {
        attr(pl1, "diagnostics")[[i]]$runtime <- NULL
        attr(pl2, "diagnostics")[[i]]$runtime <- NULL
    }

    expect_identical(pl1, pl2)
})
