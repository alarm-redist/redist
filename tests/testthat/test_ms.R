test_that("redist_mergesplit works", {
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

    nsims <- 10
    out <- redist_mergesplit(fl_map, nsims, nsims %/% 2, init_plan = plans_10[, 1], silent = TRUE)
    par <- redist.parity(as.matrix(out), total_pop = pop)

    expect_equal(range(as.matrix(out)), c(1, 3))
    expect_true(all(par <= 0.1))

    out <- redist_mergesplit(fl_map, 4, 5, thin = 4, init_plan = plans_10[, 1], silent = TRUE)
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
    # subtract 1 since plans are zero indexed in c++
    expect_false(any( (as.matrix(plans)[7, -skip] - 1L) == 2))
})



test_that("Thresholding constraints work", {
    iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.05)

    # ensure that polk and story county are always in the same region
    polk_precint <- which(iowa_map$name == "Polk")
    story_precint <- which(iowa_map$name == "Story")

    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) %>%
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) %>%
        add_constr_custom_plan(1, function(plan, seats, num_regions){
            # return 0 if in the same region
            if(plan[polk_precint] == plan[story_precint]){
                return(0)
            }else{
                return(1)
            }
        }, thresh = .5)

    plans <- redist_mergesplit(iowa_map, nsims = 100, warmup = 50,
                               constraints = constr, silent = TRUE)

    expect_true(
        all(get_plans_matrix(plans)[polk_precint, ] == get_plans_matrix(plans)[story_precint, ])
    )

    # now ensure polk and story are never in the same region
    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) %>%
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) %>%
        add_constr_custom_plan(1, function(plan, seats, num_regions){

            if(num_regions == 1 || plan[polk_precint] != plan[story_precint] ){
                return(0)
            }else{
                return(1)
            }
            # return 0 if in the same region

        }, thresh = .5)

    plans <- redist_mergesplit(iowa_map, nsims = 100, warmup = 50,
                               constraints = constr, silent = TRUE)

    expect_true(
        all(get_plans_matrix(plans)[polk_precint, ] != get_plans_matrix(plans)[story_precint, ])
    )
})

test_that("redist_mergesplit with multiple chains works", {
    skip_on_os("windows")
    data(fl25)
    ndists <- 3
    fl_map <- redist_map(fl25, ndists = ndists, pop_tol = 0.1) %>% suppressMessages()

    N <- 20
    chains <- 2

    pl1 <- redist_mergesplit(fl_map, nsims = N, warmup = N/2, chains = chains,
                                      ncores = 2, silent = TRUE)
    pl2 <- redist_mergesplit(fl_map, nsims = N, chains = chains,
                                      ncores = 2, warmup = 0, init_name = FALSE, silent = TRUE)
    pl3 <- redist_mergesplit(fl_map, nsims = N, warmup = N/2, chains = chains,
                                      ncores = 2, return_all = FALSE, init_name = FALSE, silent = TRUE)

    expect_equal(get_n_ref(pl1), chains)
    expect_equal(get_n_ref(pl2), 0)
    expect_equal(get_n_ref(pl3), 0)

    # we expect N * chains + num_refs * chains (for the reference plans )
    expect_equal(nrow(pl1), ndists*chains*N + chains*ndists)
    expect_equal(nrow(pl2), ndists*chains*N)
    expect_equal(nrow(pl3), ndists*chains*N)
})

test_that("Parallel runs are reproducible", {
    set.seed(5118)
    pl1 <- redist_mergesplit(
        fl_map, 100, warmup = 50, chains = 2, silent = TRUE,
        control = list(init_ncores = 1L))

    set.seed(5118)
    pl2 <- redist_mergesplit(
        fl_map, 100, warmup = 50, chains = 2, silent = TRUE,
        control = list(init_ncores = 1L))

    # runtime is the only thing that shouldn't be identical
    for (i in 1:2) {
        attr(pl1, "diagnostics")[[i]]$runtime <- NULL
        attr(pl2, "diagnostics")[[i]]$runtime <- NULL
        attr(pl1, "total_runtime") <- NULL
        attr(pl2, "total_runtime") <- NULL
    }

    expect_identical(pl1, pl2)
})
