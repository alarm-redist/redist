skip_on_os("solaris")
skip_on_os("linux")

test_that("SMC runs without errors", {
    N <- 20
    res <- redist_smc(fl_map, N, silent = TRUE)

    expect_s3_class(res, "redist_plans")
})

test_that("County constraint works", {
    iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.05)
    plans <- redist_smc(iowa_map, 50, counties = region, silent = TRUE)
    splits <- redistmetrics::splits_admin(plans, iowa_map, region)
    expect_true(all(splits <= 3L))
    expect_true(all(apply(get_plans_matrix(plans), 2,
        function(x) all(contiguity(iowa_map$adj, x) == 1))))

    region2 <- iowa$region
    region2[25] <- NA
    expect_error(redist_smc(iowa_map, 50, counties = region2, silent = TRUE),
        "missing values")
})

test_that("Single-precinct counties work", {
    bb <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0)))))
    grid <- sf::st_make_grid(bb, n = 4)
    tb <- sf::st_as_sf(grid) %>%
        rename(geometry = x) %>%
        mutate(pop = 1, counties = row_number())
    box_map <- redist_map(tb, ndists = 4, pop_tol = 0.25, total_pop = pop) %>%
        suppressWarnings() %>%
        suppressMessages()

    test <- redist_smc(box_map, 50, counties = counties, silent = TRUE)
    expect_s3_class(test, "redist_plans")
})

test_that("Not egregiously incorrect sampling accuracy (5-prec)", {
    skip_on_cran()
    set.seed(1935)

    g <- list(c(1L, 4L), c(0L, 2L, 4L), c(1L, 3L, 4L), c(2L, 4L), c(0L, 1L, 2L, 3L))
    g_pop <- c(2, 1, 1, 1, 1)
    map <- redist_map(pop = g_pop, ndists = 2, pop_tol = 0.5, adj = g)
    out <- redist_smc(map, 20e3, compactness = 0,
                      split_params = list(adapt_k_thresh = 0.99), ncores = 1L,
                      resample = FALSE, silent = TRUE)
    types <- apply(as.matrix(out), 2, function(x) 1L + (x[1] == x[2]))
    wgts <- weights(out)
    avg <- weighted.mean(types, wgts)
    se <- sqrt(sum((types - avg)^2*(wgts/sum(wgts))^2))
    zscores <- (avg - 1.5) / se
    expect_true(abs(zscores) <= 3)
})

test_that("Not egregiously incorrect sampling accuracy (25-prec)", {
    skip_on_cran()
    set.seed(1935)

    ref_plans <- plans_10[, redist.parity(plans_10, pop) <= 0.01]
    log_st_ref <- round(log_st_map(adj, ref_plans, rep(1L, 25), 3L), 5)

    out <- redist_smc(set_pop_tol(fl_map, 0.01), 6000, compactness = 0,
                      ncores = 1L, control = list(weight_type = "simple"),
                      split_params = list(adapt_k_thresh = .99999), seq_alpha = 1L,
                      resample = FALSE, silent = TRUE) %>%
        suppressWarnings() # efficiency
    log_st <- round(log_st_map(adj, as.matrix(out), rep(1L, 25), 3L), 5)
    types <- match(log_st, log_st_ref)

    wgts <- weights(out)
    avgs <- sapply(seq_along(log_st_ref), function(i) weighted.mean(types == i, wgts))
    ses <- sapply(seq_along(log_st_ref), function(i) {
        sqrt(sum(((types == i) - avgs[i])^2*(wgts/sum(wgts))^2))
    })
    zscores <- (avgs - (1/length(log_st_ref)))/ses
    expect_true(all(abs(zscores) <= 5))
})

test_that("Labeling accounted for", {
    skip_on_cran()
    set.seed(1935)

    g <- list(1:2, c(0L, 3L), c(0L, 3L, 4L), c(1L, 2L, 5L), c(2L, 5L, 6L),
              c(3L, 4L, 7L), c(4L, 7L), 5:6)
    map <- redist_map(pop = rep(1, 8), ndists = 4, pop_tol = 0.05, adj = g)
    out <- redist_smc(map, 10e3, ncores = 1L,
                      split_params = list(adapt_k_thresh = 1L),
                      resample = FALSE, silent = TRUE)
    types = apply(as.matrix(out), 2, function(x) {
        paste(vctrs::vec_group_id(x), collapse="")
    }) |>
        vctrs::vec_group_id()
    wgts <- weights(out)
    avgs <- sapply(1:5, function(i) weighted.mean(types == i, wgts))
    ses <- sapply(1:5, function(i) {
        sqrt(sum(((types == i) - 0.2)^2*(wgts/sum(wgts))^2))
    })
    zscores <- (avgs - 0.2) / ses
    expect_true(all(abs(zscores) <= 3))
})



test_that("Partial sampling works accurately", {
    skip_on_cran()
    set.seed(1935)
    nsims <- 1000
    compactness <- 1L

    # check that sampling all the way and partial sampling return about the same

    # sample all the way
    out_all_the_way <- redist_smc(set_pop_tol(fl_map, 0.01), nsims, compactness = compactness,
                                  split_params = list(adapt_k_thresh = .99), seq_alpha = 1L,
                                  ncores = 1L, # control = list(weight_type = "simple"),
                                  resample = F, silent = TRUE) %>%
        suppressWarnings() # efficiency

    # partially sample
    out1 <- redist_smc(set_pop_tol(fl_map, 0.01), nsims, compactness = compactness,
                       n_steps = 1, split_params = list(adapt_k_thresh = .99), seq_alpha = 1L,
                       ncores = 1L, # control = list(weight_type = "simple"),
                       resample = F, silent = TRUE) %>%
        suppressWarnings() # efficiency
    out2 <- redist_smc(set_pop_tol(fl_map, 0.01), nsims, compactness = compactness,
                       init_particles = out1,
                       ncores = 1L, #control = list(weight_type = "simple"),
                       split_params = list(adapt_k_thresh = .99),
                       seq_alpha = 1L, resample = F, silent = TRUE) %>%
        suppressWarnings() # efficiency

    # get the counts of the plans in the two samples
    all_the_way_counts_df <- get_plan_counts(
        get_plans_matrix(out_all_the_way), attr(fl_map, "ndists")
    ) |>
        arrange(plan_string)

    partial_counts_df <- get_plan_counts(
        get_plans_matrix(out2), attr(fl_map, "ndists")
    ) |>
        arrange(plan_string)

    test_result <- chisq.test(
        all_the_way_counts_df$count,
        partial_counts_df$count,
        simulate.p.value = TRUE, B = 10000
    ) |>
        suppressWarnings() # efficiency

    expect_true(test_result$p.value > .001)

})



test_that("Additional constraints work", {
    iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.05)

    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) %>%
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) %>%
        add_constr_custom(1e2, function(plan, distr) plan[7] == 2)

    plans <- redist_smc(iowa_map, 100, constraints = constr, silent = TRUE)
    expect_false(any( (as.matrix(plans)[7, ] - 1L) == 2))
})



test_that("Precise population bounds are enforced", {
    map2 <- fl_map
    attr(map2, "pop_bounds") <- c(52e3, 58e3, 60e3)
    res <- redist_smc(map2, 20, silent = TRUE)
    distr_pop <- apply(as.matrix(res), 2, function(x) tapply(pop, x, sum))
    expect_true(all(apply(distr_pop, 2, max) <= 60e3))
    expect_true(all(apply(distr_pop, 2, min) >= 52e3))
})

test_that("SMC checks arguments", {
    expect_error(redist_smc(fl_map, 10, compactness = -1), "non-negative")
    expect_error(redist_smc(fl_map, 10, seq_alpha = 1.5), "0, 1")
    expect_error(redist_smc(fl_map, 10, split_params = list(adapt_k_thresh = 1.5)), "0, 1")
    expect_error(redist_smc(fl_map, 0), "positive")
})

test_that("Parallel runs are reproducible", {
    # need to make sure 1 thread and process
    set.seed(5118)
    pl1 <- redist_smc(fl_map, 100, runs = 2, silent = TRUE,
                      ncores = 1L, control = list(nproc = 1L))
    set.seed(5118)
    pl2 <- redist_smc(fl_map, 100, runs = 2, silent = TRUE,
                      ncores = 1L, control = list(nproc = 1L))

    # runtime related is the only thing that shouldn't be identical
    for (i in 1:2) {
        attr(pl1, "diagnostics")[[i]]$runtime <- NULL
        attr(pl2, "diagnostics")[[i]]$runtime <- NULL
        attr(pl1, "total_runtime") <- NULL
        attr(pl2, "total_runtime") <- NULL
    }

    expect_identical(pl1, pl2)
})
