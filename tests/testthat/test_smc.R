skip_on_os("solaris")
skip_on_os("linux")

test_that("SMC runs without errors", {
    N <- 20
    res <- redist_smc(fl_map, N, silent = TRUE)

    expect_s3_class(res, "redist_plans")
})

test_that("County constraint works", {
    iowa_map <- suppressMessages(redist_map(iowa, ndists = 4, pop_tol = 0.05))
    plans <- redist_smc(iowa_map, 50, counties = region, silent = TRUE)
    splits <- redistmetrics::splits_admin(plans, iowa_map, region)
    expect_true(all(splits <= 3L))
    expect_true(all(apply(get_plans_matrix(plans), 2, function(x) {
        all(contiguity(iowa_map$adj, x) == 1)
    })))

    region2 <- iowa$region
    region2[25] <- NA
    expect_error(redist_smc(iowa_map, 50, counties = region2, silent = TRUE), "missing values")
})

test_that("Single-precinct counties work", {
    bb <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0)))))
    grid <- sf::st_make_grid(bb, n = 4)
    tb <- sf::st_as_sf(grid) |> rename(geometry = x) |> mutate(pop = 1, counties = row_number())
    box_map <- redist_map(tb, ndists = 4, pop_tol = 0.25, total_pop = pop) |>
        suppressWarnings() |>
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
    out <- redist_smc(
        map,
        20e3,
        compactness = 0,
        ncores = 1L,
        resample = FALSE,
        silent = TRUE
    )
    types <- apply(as.matrix(out), 2, function(x) 1L + (x[1] == x[2]))
    wgts <- weights(out)
    avg <- weighted.mean(types, wgts)
    se <- sqrt(sum((types - avg)^2 * (wgts / sum(wgts))^2))
    zscores <- (avg - 1.5) / se
    expect_true(abs(zscores) <= 3)
})

test_that("Not egregiously incorrect sampling accuracy (25-prec)", {
    skip_on_cran()
    set.seed(1935)

    ref_plans <- plans_10[, redist.parity(plans_10, pop) <= 0.01]
    log_st_ref <- round(log_st_map(adj, ref_plans, rep(1L, 25), 3L), 5)

    out <- redist_smc(
        set_pop_tol(fl_map, 0.01),
        6000,
        compactness = 0,
        ncores = 1L,
        seq_alpha = 1L,
        resample = FALSE,
        silent = TRUE
    ) |>
        suppressWarnings() # efficiency
    log_st <- round(log_st_map(adj, as.matrix(out), rep(1L, 25), 3L), 5)
    types <- match(log_st, log_st_ref)

    wgts <- weights(out)
    avgs <- sapply(seq_along(log_st_ref), function(i) weighted.mean(types == i, wgts))
    ses <- sapply(seq_along(log_st_ref), function(i) {
        sqrt(sum(((types == i) - avgs[i])^2 * (wgts / sum(wgts))^2))
    })
    zscores <- (avgs - (1 / length(log_st_ref))) / ses
    expect_true(all(abs(zscores) <= 5))
})

test_that("Labeling accounted for", {
    skip_on_cran()
    set.seed(1935)

    g <- list(1:2, c(0L, 3L), c(0L, 3L, 4L), c(1L, 2L, 5L), c(2L, 5L, 6L),
              c(3L, 4L, 7L), c(4L, 7L), 5:6)
    map <- redist_map(pop = rep(1, 8), ndists = 4, pop_tol = 0.05, adj = g)
    out <- redist_smc(
        map,
        10e3,
        ncores = 1L,
        resample = FALSE,
        silent = TRUE
    )
    types <- apply(as.matrix(out), 2, function(x) {
        paste(vctrs::vec_group_id(x), collapse="")
    }) |>
        vctrs::vec_group_id()
    wgts <- weights(out)
    avgs <- sapply(1:5, function(i) weighted.mean(types == i, wgts))
    ses <- sapply(1:5, function(i) {
        sqrt(sum(((types == i) - 0.2)^2 * (wgts / sum(wgts))^2))
    })
    zscores <- (avgs - 0.2) / ses
    expect_true(all(abs(zscores) <= 3))
})


test_that("Partial sampling works accurately", {
    skip_on_cran()
    set.seed(1935)
    nsims <- 10000
    compactness <- 1L

    # check that sampling all the way and partial sampling return about the same

    # sample all the way
    out_all_the_way <- redist_smc(
        set_pop_tol(fl_map, 0.01),
        nsims,
        compactness = compactness,
        seq_alpha = 1L,
        ncores = 1L,
        resample = FALSE,
        silent = TRUE
    ) |>
        suppressWarnings() # efficiency

    # partially sample
    out1 <- redist_smc(
        set_pop_tol(fl_map, 0.01),
        nsims,
        compactness = compactness,
        n_steps = 1,
        seq_alpha = 1L,
        ncores = 1L,
        resample = FALSE,
        silent = TRUE
    ) |>
        suppressWarnings() # efficiency
    out2 <- redist_smc(
        set_pop_tol(fl_map, 0.01),
        nsims,
        compactness = compactness,
        init_particles = out1,
        ncores = 1L,
        seq_alpha = 1L,
        resample = FALSE,
        silent = TRUE
    ) |>
        suppressWarnings() # efficiency

    # get the counts of the plans in the two samples
    all_the_way_counts_df <- get_plan_counts(
        get_plans_matrix(out_all_the_way),
        attr(fl_map, "ndists"),
        return_plan_strings = TRUE
    ) |>
        arrange(plan_string)

    partial_counts_df <- get_plan_counts(
        get_plans_matrix(out2),
        attr(fl_map, "ndists"),
        return_plan_strings = TRUE
    ) |>
        arrange(plan_string)

    test_result <- chisq.test(
        all_the_way_counts_df$count,
        partial_counts_df$count,
        simulate.p.value = TRUE,
        B = 10000
    ) |>
        suppressWarnings() # efficiency

    expect_true(test_result$p.value > 0.001)
})


test_that("Additional constraints work", {
    iowa_map <- suppressMessages(redist_map(iowa, ndists = 4, pop_tol = 0.05))
    constr <- redist_constr(iowa_map) |>
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) |>
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) |>
        add_constr_custom(1e2, function(plan, distr) plan[7] == 2)

    # do district only splitting to avoid putting precinct 7 in district 2
    plans <- redist_smc(iowa_map, 100, constraints = constr, silent = TRUE,
                        split_params = list(splitting_schedule = "split_district_only"))
    expect_false(any((as.matrix(plans)[7, ] - 1L) == 2))
})


test_that("Thresholding constraints work", {
    iowa_map <- suppressMessages(redist_map(iowa, ndists = 4, pop_tol = 0.05))
    # ensure that polk and story county are always in the same region
    polk_precint <- which(iowa_map$name == "Polk")
    story_precint <- which(iowa_map$name == "Story")

    constr <- redist_constr(iowa_map) |>
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) |>
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) |>
        add_constr_custom_plan(
            1,
            function(plan, seats, num_regions) {
                # return 0 if in the same region
                if (plan[polk_precint] == plan[story_precint]) {
                    return(0)
                } else {
                    return(1)
                }
            },
            thresh = 0.5
        )

    plans <- redist_smc(iowa_map, 100, constraints = constr, silent = TRUE,
                        sampling_space = "spanning_forest")

    expect_true(
        all(get_plans_matrix(plans)[polk_precint, ] == get_plans_matrix(plans)[story_precint, ])
    )

    # now ensure polk and story are never in the same region
    constr <- redist_constr(iowa_map) |>
        add_constr_grp_hinge(5, dem_08, tot_08, c(0.5, 0.6)) |>
        add_constr_grp_hinge(5, bvap + hvap, vap, c(0.5, 0)) |>
        add_constr_plan_incumbency(1, c(polk_precint, story_precint), thresh = 1)

    plans <- redist_smc(iowa_map, 100, constraints = constr, silent = TRUE,
                        sampling_space = "spanning_forest")

    expect_true(
        all(get_plans_matrix(plans)[polk_precint, ] != get_plans_matrix(plans)[story_precint, ])
    )
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
    expect_error(redist_smc(fl_map, 10, sampling_space = "graph_plan",
                            split_params = list(adapt_k_thresh = 1.5)), "0, 1")
    expect_error(redist_smc(fl_map, 0), "positive")
})

test_that("Parallel runs are reproducible", {
    # need to make sure 1 thread and process
    set.seed(5118)
    pl1 <- redist_smc(
        fl_map,
        100,
        runs = 2,
        silent = TRUE,
        ncores = 1L,
        control = list(nproc = 1L)
    )
    set.seed(5118)
    pl2 <- redist_smc(
        fl_map,
        100,
        runs = 2,
        silent = TRUE,
        ncores = 1L,
        control = list(nproc = 1L)
    )

    # runtime related is the only thing that shouldn't be identical
    for (i in 1:2) {
        attr(pl1, "diagnostics")[[i]]$runtime <- NULL
        attr(pl2, "diagnostics")[[i]]$runtime <- NULL
        attr(pl1, "internal_diagnostics")[[i]]$granular_times <- NULL
        attr(pl2, "internal_diagnostics")[[i]]$granular_times <- NULL
        attr(pl1, "internal_diagnostics")[[i]]$time_breakdowns <- NULL
        attr(pl2, "internal_diagnostics")[[i]]$time_breakdowns <- NULL
        attr(pl1, "total_runtime") <- NULL
        attr(pl2, "total_runtime") <- NULL
    }

    expect_identical(pl1, pl2)
})


# Normalizing constant estimation -----------------------------------------
#
# `fl25_enum` holds every 3-district plan of `fl25` along with its population
# deviation, so the exact number of valid plans at a given tolerance is known.
# With `compactness = 0` the SMC target is uniform over those plans, which
# makes that count the true normalizing constant to check estimates against.

exact_log_norm <- function(pop_tol) {
    log(sum(fl25_enum$pop_dev <= pop_tol))
}

test_that("est_norm_unbiased recovers the known normalizing constant", {
    skip_on_cran()
    set.seed(1935)

    map <- suppressMessages(
        redist_map(fl25, ndists = 3, pop_tol = 0.2, total_pop = pop)
    )
    # low resampling efficiency is expected on a map this small and does not
    # affect the estimate
    plans <- suppressWarnings(
        redist_smc(map, 2000, runs = 5, compactness = 0,
                   ncores = 1, silent = TRUE)
    )

    est <- est_norm_unbiased(map, plans)

    # one estimate per run
    expect_length(est, 5)
    expect_true(all(is.finite(est)))

    # the estimator is unbiased on the natural scale, so average there rather
    # than averaging the logs
    pooled <- log(mean(exp(est)))
    expect_equal(pooled, exact_log_norm(0.2), tolerance = 0.05)
})

test_that("est_norm_unbiased tracks the constant across population tolerances", {
    skip_on_cran()

    # A looser tolerance admits strictly more plans, so the normalizing
    # constant must rise with it. This checks the estimator follows a change in
    # the target rather than happening to land on one value.
    pooled <- vapply(c(0.1, 0.2), function(tol) {
        set.seed(1935)
        map <- suppressMessages(
            redist_map(fl25, ndists = 3, pop_tol = tol, total_pop = pop)
        )
        # the tighter tolerance warns about low resampling efficiency, which is
        # expected on a map this small and does not affect the estimate
        plans <- suppressWarnings(
            redist_smc(map, 2000, runs = 5, compactness = 0,
                       ncores = 1, silent = TRUE)
        )
        log(mean(exp(est_norm_unbiased(map, plans))))
    }, numeric(1))

    expect_equal(pooled[1], exact_log_norm(0.1), tolerance = 0.05)
    expect_equal(pooled[2], exact_log_norm(0.2), tolerance = 0.05)
    expect_lt(pooled[1], pooled[2])
})

test_that("est_norm_unbiased accounts for rejected draws", {
    skip_on_cran()
    set.seed(1935)

    map <- suppressMessages(
        redist_map(fl25, ndists = 3, pop_tol = 0.2, total_pop = pop)
    )
    plans <- suppressWarnings(
        redist_smc(map, 1000, runs = 3, compactness = 0,
                   ncores = 1, silent = TRUE)
    )

    diags <- attr(plans, "internal_diagnostics")

    # the rejection bookkeeping the estimator relies on is recorded
    for (d in diags) {
        expect_true("draw_tries_mat" %in% names(d))
        expect_true(is.numeric(d$tries_before_extra_particle))

        # more spanning trees were drawn than particles kept, i.e. draws were
        # actually rejected on the way to a valid plan
        expect_gt(sum(d$draw_tries_mat), 1000 * ncol(d$draw_tries_mat))
    }

    truth <- exact_log_norm(0.2)
    pooled <- log(mean(exp(est_norm_unbiased(map, plans))))
    expect_equal(pooled, truth, tolerance = 0.05)

    # Dropping the attempt counts is not a small correction: without dividing
    # by the number of tries the estimate is off by orders of magnitude. This
    # pins down that the rejection accounting is doing real work.
    ignoring_tries <- vapply(diags, function(d) {
        d$log_blank_map_target_density +
            sum(log(colSums(exp(d$log_incremental_weights_mat))))
    }, numeric(1))

    expect_gt(log(mean(exp(ignoring_tries))) - truth, 5)
})

test_that("est_norm_biased agrees with est_norm_unbiased and needs no extra tries", {
    skip_on_cran()
    set.seed(1935)

    map <- suppressMessages(
        redist_map(fl25, ndists = 3, pop_tol = 0.2, total_pop = pop)
    )
    plans <- suppressWarnings(
        redist_smc(map, 2000, runs = 5, compactness = 0,
                   ncores = 1, silent = TRUE)
    )

    biased <- est_norm_biased(map, plans)
    unbiased <- est_norm_unbiased(map, plans)

    # the biased version leaves the discarded extra particle out of the
    # denominator, so it divides by less and always comes out larger
    expect_true(all(biased > unbiased))

    # that extra particle is one draw among many, so the two barely differ
    expect_equal(log(mean(exp(biased))), exact_log_norm(0.2), tolerance = 0.05)

    # the biased estimator needs no extra bookkeeping, so it still works when
    # the sampler was told not to draw the extra particle
    untracked <- suppressWarnings(
        redist_smc(map, 500, runs = 2, compactness = 0, ncores = 1,
                   control = list(est_norm_unbiased = FALSE), silent = TRUE)
    )
    expect_true(all(is.finite(est_norm_biased(map, untracked))))
    expect_error(
        est_norm_unbiased(map, untracked),
        "tries_before_extra_particle"
    )
})

test_that("est_norm_unbiased handles a hard thresholding constraint", {
    skip_on_cran()
    set.seed(1935)

    map <- suppressMessages(
        redist_map(fl25, ndists = 3, pop_tol = 0.2, total_pop = pop)
    )

    # A 0/1 score with `thresh = 0.5` rejects outright, since any plan scoring
    # at or above the threshold is thrown out at the splitting stage. That
    # restricts the target to plans keeping precincts 1 and 5 together.
    constr <- redist_constr(map) |>
        add_constr_custom_plan(
            1,
            function(plan, seats, num_regions) {
                if (plan[1] == plan[5]) 0 else 1
            },
            thresh = 0.5
        )

    # At the default compactness each plan is weighted by its number of
    # spanning forests, so the normalizing constant is the summed spanning tree
    # count over the enumerated plans that clear both the population bound and
    # the constraint, rather than a plain count of them.
    enum <- fl25_enum$plans
    keep <- fl25_enum$pop_dev <= 0.2 & (enum[1, ] == enum[5, ])
    log_st_kept <- log_st_map(adj, enum[, keep, drop = FALSE], rep(1L, 25), 3L)
    truth <- log(sum(exp(log_st_kept)))

    plans <- suppressWarnings(
        redist_smc(map, 2000, runs = 5, compactness = 1, ncores = 1,
                   sampling_space = "spanning_forest", constraints = constr,
                   silent = TRUE)
    )

    plan_m <- get_plans_matrix(plans)

    # the threshold really did reject: every sampled plan obeys it
    expect_true(all(plan_m[1, ] == plan_m[5, ]))

    # and rejection was a large part of the cost, with several times more trees
    # drawn than particles kept
    tries <- vapply(attr(plans, "internal_diagnostics"), function(d) {
        sum(d$draw_tries_mat) / (2000 * ncol(d$draw_tries_mat))
    }, numeric(1))
    expect_true(all(tries > 3))

    # the estimate tracks the constrained target, not the unconstrained one
    log_st_all <- log_st_map(
        adj, enum[, fl25_enum$pop_dev <= 0.2, drop = FALSE], rep(1L, 25), 3L
    )
    pooled <- log(mean(exp(est_norm_unbiased(map, plans))))

    expect_equal(pooled, truth, tolerance = 0.01)
    expect_lt(pooled, log(sum(exp(log_st_all))))
})

test_that("initial plan weights are normalized and default to equal", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    init_weights_of <- function(init_particles, init_weights = NULL) {
        get_init_plan_params(
            10, attr(iowa_map, "nseats"), iowa_map$pop, attr(iowa_map, "pop_bounds"),
            init_particles, NULL, init_weights
        )$initial_weights
    }

    # no initial plans or weights means all plans start equally weighted
    expect_null(init_weights_of(NULL))

    resampled <- redist_smc(iowa_map, 10, n_steps = 2, silent = TRUE)
    expect_null(init_weights_of(resampled))

    # the weights of plans that were not resampled are used, normalized
    not_resampled <- redist_smc(iowa_map, 10, n_steps = 2, silent = TRUE, resample = FALSE)
    stored_weights <- as.vector(get_plans_weights(not_resampled))
    expect_equal(init_weights_of(not_resampled), stored_weights / sum(stored_weights))

    # passed in weights are used, normalized
    expect_equal(init_weights_of(not_resampled, 1:10), 1:10 / sum(1:10))
})

test_that("initial plan weights are only used to sample the first step's parents", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    init <- redist_smc(iowa_map, 20, n_steps = 2, silent = TRUE)

    # nearly all the weight on plan 3 means every first step parent is plan 3
    restarted <- redist_smc(
        iowa_map, 20, init_particles = init, init_weights = c(1e-12, 1e-12, 1, rep(1e-12, 17)),
        n_steps = 1, silent = TRUE, resample = FALSE
    )
    parents <- attr(restarted, "internal_diagnostics")[[1]]$parent_index_mat[, 1]
    expect_true(all(parents == 3))

    # with seq_alpha the initial weights are not carried into the final weights,
    # which are then just the incremental weights of the one step
    restarted <- redist_smc(
        iowa_map, 20, init_particles = init, init_weights = seq_len(20),
        n_steps = 1, seq_alpha = 0.5, silent = TRUE, resample = FALSE
    )
    log_incremental_weights <- attr(restarted, "internal_diagnostics")[[1]]$log_incremental_weights_mat[, 1]
    log_final_weights <- log(get_plans_weights(subset_sampled(restarted)))
    expect_equal(
        as.vector(log_final_weights - mean(log_final_weights)),
        log_incremental_weights - mean(log_incremental_weights)
    )
})

test_that("mergesplit rounds are sized off the previous mergesplit round", {
    iowa_map <- suppressMessages(redist_map(iowa, ndists = 6, pop_tol = 0.2))
    # mergesplit after every second split, so a round is two steps after the
    # split before it and three after the previous round
    plans <- redist_smc(
        iowa_map, 20, silent = TRUE,
        ms_params = list(frequency = 2, mh_accept_per_smc = 2)
    )
    diag <- attr(plans, "diagnostics")[[1]]
    step_types <- attr(plans, "run_information")[[1]]$step_types
    expect_equal(step_types, c("smc", "smc", "ms", "smc", "smc", "ms", "smc"))

    expected_steps <- function(rate) ceiling(2 * ceiling(1 / ifelse(rate > 0, rate, 0.1)))
    # the first round is sized off the split before it, later ones off the
    # previous mergesplit round
    expect_equal(diag$ms_step_counts[1], expected_steps(diag$accept_rate[2]))
    expect_equal(diag$ms_step_counts[2], expected_steps(diag$accept_rate[3]))
})

test_that("mergesplit runs by default for maps with 10 or more districts", {
    # fewer than 10 districts runs no mergesplit unless asked for
    expect_false(extract_ms_params(list(), 8, ndists = 9)$run_ms)
    small <- extract_ms_params(list(frequency = 1), 4, ndists = 5)
    expect_true(small$run_ms)
    expect_equal(small$mh_accept_per_smc, 1)

    # 10 or more districts defaults to every step with a third of the districts
    default <- extract_ms_params(list(), 9, ndists = 10)
    expect_true(default$run_ms)
    expect_equal(default$mh_accept_per_smc, 3)
    expect_equal(default$ms_frequency, 1)
    expect_equal(default$merge_split_step_vec, rep(c(FALSE, TRUE), 9))

    # the defaults fill in anything not passed
    expect_equal(extract_ms_params(list(frequency = 2), 12, ndists = 13)$mh_accept_per_smc, 4)
    expect_equal(
        extract_ms_params(list(mh_accept_per_smc = 2), 12, ndists = 13)$mh_accept_per_smc, 2
    )

    # NULL turns mergesplit off
    expect_false(extract_ms_params(NULL, 12, ndists = 13)$run_ms)
})

test_that("normalizing constant estimates require seq_alpha = 1", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    plans <- redist_smc(iowa_map, 30, silent = TRUE, seq_alpha = 0.5, ms_params = NULL)
    expect_error(est_norm_unbiased(iowa_map, plans), "seq_alpha = 1")
    expect_error(est_norm_biased(iowa_map, plans), "seq_alpha = 1")
})
