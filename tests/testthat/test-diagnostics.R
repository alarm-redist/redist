# Tests for chain-diagnostic helpers and multi-chain MH acceptance.

test_that("diag_fold and diag_ranknorm behave as documented", {
    set.seed(1)
    x <- c(-2, -1, 0, 1, 2)
    expect_equal(diag_fold(x), abs(x - median(x)))

    rn <- diag_ranknorm(x)
    expect_length(rn, length(x))
    expect_true(all(is.finite(rn)))
    # rank-normalized values are increasing in their input rank
    expect_equal(order(rn), order(x))
})

test_that("diag_calc_rhat ~ 1 for iid chains and > 1 when chains drift apart", {
    set.seed(42)
    n_per <- 200
    n_chain <- 4
    grp <- rep(seq_len(n_chain), each = n_per)

    # iid: rhat ~ 1
    x_iid <- rnorm(n_per * n_chain)
    expect_lt(diag_calc_rhat(x_iid, grp), 1.1)

    # drifted: shift each chain mean far apart -> large rhat
    x_drift <- x_iid + rep(c(0, 5, 10, 15), each = n_per)
    expect_gt(diag_calc_rhat(x_drift, grp), 2)
})

test_that("diag_rhat returns the max over folded and ranknormed Rhat", {
    set.seed(7)
    x <- rnorm(800)
    grp <- rep(1:4, each = 200)
    r <- diag_rhat(x, grp)
    expect_true(is.finite(r))
    expect_gt(r, 0)

    # `split = TRUE` splits each chain in half
    r_split <- diag_rhat(x, grp, split = TRUE)
    expect_true(is.finite(r_split))
})

test_that("redist_mergesplit reports a valid MH acceptance rate", {
    set.seed(5118)
    pl <- redist_mergesplit(fl_map, nsims = 30, warmup = 10, silent = TRUE)
    mh <- get_mh_acceptance_rate(pl)
    expect_true(is.numeric(mh))
    expect_true(is.finite(mh))
    expect_true(mh >= 0 && mh <= 1)
})

test_that("redist_mergesplit with chains > 1 yields one diagnostics entry per chain", {
    skip_on_os("windows")
    set.seed(5118)
    pl <- redist_mergesplit(
        fl_map,
        nsims = 30,
        warmup = 10,
        chains = 2,
        ncores = 2,
        silent = TRUE
    )
    diag <- attr(pl, "diagnostics")
    expect_length(diag, 2)
})

test_that("redist_smc with runs > 1 yields one diagnostics entry per run", {
    set.seed(5118)
    pl <- redist_smc(fl_map, 50, runs = 2, silent = TRUE)
    diag <- attr(pl, "diagnostics")
    expect_length(diag, 2)
})

test_that("summary shows unique plans when stored and computes a missing resample count", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    plans <- subset_sampled(redist_smc(iowa_map, 50, silent = TRUE))
    diagn <- attr(plans, "diagnostics")[[1]]
    run_info <- attr(plans, "run_information")[[1]]
    plans_m <- get_plans_matrix(plans)

    stored <- get_run_nunique_plans(diagn, run_info, TRUE, plans_m)
    expect_equal(stored, diagn$nunique_plans)

    # older runs (without all_step_types) that didn't store the count after
    # resampling have it computed
    no_resample_count <- diagn
    no_resample_count$nunique_plans <- head(diagn$nunique_plans, -1)
    old_run_info <- run_info
    old_run_info$all_step_types <- NULL
    expect_equal(get_run_nunique_plans(no_resample_count, old_run_info, TRUE, plans_m), stored)

    # runs that didn't store any counts don't get a column
    no_counts <- diagn
    no_counts$nunique_plans <- NULL
    expect_null(get_run_nunique_plans(no_counts, run_info, TRUE, plans_m))
    smc_df <- get_smc_summary_df(no_counts, run_info, TRUE, FALSE)$smc_summary_df
    expect_false("unique_plans" %in% names(smc_df))

    smc_df <- get_smc_summary_df(diagn, run_info, TRUE, FALSE, stored)$smc_summary_df
    expect_equal(smc_df$unique_plans, diagn$nunique_plans)
})

test_that("summary marks mergesplit rounds with acceptance below 5%", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    plans <- redist_smc(iowa_map, 50, silent = TRUE, ms_params = list(frequency = 1))
    run_info <- attr(plans, "run_information")[[1]]
    ms_steps <- run_info$step_types == "ms"

    diagn <- attr(plans, "diagnostics")[[1]]
    diagn$accept_rate[ms_steps] <- c(0.5, 0.01, 0.2)
    ms_result <- get_smc_ms_summary_df(diagn, run_info, TRUE)
    expect_true(ms_result$warn_low_ms_accept)
    # only the round below 5% is marked
    marker_col <- ms_result$smc_ms_print_tbl[[ncol(ms_result$smc_ms_print_tbl)]]
    expect_equal(nzchar(trimws(marker_col)), c(FALSE, TRUE, FALSE))

    attr(plans, "diagnostics")[[1]] <- diagn
    expect_match(
        paste(utils::capture.output(summary(plans)), collapse = "\n"),
        "Low merge-split acceptance"
    )

    diagn$accept_rate[ms_steps] <- 0.5
    expect_false(get_smc_ms_summary_df(diagn, run_info, TRUE)$warn_low_ms_accept)
})

test_that("summary flags low SMC acceptance (A) and low effective sample size (W)", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    constr <- redist_constr(iowa_map) |>
        add_constr_status_quo(strength = 1, current = iowa_map$cd_2010) |>
        add_constr_pop_dev(strength = 1, thresh = 1e6)
    expect_equal(
        count_soft_and_hard_constraints(attr(
            redist_smc(iowa_map, 10, constraints = constr, sampling_space = "spanning_forest",
                       silent = TRUE),
            "constraints"
        )),
        c(soft = 1, hard = 1)
    )
    expect_equal(count_soft_and_hard_constraints(list()), c(soft = 0, hard = 0))
    expect_null(count_soft_and_hard_constraints(NA))

    plans <- redist_smc(iowa_map, 50, silent = TRUE)
    run_info <- attr(plans, "run_information")[[1]]
    diagn <- attr(plans, "diagnostics")[[1]]
    diagn$accept_rate[1:3] <- c(0.005, 0.5, 0.5)
    diagn$step_n_eff[1:3] <- c(40, 2, 40)
    diagn$sd_lp[1:3] <- c(0.5, 0.5, 4)

    result <- get_smc_summary_df(diagn, run_info, TRUE, FALSE)
    expect_true(result$warn_low_accept)
    expect_true(result$warn_low_ess)
    flags <- trimws(result$smc_print_tbl[[ncol(result$smc_print_tbl)]])
    # split 1 has low acceptance, split 2 low effective samples, split 3 a
    # large log weight sd
    expect_true(grepl("A", flags[1]) && !grepl("W", flags[1]))
    expect_true(grepl("W", flags[2]) && !grepl("A", flags[2]))
    expect_true(grepl("W", flags[3]))

    attr(plans, "diagnostics")[[1]] <- diagn
    # whitespace is collapsed since cli wraps the text
    out <- gsub("\\s+", " ", paste(utils::capture.output(summary(plans)), collapse = " "))
    expect_match(out, "Low acceptance")
    expect_match(out, "No hard \\(thresholded\\) constraints were used")
    expect_match(out, "Low effective sample size")

    diagn$accept_rate[1:3] <- 0.5
    diagn$step_n_eff[1:3] <- 40
    diagn$sd_lp[1:3] <- 0.5
    diagn$n_eff <- 40
    attr(plans, "diagnostics")[[1]] <- diagn
    # whitespace is collapsed since cli wraps the text
    out <- gsub("\\s+", " ", paste(utils::capture.output(summary(plans)), collapse = " "))
    expect_no_match(out, "Low acceptance")
    expect_no_match(out, "Low effective sample size")
})

test_that("summary flags bottlenecks from the unique parents only", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    plans <- redist_smc(iowa_map, 200, silent = TRUE)
    run_info <- attr(plans, "run_information")[[1]]
    diagn <- attr(plans, "diagnostics")[[1]]
    # 200 plans flags a step with fewer than 80 unique parents
    diagn$unique_survive[1:3] <- c(150, 20, 150)
    # low effective samples alone is W, not a bottleneck
    diagn$step_n_eff[1:3] <- c(150, 150, 5)
    diagn$sd_lp[1:3] <- 0.5

    result <- get_smc_summary_df(diagn, run_info, FALSE, FALSE)
    flags <- trimws(result$smc_print_tbl[[ncol(result$smc_print_tbl)]])
    expect_true(result$warn_bottlenecks)
    expect_equal(grepl("\\*", flags), c(FALSE, TRUE, FALSE))
    expect_equal(grepl("W", flags), c(FALSE, FALSE, TRUE))
})

test_that("summary works for plans that were not resampled", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    plans <- redist_smc(iowa_map, 50, silent = TRUE, resample = FALSE)
    expect_no_error(utils::capture.output(summary(plans)))
    tbl <- get_smc_summary_df(
        attr(plans, "diagnostics")[[1]], attr(plans, "run_information")[[1]], FALSE, FALSE
    )$smc_print_tbl
    expect_equal(rownames(tbl), paste("Split", 1:3))
})

test_that("summary reads unique plans in the order the steps were done", {
    iowa_map <- suppressMessages(redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05))
    plans <- redist_smc(iowa_map, 60, silent = TRUE, ms_params = list(frequency = 1))
    diagn <- attr(plans, "diagnostics")[[1]]
    run_info <- attr(plans, "run_information")[[1]]
    expect_equal(run_info$all_step_types, c("smc", "ms", "smc", "ms", "smc", "resample", "ms"))

    # give each step a distinct count to check which one each row shows
    diagn$nunique_plans <- c(11, 12, 13, 14, 15, 16, 17)
    smc_df <- get_smc_summary_df(diagn, run_info, TRUE, FALSE, diagn$nunique_plans)$smc_summary_df
    expect_equal(smc_df$unique_plans, c(11, 13, 15, 16))
    ms_df <- get_smc_ms_summary_df(diagn, run_info, TRUE, diagn$nunique_plans)$smc_ms_summary_df
    expect_equal(ms_df$unique_plans, c(12, 14, 17))

    # the unique plans after the final mergesplit round are reported separately
    attr(plans, "diagnostics")[[1]] <- diagn
    out <- gsub("\\s+", " ", paste(utils::capture.output(summary(plans)), collapse = " "))
    expect_match(out, "merge-split round there are 17 \\(28%\\) unique plans")

    not_resampled <- redist_smc(
        iowa_map, 60, silent = TRUE, ms_params = list(frequency = 1), resample = FALSE
    )
    out <- paste(utils::capture.output(summary(not_resampled)), collapse = " ")
    expect_no_match(out, "merge-split round there are")

    # older runs without all_step_types resampled after every other step
    old_info <- run_info
    old_info$all_step_types <- NULL
    expect_equal(
        get_all_step_types(old_info, TRUE), c("smc", "ms", "smc", "ms", "smc", "ms", "resample")
    )
    smc_df <- get_smc_summary_df(diagn, old_info, TRUE, FALSE, diagn$nunique_plans)$smc_summary_df
    expect_equal(smc_df$unique_plans, c(11, 13, 15, 17))
})
