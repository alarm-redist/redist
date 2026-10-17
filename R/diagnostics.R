#' Diagnostic information on sampled plans
#'
#' Prints diagnostic information, which varies by algorithm. All algorithms
#' compute the [plans_diversity()] of the samples.
#'
#' For SMC and MCMC, if there are multiple runs/chains, R-hat values will be
#' computed for each summary statistic. These values should be close to 1.
#' If they are not, then there is too much between-chain variation, indicating
#' that there are not enough samples. R-hat values are calculated after
#' rank-normalization and folding.  MCMC chains are split in half before R-hat
#' is computed. For summary statistics that vary across districts, R-hat is
#' calculated for the first district only.
#'
#' For SMC, diagnostics statistics include:
#'
#' * **Effective samples**: the effective sample size at each split, computed
#' using the SMC weights. Larger is better. The percentage in parentheses is the
#' ratio of the effective samples to the total samples.
#' * **Acceptance rate**: the fraction of drawn spanning trees that gave a valid
#' split. Low values make sampling slower but don't by themselves make the
#' sample less accurate.
#' * **Standard deviation of the log weights**: More variable weights (larger s.d.)
#' indicate less efficient sampling.
#' * **Maximum unique plans:** The number of distinct plans from the previous
#' step that became parents of the plans at this step. The percentage in
#' parentheses is relative to the number expected if the weights were even,
#' about 63% of the samples, so it can be above 100%.
#' * **Unique plans**: The number of distinct plans after each split, and after
#' the final resampling. Only shown for plans that stored it.
#' * **Estimated `k` parameter**: For graph space plans displays how many
#' spanning tree edges were considered for cutting at each split. Mostly
#' informational, though large jumps may indicate a need to increase
#' `adapt_k_thresh`.
#' * **Flags**: The right column flags potential problems at each split:
#'     * `*` (bottleneck): the maximum unique plans is below 5% of the samples,
#'       or below the smaller of 40% of the samples and 100. Few plans from the
#'       previous step became parents, which lowers sample diversity and can
#'       increase the variance of estimates across runs. This is often driven
#'       by uneven weights at the previous step or by hard constraints.
#'     * `A` (low acceptance): the acceptance rate is below 1%, usually because
#'       of a tight population tolerance, hard (thresholded) constraints, or
#'       splitting within administrative units.
#'     * `W` (low effective sample size): the effective samples are below 10%
#'       or the standard deviation of the log weights is above 3, usually
#'       because of strong soft constraints or a `compactness` far from 1.
#'
#' For SMC with merge-split steps, a second table shows the acceptance rate of
#' each merge-split step, the number of merge-split steps run per plan, and the
#' number of unique plans afterward. Merge-split steps after a bottleneck can
#' increase the number of unique plans. If the plans were resampled and the run
#' ends with a merge-split step, the resampling is done before that step, so the
#' "Resample" row of the SMC table shows the unique plans before it and the
#' number of unique plans after it is reported separately. Steps with an
#' acceptance rate below 5%
#' are flagged with a dagger (†), since mixing will be slow and increasing
#' `mh_accept_per_smc` will be expensive.
#'
#' For CycleWalk MCMC, additional diagnostic statistics include:
#'
#' * **Cycle walk percentage**: the fraction of iterations that attempted a cycle walk
#' proposal (vs. forest walk). Typically around 10%.
#' * **Success rate**: the fraction of cycle walk attempts that found a valid cycle
#' and cut. Low values (< 30%) may indicate tight constraints.
#' * **Mean acceptance probability**: average MH acceptance probability for successful
#' cycle walk proposals.
#' * **Mean cycle length**: average length of cycles found (typically 3-16).
#' * **Mean valid cuts**: average number of valid cut pairs found per cycle.
#' * **Failures**: breakdown of why cycle walk proposals failed (no adjacent districts,
#' too few boundary edges, no path found, or no valid cuts).
#'
#' In the event of problematic diagnostics, the function will provide
#' suggestions for improvement.
#'
#' @param object a [redist_plans] object
#' @param district For R-hat values, which district to use for district-level
#' summary statistics. We strongly recommend calling `match_numbers()` or
#' `number_by()` before examining these district-level statistics. A value of
#' `FALSE` indicates to compute the rhats for all districts
#' @param all_runs When there are multiple SMC runs, whether to show detailed
#' summary statistics for all runs, or only the first run (the default).
#' @param vi_max The maximum number of plans to sample in computing the pairwise
#' variation of information distance (sample diversity).
#' @param order_stats Whether or not to compute rhats on the ordered district
#' statistics.
#' @param rhat_thresh What values to use when checking convergence. It should
#' consist of two named values: `q99` and `max`, where `q99` is the tighter of
#' the two. Any name that is omitted falls back to its default.
#' Strong convergence is when all rhats are less than or equal to the `q99`
#' value, and weak convergence is when the 99th quantile is less than or equal
#' to `q99` and all rhats are less than or equal to `max`. Chains that meet
#' neither are reported as not converged.
#'
#' @param ncores Number of cores to use when computing R-hat values. Defaults
#' to 1, which runs serially. Raising this is worthwhile for large ensembles,
#' where the R-hat computation dominates the cost of this function.
#'
#' @param \dots additional arguments (ignored)
#'
#' @return A list of diagnostic information, invisibly. For SMC this has a
#' `smc_diagnostic_dfs` data frame with the statistics for each split, and a
#' `ms_diagnostic_dfs` data frame for any merge-split steps. For MCMC it has a
#' `mcmc_diagnostic_dfs` data frame. If R-hat values were computed, it also has
#' `rhats_df` and `rhat`.
#'
#' @examples
#' data(iowa)
#' iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.1)
#' plans <- redist_smc(iowa_map, 100)
#' summary(plans)
#'
#' @method summary redist_plans
#' @concept analyze
#' @md
#' @export
summary.redist_plans <- function(
    object,
    district = FALSE,
    all_runs = FALSE,
    vi_max = 100,
    order_stats = TRUE,
    rhat_thresh = getOption("redist.rhat_thresh", c(q99=1.05, max=1.1)),
    ncores = 1,
    ...
) {
    cli::cli_process_done(done_class = "") # in case an earlier

    algo <- attr(object, "algorithm")
    name <- deparse(substitute(object))

    object <- subset_sampled(object)
    all_diagn <- attr(object, "diagnostics")
    all_run_info <- attr(object, "run_information")
    plans_m <- get_plans_matrix(object)
    n_samp <- ncol(plans_m)
    n_distr <- attr(object, "ndists")
    resampled <- attr(object, "resampled")
    if (is.null(n_distr)) {
        n_distr <- max(plans_m[, 1])
    }

    out_list <- list()

    if (n_distr == 1 || nrow(plans_m) == 1) {
        cat_cli("{fmt_comma(n_samp)}{cli::qty(n_samp)} sampled plan{?s} of
                 {n_distr} district{?s} on
                 {fmt_comma(nrow(plans_m))}{cli::qty(nrow(plans_m))} unit{?s}")
        return(invisible(out_list))
    }

    prec_pop <- attr(object, "prec_pop")
    if (is.null(prec_pop)) {
        cli::cli_warn(c("Precinct population missing; plan diversity estimates may be misleading.",
                   ">" = 'Run `attr({name}, "prec_pop") <- <map object>$<pop column>` to fix.'))
        prec_pop <- rep(1, nrow(plans_m))
    }
    est_div <- plans_diversity(object, total_pop = prec_pop, n_max = vi_max)
    div_rg <- format(quantile(est_div, c(0.1, 0.9)), digits = 2)
    div_bad <- (mean(est_div) <= 0.35) || (mean(est_div <= 0.25) > 0.1)

    summary_supported_algs <- c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE, MCMC_ALG_TYPE, "flip", "cyclewalk")

    # ignore if not a supported algorithm
    if (!algo %in% summary_supported_algs) {
        cli::cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
    }

    # bool for not flip
    revamped_alg <- algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE, MCMC_ALG_TYPE)

    # checks if plans are earlier than version 5.0
    pre_v5_plans <- is.null(attr(object, "version")) || attr(object, "version") < "5"

    # get the display name
    alg_display_name <- dplyr::case_when(
        algo == SMC_ALG_TYPE ~ "SMC",
        algo == MS_SMC_ALG_TYPE ~ "SMC with Merge-Split MCMC Steps",
        algo == MCMC_ALG_TYPE ~ "Merge-Split MCMC",
        algo == "flip" ~ "Flip MCMC",
        algo == "cyclewalk" ~ "CycleWalk MCMC"
    )
    cat_cli("{.strong {alg_display_name}:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")

    if (pre_v5_plans) {
        display_sampling_space <- GRAPH_PLAN_SPACE_SAMPLING
        display_splitting_method <- NAIVE_K_SPLITTING
        cat_cli("Plans sampled on {display_sampling_space} using the {display_splitting_method} forward kernel.")
        cat_cli("Forward kernel parameters: {.arg adapt_k_thresh}={format(all_diagn[[1]]$adapt_k_thresh, digits=3)}")
    } else if (revamped_alg && !pre_v5_plans) {
        # if revamped alg check that sampling space and splitting methods all the same
        # check same sampling space
        all_sampling_spaces <- sapply(all_run_info, function(x) x$sampling_space)
        if (length(unique(all_sampling_spaces)) != 1) {
            cli::cli_abort("{.fn summary} is not supported for plans sampled using different sampling spaces.")
        }
        sampling_space <- all_sampling_spaces[1]
        display_sampling_space <- dplyr::case_when(
            sampling_space == GRAPH_PLAN_SPACE_SAMPLING ~ "Graph Space",
            sampling_space == FOREST_SPACE_SAMPLING ~ "Spanning Forest Space",
            sampling_space == LINKING_EDGE_SPACE_SAMPLING ~ "Linking Edge Forest Space"
        )
        # check same splitting method
        all_splitting_methods <- sapply(all_run_info, function(x) x$split_method)
        if (length(unique(all_splitting_methods)) != 1) {
            cli::cli_abort("{.fn summary} is not supported for plans sampled using different splitting methods")
        }
        split_method <- all_splitting_methods[1]
        display_splitting_method <- dplyr::case_when(
            split_method == NAIVE_K_SPLITTING ~ "Naive Top K",
            split_method == UNIF_VALID_EDGE_SPLITTING ~ "Uniform Valid Edge",
            split_method == EXP_BIGGER_ABS_DEV_SPLITTING ~ "Exponential Absolute Deviance"
        )

        # check the splitting parameters are all the same
        all_forward_kernel_params <- lapply(all_diagn, function(x) x$forward_kernel_params)
        for (i in seq_len(length(all_forward_kernel_params))) {
            # don't want to compare cut k used since won't be the same
            all_forward_kernel_params[[i]]$cut_k_used <- NULL
            if (!identical(all_forward_kernel_params[[1]], all_forward_kernel_params[[i]])) {
                cli::cli_abort("{.fn summary} is not supported for plans sampled using different splitting parameters")
            }
        }
        forward_kernel_params <- all_forward_kernel_params[[1]]

        cat_cli("Plans sampled on {display_sampling_space} using the {display_splitting_method} forward kernel.")

        if (split_method == NAIVE_K_SPLITTING) {
            # only display adapt k threshold if k values estimated
            if (forward_kernel_params$estimate_cut_k) {
                cat_cli("Forward kernel parameters: {.arg adapt_k_thresh}={format(forward_kernel_params$adapt_k_thresh, digits=3)}")
            } else {
                cat_cli("Forward kernel parameters: {.arg manual_k_params}={forward_kernel_params$manual_k_params}")
            }
        } else if (split_method == UNIF_VALID_EDGE_SPLITTING) {} else if (
            split_method == EXP_BIGGER_ABS_DEV_SPLITTING
        ) {
            cat_cli("Forward Kernel Parameters: {.arg splitting_alpha}={format(forward_kernel_params$splitting_alpha, digits=2)}")
        }
    }

    # print algorithm specific parameters
    print_algo_specific_params(algo, all_diagn, all_run_info, pre_v5_plans)

    cat_cli("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
    if (div_bad) {
        cat_cli_fmt(cli::cli_alert_danger("{.strong WARNING:} Low plan diversity"))
    }

    # now compute rhats if more than 1 chain
    rhat_cols <- select_rhat_columns(object, district, n_distr, n_samp)

    warn_converge <- FALSE
    # do nothing if no additional columns or no chain column
    rhats_computed <- length(rhat_cols) > 0

    if (rhats_computed) {
        district <- validate_rhat_district(district, n_distr)

        split_rhat <- algo %in% c(MCMC_ALG_TYPE, "flip")

        # get rhats
        rhats_df <- compute_all_rhats(
            # drop everything but the columns to save size
            as.data.frame(object)[c("chain", "district", rhat_cols)],
            rhat_cols,
            order_stats,
            district,
            n_distr,
            split_rhat,
            ncores
        )

        # remove NA rhats
        rhats_df <- drop_na_rhats(rhats_df)

        out_list[["rhats_df"]] <- rhats_df
        out_list[["rhat"]] <- rhats_df$rhat

        thresh <- resolve_rhat_thresh(rhat_thresh)
        convergence <- classify_rhat_convergence(rhats_df$rhat, thresh)

        print_rhat_summary(rhats_df, convergence, thresh, order_stats)

        warn_converge <- convergence$status == "not_converged"
    }

    # Now print algorithm specific diagnostics
    if (algo == "smc" && pre_v5_plans) {
        out <- legacy_print_smc_information(
            name,
            all_runs,
            object,
            algo,
            div_rg,
            all_diagn,
            warn_converge,
            div_bad
        )
    } else if (algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE)) {
        smc_run_dfs <- list()
        smc_ms_run_dfs <- list()
        n_runs <- length(all_diagn)
        warn_bottlenecks <- FALSE
        warn_low_accept <- FALSE
        warn_low_ess <- FALSE
        warn_low_ms_accept <- FALSE

        for (i in seq_len(n_runs)) {
            # each run's plans are its block of columns in run order
            run_nsims <- all_run_info[[i]]$nsims
            run_plans_m <- plans_m[, (i - 1) * run_nsims + seq_len(run_nsims), drop = FALSE]
            run_nunique_plans <- get_run_nunique_plans(
                all_diagn[[i]], all_run_info[[i]], resampled, run_plans_m
            )

            smc_summary_result_list <- get_smc_summary_df(
                all_diagn[[i]],
                all_run_info[[i]],
                resampled,
                warn_bottlenecks,
                run_nunique_plans
            )

            smc_run_dfs[[i]] <- smc_summary_result_list$smc_summary_df
            smc_tbl_print <- smc_summary_result_list$smc_print_tbl
            warn_bottlenecks <- smc_summary_result_list$warn_bottlenecks
            warn_low_accept <- warn_low_accept || smc_summary_result_list$warn_low_accept
            warn_low_ess <- warn_low_ess || smc_summary_result_list$warn_low_ess

            if (i == 1 || isTRUE(all_runs)) {
                cat_cli("Sampling diagnostics for SMC run {i} of {n_runs} ({fmt_comma(run_nsims)} samples)")
                print(smc_tbl_print, digits = 2)
                cat("\n")
            }

            if (algo == MS_SMC_ALG_TYPE) {
                smc_ms_summary_result_list <- get_smc_ms_summary_df(
                    all_diagn[[i]],
                    all_run_info[[i]],
                    resampled,
                    run_nunique_plans
                )
                smc_ms_run_dfs[[i]] <- smc_ms_summary_result_list$smc_ms_summary_df
                smc_ms_tbl_print <- smc_ms_summary_result_list$smc_ms_print_tbl
                warn_low_ms_accept <- warn_low_ms_accept ||
                    smc_ms_summary_result_list$warn_low_ms_accept

                if (i == 1 || isTRUE(all_runs)) {
                    # If the plans were resampled before the final mergesplit
                    # round then the Resample row above is from before it, so
                    # also report the unique plans after it
                    run_all_step_types <- get_all_step_types(all_run_info[[i]], resampled)
                    if (resampled && !is.null(run_nunique_plans) &&
                        utils::tail(run_all_step_types, 1) == "ms") {
                        final_unique <- utils::tail(run_nunique_plans, 1)
                        final_unique_pct <- sprintf("%0.0f%%", 100 * final_unique / run_nsims)
                        cat_cli("After the post-resampling merge-split round there are {fmt_comma(final_unique)} ({final_unique_pct}) unique plans.")
                        cat("\n")
                    }
                    cat_cli("Sampling diagnostics for merge-split steps of SMC run {i} of {n_runs} ({fmt_comma(run_nsims)} samples)")
                    print(smc_ms_tbl_print, digits = 2)
                    cat("\n")
                }
            }
        }

        out_list[["smc_diagnostic_dfs"]] <- bind_rows(smc_run_dfs)
        if (algo == MS_SMC_ALG_TYPE) {
            out_list[["ms_diagnostic_dfs"]] <- bind_rows(smc_ms_run_dfs)
        }

        cat_cli_fmt(cli::cli_li(cli::col_grey(
            "
            Watch out for low effective samples, very low acceptance rates (less than 1%),
            large std. devs. of the log weights (more than 3 or so),
            and low numbers of unique plans.
            R-hat values for summary statistics should be between 1 and 1.05."
        )))

        if (div_bad) {
            cat_cli_fmt(cli::cli_li(
                "{.strong Low diversity:} Check for potential bottlenecks.
                        Increase the number of samples or `mh_accept_per_smc`.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance."
            ))
        }
        if (warn_converge) {
            cat_cli_fmt(cli::cli_li(
                "{.strong SMC convergence:} Increase the number of samples or `mh_accept_per_smc`.
                        If you are experiencing low plan diversity or bottlenecks as well,
                        address those issues first."
            ))
        }
        if (warn_bottlenecks) {
            bottleneck_causes <- "This is often driven by uneven weights at the
                previous step (see W) or by hard constraints that only allow a
                subset of plans to be split."
            bottleneck_ms <- if (algo == MS_SMC_ALG_TYPE) {
                "Merge-split steps after a bottleneck can alleviate this somewhat by
                increasing diversity, though the new plans stay similar to their
                shared ancestors. Check how much the number of unique plans increases
                in the merge-split step after each flagged split."
            } else {
                "Running merge-split steps after each split (see `ms_params`) can
                alleviate this somewhat by increasing diversity, though the new plans
                stay similar to their shared ancestors. The number of unique plans
                after each merge-split step shows how much they help."
            }
            cat_cli_fmt(cli::cli_li(
                "(*) {.strong Bottlenecks found:} Only a small percentage of plans
                from the previous step were successfully split, potentially degrading
                sample quality. {bottleneck_causes} {bottleneck_ms}
                Consider weakening or removing
                            constraints, or increasing the population tolerance.
                            If the acceptance rate drops quickly in the final splits,
                            try increasing {.arg pop_temper} by 0.01.
                            If the weight variance (Log wgt. sd) increases steadily
                            or is particularly large for the \"Resample\" step,
                            consider increasing {.arg seq_alpha} if merge-split steps
                            are not also being used."
            ))
        }
        constraint_counts <- count_soft_and_hard_constraints(attr(object, "constraints"))
        # "There is 1 ... constraint" / "There are n ... constraints"
        describe_constraint_count <- function(n, type) {
            sprintf(
                "There %s %d %s constraint%s that could be contributing.",
                if (n == 1) "is" else "are", n, type, if (n == 1) "" else "s"
            )
        }
        if (warn_low_accept) {
            accept_causes <- if (is.null(constraint_counts)) {
                "This is usually caused by a tight population tolerance or by hard
                (thresholded) constraints rejecting splits."
            } else if (constraint_counts[["hard"]] > 0) {
                paste(
                    "This is usually caused by a tight population tolerance or by hard
                    (thresholded) constraints rejecting splits.",
                    describe_constraint_count(constraint_counts[["hard"]], "hard (thresholded)")
                )
            } else {
                "No hard (thresholded) constraints were used, so the population
                tolerance or the geography of the map is the likely cause."
            }
            if (isTRUE(attr(object, "num_admin_units") > 1)) {
                accept_causes <- paste(
                    accept_causes,
                    "Splitting within administrative units can also lower the acceptance rate."
                )
            }
            # only suggest loosening thresholds if there may be hard constraints
            accept_fixes <- if (isTRUE(constraint_counts[["hard"]] == 0)) {
                "Consider increasing the population tolerance."
            } else {
                "Consider increasing the population tolerance or loosening constraint
                thresholds."
            }
            low_accept_pct <- sprintf("%g%%", 100 * SMC_LOW_ACCEPT_RATE)
            cat_cli_fmt(cli::cli_li(
                "(A) {.strong Low acceptance:} Fewer than {low_accept_pct} of tree
                        draws gave a valid split. This makes sampling slower but
                        doesn't by itself make the sample less accurate.
                        {accept_causes} {accept_fixes} If the acceptance rate only
                        drops sharply in the final splits, try increasing
                        {.arg pop_temper} by 0.01."
            ))
        }
        if (warn_low_ess) {
            ess_causes <- if (is.null(constraint_counts)) {
                ""
            } else if (constraint_counts[["soft"]] > 0) {
                describe_constraint_count(constraint_counts[["soft"]], "soft")
            } else {
                "No soft constraints were used."
            }
            compactness <- attr(object, "compactness")
            if (is.numeric(compactness) && isTRUE(compactness != 1)) {
                ess_causes <- paste(
                    ess_causes,
                    sprintf("These plans were sampled with `compactness` = %g.", compactness)
                )
            }
            pop_temper <- all_diagn[[1]]$pop_temper
            if (isTRUE(pop_temper > 0)) {
                ess_causes <- paste(
                    ess_causes,
                    sprintf("Population tempering (`pop_temper` = %g) trades acceptance
                            in the final splits for more uneven weights.", pop_temper)
                )
            }
            if (any(vapply(all_diagn, function(d) isTRUE(d$seq_alpha < 1), logical(1)))) {
                ess_causes <- paste(
                    ess_causes,
                    "If the weight variance grows steadily across the splits,
                    increasing `seq_alpha` may help."
                )
            }
            low_eff_pct <- sprintf("%g%%", 100 * SMC_LOW_EFF_SAMPLE_FRAC)
            cat_cli_fmt(cli::cli_li(
                "(W) {.strong Low effective sample size:} The weights are very uneven
                        (effective samples below {low_eff_pct} or a log weight sd above
                        {SMC_HIGH_LOG_WGT_SD}), which is usually caused by strong soft
                        constraints or a {.arg compactness} far from 1.
                        {ess_causes} Consider weakening the strongest
                        constraints."
            ))
        }
        if (warn_low_ms_accept) {
            cat_cli_fmt(cli::cli_li(
                "(\u2020) {.strong Low merge-split acceptance:} An acceptance rate
                        below 5% means mixing will be slow, and increasing the
                        expected number of successful steps
                        ({.arg mh_accept_per_smc}) will be expensive."
            ))
        }
    } else if (algo == "cyclewalk") {
        mh <- attr(object, "mh_acceptance")
        if (!is.null(mh)) {
            accept_rate <- sprintf("%0.1f%%", 100 * mh)
            cat_cli("Chain acceptance rate{?s}: {accept_rate}")
        }

        out <- tibble(
            accept_rate = mh %||% NA_real_,
            div_q10 = div_rg[1],
            div_q90 = div_rg[2]
        )

        out_list[["mcmc_diagnostic_dfs"]] <- out
    } else if (algo %in% c("mergesplit", "flip")) {
        accept_rate <- sprintf("%0.1f%%", 100 * attr(object, "mh_acceptance"))
        cat_cli("Chain acceptance rate{?s}: {accept_rate}")

        out <- tibble(
            accept_rate = attr(object, "mh_acceptance"),
            div_q10 = div_rg[1],
            div_q90 = div_rg[2]
        )

        out_list[["mcmc_diagnostic_dfs"]] <- out

        cat_cli_fmt(cli::cli_li(cli::col_grey(
            "
            Watch out for low acceptance rates (less than 10%).
            R-hat values for summary statistics should be between 1 and 1.05.
            For district-level statistics (like district partisan leans), you
            should call `match_numbers()` or `number_by()` before examining
            the R-hat values."
        )))

        if (div_bad) {
            cat_cli_fmt(cli::cli_li(
                "{.strong Low diversity:} Increase the number of samples.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance."
            ))
        }
        if (warn_converge) {
            cat_cli_fmt(cli::cli_li(
                "{.strong Chain convergence:} Increase the number of samples.
                        If you are experiencing low plan diversity, address that issue first."
            ))
        }
    } else {
        cli::cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
    }

    invisible(out_list)
}


fmt_comma <- function(x) {
    format(x, nsmall = 0, digits = 1, big.mark = ",")
}


#' Pretty prints relevant parameters for an algorithm type
#'
#' Given a dataframe with `chain` and `district` columns and all other columns
#' being statistics of interest this computes the rhats for each of the districts.
#' If `order_stats` is true then computes the rhats on the order statistics for
#' each plan ie `district == 1` represents the smallest value for a plan.
#'
#' @param algo The algorithm type
#' @param all_diagn List of `diagnostics` for each run of the plan.
#' @param all_run_info The list of `run_information` outputs for each run of
#' the plan.
#' @param pre_v5_plans Boolean for whether or not the plans are legacy plans
#' (ie pre Redist 5.0 meaning it was only graph space, simple weights)
#'
#' @returns A long dataframe of rhats
#' @noRd
print_algo_specific_params <- function(algo, all_diagn, all_run_info, pre_v5_plans) {
    if (algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE)) {
        if (!pre_v5_plans) {
            weight_type <- all_run_info[[1]]$weight_type
        } else {
            weight_type <- "simple"
        }
        cat_cli_fmt(cli::cli_bullets(c(
            "SMC Parameters:",
            "*"="{.arg weight_type} = {weight_type}",
            "*"="{.arg seq_alpha} = {format(all_diagn[[1]]$seq_alpha, digits=2)}",
            "*"="{.arg pop_temper} = {format(all_diagn[[1]]$pop_temper, digits=3)}"
        )))
    } else if (algo == MCMC_ALG_TYPE) {
        cat_cli_fmt(cli::cli_bullets(c(
            "MCMC Parameters:",
            "*"="{.arg warmup} = {format(all_diagn[[1]]$warmup)}",
            "*"="{.arg thin} = {format(all_diagn[[1]]$thin)}",
            "*"="{.arg total_steps} = {format(all_diagn[[1]]$total_steps)}"
        )))
    }
    if (algo == MS_SMC_ALG_TYPE) {
        total_ms_steps <- sum(all_run_info[[1]]$step_types == "ms")

        cat_cli_fmt(cli::cli_bullets(c(
            "Mergesplit Parameters:",
            "*"="{.arg total_ms_steps} = {format(total_ms_steps, digits=2)}",
            "*"="{.arg mh_accept_per_smc} = {format(all_run_info[[1]]$mh_accept_per_smc, digits=2)}",
            "*"="{.arg pair_rule} = {format(all_run_info[[1]]$pair_rule, digits=2)}"
        )))
    }
}


#' Get the number of unique plans after each step of an SMC run
#'
#' Not every run stores these, and some runs only stored them for the
#' splitting and mergesplit steps and not the final resampling step. In that
#' case the count for the resampling step is computed from the run's final
#' plans.
#'
#' @param diagn The run's diagnostics
#' @param run_info The run's `run_information`
#' @param resampled Whether the plans were resampled after the last step
#' @param run_plans_m The run's final plans matrix
#'
#' @returns A vector with the number of unique plans after each step, and after
#' the resampling step last if resampled, or `NULL` if they weren't stored.
#' @noRd
get_run_nunique_plans <- function(diagn, run_info, resampled, run_plans_m) {
    nunique_plans <- diagn$nunique_plans
    n_steps <- length(run_info$step_types)
    if (is.null(nunique_plans) || n_steps == 0) {
        return(NULL)
    }

    # Compute the count after resampling if it wasn't stored. Only older runs
    # without `all_step_types` didn't store it, and they resampled after every
    # other step so their final plans are the resampled plans.
    if (resampled && length(nunique_plans) == n_steps &&
        is.null(run_info$all_step_types)) {
        num_regions <- dplyr::n_distinct(run_plans_m[, 1])
        nunique_plans <- c(nunique_plans, nrow(get_plan_counts(run_plans_m, num_regions)))
    }

    if (length(nunique_plans) != n_steps + resampled) {
        return(NULL)
    }
    nunique_plans
}

# SMC steps are flagged in `summary()` with "A" if the acceptance rate is below
# SMC_LOW_ACCEPT_RATE and with "W" if the effective sample fraction is below
# SMC_LOW_EFF_SAMPLE_FRAC or the log weight sd is above SMC_HIGH_LOG_WGT_SD
SMC_LOW_ACCEPT_RATE <- 0.01
SMC_LOW_EFF_SAMPLE_FRAC <- 0.1
SMC_HIGH_LOG_WGT_SD <- 3

#' Counts the soft and hard (thresholded) constraints plans were sampled with
#'
#' @param constraints The `constraints` attribute of a `redist_plans` object
#'
#' @returns A vector with the number of `soft` and `hard` constraints, or
#' `NULL` if the constraints aren't known
#' @noRd
count_soft_and_hard_constraints <- function(constraints) {
    if (!is.list(constraints)) {
        return(NULL)
    }
    # as.logical so no constraints gives logical(0) rather than NULL
    is_hard <- as.logical(unlist(lapply(constraints, function(constraint_type) {
        vapply(constraint_type, function(x) isTRUE(x$hard_constraint), logical(1))
    })))
    c(soft = sum(!is_hard), hard = sum(is_hard))
}

#' Formats unique plan counts as the count and its percent of the samples
#' @noRd
format_unique_plans <- function(unique_plans, n_samp) {
    str_glue("{fmt_comma(unique_plans)} ({sprintf('%3.0f%%', 100 * unique_plans / n_samp)})")
}

#' Get Summary and Print Dataframe for SMC steps
#'
#' @noRd
get_smc_summary_df <- function(diagn, run_info, resampled, warn_bottlenecks,
                               nunique_plans = NULL) {
    n_samp <- run_info$nsims
    run_sampling_space <- run_info$sampling_space

    smc_accept_rate <- diagn$accept_rate[run_info$step_types == "smc"]

    smc_step_mask <- run_info$step_types == "smc"

    run_summary_df <- tibble(
        n_eff = diagn$step_n_eff,
        eff = diagn$step_n_eff / n_samp,
        accept_rate = smc_accept_rate,
        sd_log_wgt = diagn$sd_lp,
        max_unique = diagn$unique_survive[seq_len(length(smc_accept_rate))]
    )

    run_summary_names <- c("Eff. samples (%)", "Acc. rate",
                           "Log wgt. sd", " Max. unique")

    # if graph space then add the k
    if (run_sampling_space == GRAPH_PLAN_SPACE_SAMPLING) {
        run_summary_df$est_k <- diagn$forward_kernel_params$cut_k_used[smc_step_mask]
        if (diagn$forward_kernel_params$estimate_cut_k) {
            run_summary_names <- c(run_summary_names, "Est. k")
        } else {
            run_summary_names <- c(run_summary_names, "Manual k")
        }
    }

    # add the number of unique plans after each split if stored
    # the steps in the order they were done, which `nunique_plans` follows
    all_step_types <- get_all_step_types(run_info, resampled)
    if (!is.null(nunique_plans)) {
        run_summary_df$unique_plans <- nunique_plans[all_step_types == "smc"]
        run_summary_names <- c(run_summary_names, "Unique plans")
    }

    # add row if resampled
    if (resampled) {
        run_summary_df[nrow(run_summary_df) + 1, ] <- NA
        # add eff sample size and how many survived resample
        run_summary_df[nrow(run_summary_df), "n_eff"] <- diagn$n_eff
        run_summary_df[nrow(run_summary_df), "eff"] <- diagn$n_eff / n_samp
        run_summary_df[nrow(run_summary_df), "max_unique"] <- tail(diagn$unique_survive, 1)
        if (!is.null(nunique_plans)) {
            # the resampling may be before the final mergesplit round
            run_summary_df[nrow(run_summary_df), "unique_plans"] <-
                nunique_plans[all_step_types == "resample"]
        }
    }

    # add extra column name for the flags (*, A, W)
    run_summary_names <- c(run_summary_names, "")

    tbl_print <- as.data.frame(run_summary_df)
    min_n <- max(0.05 * n_samp, min(0.4 * n_samp, 100))
    # * bottleneck: few distinct plans from the previous step became parents of
    # this step's plans. This only uses the unique parents since uneven weights
    # are flagged separately with W
    bottlenecks <- dplyr::coalesce(tbl_print$max_unique < min_n, FALSE)
    warn_bottlenecks <- warn_bottlenecks || any(bottlenecks)
    # A low acceptance: few tree draws give a valid split, usually from the
    # population tolerance or hard constraints
    low_accept <- dplyr::coalesce(tbl_print$accept_rate < SMC_LOW_ACCEPT_RATE, FALSE)
    # W low effective sample size: very uneven weights, usually from strong
    # soft constraints
    low_ess <- dplyr::coalesce(
        tbl_print$eff < SMC_LOW_EFF_SAMPLE_FRAC | tbl_print$sd_log_wgt > SMC_HIGH_LOG_WGT_SD,
        FALSE
    )
    tbl_print$flags <- paste0(
        " ", ifelse(bottlenecks, "*", ""), ifelse(low_accept, "A", ""),
        ifelse(low_ess, "W", ""), " "
    )
    tbl_print$n_eff <- with(
        tbl_print,
        str_glue("{fmt_comma(n_eff)} ({sprintf('%0.1f%%', 100*eff)})")
    )
    tbl_print$eff <- NULL
    tbl_print$accept_rate <- with(tbl_print, sprintf("%0.1f%%", 100 * accept_rate))
    max_pct <- with(tbl_print, max_unique / (-n_samp * expm1(-1)))
    tbl_print$max_unique <- with(
        tbl_print,
        str_glue("{fmt_comma(max_unique)} ({sprintf('%3.0f%%', 100*max_pct)})")
    )
    if (!is.null(nunique_plans)) {
        tbl_print$unique_plans <- format_unique_plans(tbl_print$unique_plans, n_samp)
    }

    names(tbl_print) <- run_summary_names
    # the last row is the resampling step only if the plans were resampled
    new_row_names <- paste("Split", seq_len(nrow(tbl_print) - resampled))

    if (resampled) {
        new_row_names <- c(new_row_names, "Resample")
    }

    rownames(tbl_print) <- new_row_names

    list(
        smc_summary_df = run_summary_df,
        smc_print_tbl = tbl_print,
        warn_bottlenecks = warn_bottlenecks,
        warn_low_accept = any(low_accept),
        warn_low_ess = any(low_ess)
    )
}


#' Get Summary and Print Dataframe for Mergesplit within SMC steps
#'
#' @noRd
get_smc_ms_summary_df <- function(diagn, run_info, resampled, nunique_plans = NULL) {
    n_samp <- run_info$nsims
    run_sampling_space <- run_info$sampling_space

    ms_step_mask <- run_info$step_types == "ms"
    ms_accept_rate <- diagn$accept_rate[ms_step_mask]
    ms_steps_per_plan <- diagn$ms_step_counts

    run_summary_df <- tibble(
        accept_rate = ms_accept_rate,
        ms_steps = ms_steps_per_plan
    )

    run_summary_names <- c("Acc. rate", "MS Steps per Plan")

    # add the number of unique plans after each mergesplit round if stored
    if (!is.null(nunique_plans)) {
        run_summary_df$unique_plans <-
            nunique_plans[get_all_step_types(run_info, resampled) == "ms"]
        run_summary_names <- c(run_summary_names, "Unique plans")
    }

    tbl_print <- as.data.frame(run_summary_df)

    # mark rounds with low acceptance, like the bottleneck asterisk in the SMC table
    low_accept <- dplyr::coalesce(tbl_print$accept_rate < 0.05, FALSE)
    tbl_print$low_accept <- ifelse(low_accept, " \u2020 ", "")
    run_summary_names <- c(run_summary_names, "")

    tbl_print$accept_rate <- with(tbl_print, sprintf("%0.1f%%", 100 * accept_rate))
    if (!is.null(nunique_plans)) {
        tbl_print$unique_plans <- format_unique_plans(tbl_print$unique_plans, n_samp)
    }

    names(tbl_print) <- run_summary_names
    new_row_names <- paste("MS Step", seq_len(nrow(tbl_print)))

    rownames(tbl_print) <- new_row_names

    list(
        smc_ms_summary_df = run_summary_df,
        smc_ms_print_tbl = tbl_print,
        warn_low_ms_accept = any(low_accept)
    )
}

diag_fold <- function(x) {
    abs(x - median(x))
}

diag_ranknorm <- function(x) {
    qnorm(rank(x) / (length(x) + 1))
}

#' Precompute the chain grouping that the rhat calculation needs
#'
#' The grouping depends only on `grp`, never on the statistic being
#' summarised, so it can be built once and reused across every statistic and
#' across both halves of `diag_rhat()`. Splitting the row indices once is much
#' cheaper than letting `tapply()` and `table()` rebuild the factor on every
#' call.
#'
#' @param grp vector of chain ids, with each chain in a contiguous run
#' @param split whether to split each chain in half first
#'
#' @returns A list with `idx`, a list of index vectors in factor level order,
#' and `n`, the mean group size.
#' @noRd
make_rhat_groups <- function(grp, split = FALSE) {
    if (split) {
        lengths <- rle(grp)$lengths
        grp <- grp + do.call(c, lapply(lengths, function(l) rep(c(0.0, 0.5), each = l / 2)))
    }

    idx <- split(seq_along(grp), grp)

    list(idx = idx, n = mean(lengths(idx)))
}

diag_calc_rhat <- function(x, grp = NULL, groups = make_rhat_groups(grp)) {
    n <- groups$n
    group_means <- vapply(groups$idx, function(rows) mean(x[rows]), numeric(1))
    group_vars <- vapply(groups$idx, function(rows) var(x[rows]), numeric(1))

    var_between <- n * var(group_means)
    var_within <- mean(group_vars)
    sqrt((var_between / var_within + n - 1) / n)
}

diag_rhat <- function(x, grp = NULL, split = FALSE,
                      groups = make_rhat_groups(grp, split)) {
    max(
        diag_calc_rhat(diag_ranknorm(x), groups = groups),
        diag_calc_rhat(diag_ranknorm(diag_fold(x)), groups = groups)
    )
}


#' Compute an rhat for every (statistic, district) pair
#'
#' Flattens the per-district row indices and group memberships into the
#' CSR-style layout `compute_rhats_cpp()` expects, then hands the numeric work
#' to C++, which parallelizes over the pairs.
#'
#' The groups are built in R and passed down so that the split-chain
#' convention is defined in one place only.
#'
#' @param stats_df dataframe of summary statistics
#' @param rhat_cols names of the columns to compute rhats for
#' @param district_idx list of row indices, one entry per district
#' @param district_groups list of `make_rhat_groups()` results, one per district
#' @param ncores number of cores to use; 1 runs serially
#'
#' @returns A districts by statistics matrix of rhats.
#' @noRd
rhats_over_districts <- function(stats_df, rhat_cols, district_idx, district_groups,
                                 ncores = 1) {
    row_index <- unlist(district_idx, use.names = FALSE)
    district_start <- c(0L, cumsum(lengths(district_idx)))

    # For each row of `row_index`, which group within its district it is in.
    # `make_rhat_groups()` already returns groups in the order R's `split()`
    # would produce them, so the numbering carries over directly.
    group_index <- integer(length(row_index))
    for (k in seq_along(district_groups)) {
        offset <- district_start[k]
        idx <- district_groups[[k]]$idx
        for (g in seq_along(idx)) {
            group_index[offset + idx[[g]]] <- g - 1L
        }
    }

    n_groups <- max(lengths(lapply(district_groups, `[[`, "idx")))

    compute_rhats_cpp(
        values = lapply(rhat_cols, function(col) as.numeric(stats_df[[col]])),
        row_index = as.integer(row_index) - 1L,
        district_start = as.integer(district_start),
        group_index = group_index,
        n_groups = as.integer(n_groups),
        num_threads = as.integer(ncores)
    )
}


#' Computes all rhats for a dataframe of summary statistics
#'
#' Given a dataframe with `chain` and `district` columns and all other columns
#' being statistics of interest this computes the rhats for each of the districts.
#' If `order_stats` is true then computes the rhats on the order statistics for
#' each plan ie `district == 1` represents the smallest value for a plan.
#'
#' @param stats_df A dataframe with `chain` and `district` columns along with
#' all of the `rhat_cols`
#' @param rhat_cols Vector of names of columns to compute rhats for
#' @param split_rhat Whether or not to split the chains in rhat calculations
#' @inheritParams summary.redist_plans
#'
#' @returns A long dataframe of rhats
#' @noRd
compute_all_rhats <- function(stats_df, rhat_cols, order_stats, district, ndists,
                              split_rhat, ncores = 1) {
    # order values if needed
    if (order_stats) {
        stats_df <- order_columns_by_district(
            stats_df |> filter(!is.na(chain)),
            rhat_cols,
            ndists,
            num_threads = ncores
        )
    }
    # filter district if needed
    if (!isFALSE(district)) {
        stats_df <- stats_df[stats_df$district %in% district, ]
    }
    # Split the row indices rather than the data frame itself. Splitting the
    # frame copies every column, and it used to happen once per statistic.
    district_idx <- split(seq_len(nrow(stats_df)), stats_df$district)
    chain <- stats_df$chain

    # The chain grouping depends only on the district, so build it once and
    # share it across every statistic instead of rebuilding the factor inside
    # each `tapply()` call. Building it here rather than in C++ also keeps the
    # split-chain convention in exactly one place.
    district_groups <- lapply(district_idx, function(rows) {
        make_rhat_groups(chain[rows], split = split_rhat)
    })

    rhat_matrix <- rhats_over_districts(
        stats_df, rhat_cols, district_idx, district_groups, ncores
    )

    rhat_results <- lapply(seq_along(rhat_cols), function(j) {
        stats::setNames(rhat_matrix[, j], names(district_idx))
    })
    names(rhat_results) <- rhat_cols

    rhats_df <- do.call(
        rbind,
        lapply(rhat_cols, function(col_name) {
            data.frame(
                district = rhat_results[[col_name]] |> names() |> as.integer(),
                stat_name = col_name,
                rhat = rhat_results[[col_name]],
                row.names = NULL
            )
        })
    )

    rhats_df
}


#' Pick the columns worth computing rhats for
#'
#' Returns the summary statistic columns that actually vary, dropping index
#' columns, columns that are entirely `NA`, columns that are constant, and
#' columns that are constant within any single chain. Returns an empty vector
#' when there is no `chain` column or only one chain, in which case no rhats
#' can be computed at all.
#'
#' @param object a [redist_plans] object
#' @param district `FALSE` for all districts, or a single district to restrict to
#' @param n_distr the number of districts
#' @param n_samp the number of sampled plans
#'
#' @returns A character vector of column names, possibly empty.
#' @noRd
select_rhat_columns <- function(object, district, n_distr, n_samp) {
    cols <- names(object)

    multiple_chains <- "chain" %in%
        cols &&
        dplyr::n_distinct(object[["chain"]]) > 1

    if (!multiple_chains) {
        return(character(0))
    }

    one_district_only <- 1 <= district && district <= n_distr

    addl_cols <- setdiff(
        cols,
        c("chain", "draw", "district", "total_pop", "seats", "mcmc_accept")
    )

    # When restricted to one district, look only at that district's rows.
    if (one_district_only) {
        idx <- seq_len(n_samp)
        if ("district" %in% cols) {
            idx <- as.integer(district) + (idx - 1) * n_distr
        }
        column_values <- function(col) object[[col]][idx]
        chain_ids <- object[["chain"]][idx]
    } else {
        column_values <- function(col) object[[col]]
        chain_ids <- object[["chain"]]
    }

    # `length(unique(x)) <= 1` covers both an all-NA column and a constant
    # one, and unlike `all(x == x[1])` it cannot itself be NA when the column
    # holds some but not all NAs. `na.rm` on the `any()` keeps the predicate
    # total even if a chain contributes no rows.
    const_cols <- vapply(
        addl_cols,
        function(col) {
            x <- column_values(col)
            length(unique(x)) <= 1 ||
                any(
                    tapply(x, chain_ids, FUN = function(z) length(unique(z))) ==
                        1,
                    na.rm = TRUE
                )
        },
        logical(1)
    )

    addl_cols[!const_cols]
}


#' Validate the `district` argument used for rhat computation
#'
#' @inheritParams select_rhat_columns
#' @param call environment to report the error against
#'
#' @returns `FALSE`, or `district` coerced to integer.
#' @noRd
validate_rhat_district <- function(district, n_distr, call = rlang::caller_env()) {
    if (isFALSE(district)) {
        return(district)
    }

    if (!rlang::is_integerish(district)) {
        cli::cli_abort("{.arg district} must be an integer!", call = call)
    }

    district <- as.integer(district)

    if (!all(1 <= district && district <= n_distr)) {
        cli::cli_abort(
      "{.arg district} must be between 1 and {.arg ndists}!",
      call = call
    )
    }

    district
}


#' Resolve the rhat reporting thresholds
#'
#' @inheritParams summary.redist_plans
#'
#' @returns A list with `q99` and `max` entries.
#' @noRd
resolve_rhat_thresh <- function(rhat_thresh) {
    thresh_names <- names(rhat_thresh)

    list(
    q99 = if ("q99" %in% thresh_names) rhat_thresh[["q99"]] else 1.05,
    max = if ("max" %in% thresh_names) rhat_thresh[["max"]] else 1.1
  )
}


#' Drop rhats that are `NA`, reporting how many were dropped
#'
#' @param rhats_df long dataframe of rhats
#'
#' @returns `rhats_df` with `NA` rhats removed.
#' @noRd
drop_na_rhats <- function(rhats_df) {
    if (anyNA(rhats_df$rhat)) {
        num_na_rhats <- sum(is.na(rhats_df$rhat))
        cat_cli("{num_na_rhats} rhat values are `NA`")
        rhats_df <- rhats_df[!is.na(rhats_df$rhat), ]
    }

    rhats_df
}


#' Classify chain convergence from a vector of rhats
#'
#' Strong convergence is all rhats at or below 1.05. Weak convergence is a
#' 99th percentile at or below 1.05 with all rhats at or below 1.1.
#'
#' @param rhat_vals numeric vector of rhats
#' @param thresh the result of `resolve_rhat_thresh()`
#'
#' @returns A list with `status` (one of `"strong"`, `"weak"`,
#' `"not_converged"`) and `q99`, the 99th percentile of the rhats.
#' @noRd
classify_rhat_convergence <- function(rhat_vals, thresh) {
    q99_rhat <- quantile(x = rhat_vals, probs = 0.99) |>
        unname()

    status <- if (all(rhat_vals <= thresh$q99)) {
        "strong"
    } else if (q99_rhat <= thresh$q99 && all(rhat_vals <= thresh$max)) {
        "weak"
    } else {
        "not_converged"
    }

    list(status = status, q99 = q99_rhat)
}


#' Print the rhat table, the threshold breakdown, and any convergence alert
#'
#' @param rhats_df long dataframe of rhats, with `NA`s already removed
#' @param convergence the result of `classify_rhat_convergence()`
#' @param thresh the result of `resolve_rhat_thresh()`
#' @inheritParams summary.redist_plans
#'
#' @returns `NULL`, invisibly. Called for its printed output.
#' @noRd
print_rhat_summary <- function(rhats_df, convergence, thresh, order_stats) {
    # cli tracks open containers per frame, so opening the list from this
    # frame rather than the caller's would indent it one level deeper than it
    # was when these calls lived inline in `summary.redist_plans()`.
    caller_env <- parent.frame()

    q99_rhat_thresh <- thresh$q99
    rhat_max_thresh <- thresh$max

    ordered_str <- ifelse(order_stats, "ordered ", "")
    cat_cli("Largest R-hat values for {ordered_str}summary statistics:\n")

    # get maximum rhats for each statistic
    max_rhats <- tapply(rhats_df$rhat, rhats_df$stat_name, max, na.rm = TRUE)

    rhats_p <- vapply(
        max_rhats,
        function(x) {
            ifelse(x <= q99_rhat_thresh, sprintf("%.3f", x), paste0("\U274C", round(x, 3)))
        },
        FUN.VALUE = character(1)
    )
    print(noquote(rhats_p))

    # print counts
    rhat_vals <- rhats_df$rhat

    cat_cli_fmt({
        cli::cli_ul(.envir = caller_env)
        cli::cli_li(
            "R-hat \u2264 {format(q99_rhat_thresh, digits=3)}: {sum(rhat_vals <= q99_rhat_thresh)}"
        )
        cli::cli_li(
            "{format(q99_rhat_thresh, digits=3)} < R-hat \u2264 {format(rhat_max_thresh, digits=3)}:
                        {sum(rhat_vals > q99_rhat_thresh & rhat_vals <= rhat_max_thresh)}"
        )
        cli::cli_li(
            "R-hat > {format(rhat_max_thresh, digits=3)}: {sum(rhat_vals > rhat_max_thresh)}"
        )
        cli::cli_li("Total R-hats: {length(rhat_vals)}")
        cli::cli_end()
    })

    if (convergence$status == "weak") {
        cat_cli_fmt(cli::cli_alert_info("{.strong ALERT:} Chains have weakly converged."))
    } else if (convergence$status == "not_converged") {
        cat_cli_fmt(cli::cli_alert_danger("{.strong WARNING:} Chains have not converged."))
    }

    invisible(NULL)
}


#' Legacy code to print diagnostic informaiton for old (pre Redist 5.0) plans
#'
#' @noRd
legacy_print_smc_information <- function(
    name,
    all_runs,
    object,
    algo,
    div_rg,
    all_diagn,
    warn_converge,
    div_bad
) {
    if (algo == "smc") {
        run_dfs <- list()
        n_runs <- length(all_diagn)
        warn_bottlenecks <- FALSE

        for (i in seq_len(n_runs)) {
            diagn <- all_diagn[[i]]
            n_samp <- nrow(diagn$ancestors)

            run_dfs[[i]] <- tibble(
                n_eff = c(diagn$step_n_eff, diagn$n_eff),
                eff = c(diagn$step_n_eff, diagn$n_eff) / n_samp,
                accept_rate = c(diagn$accept_rate, NA),
                sd_log_wgt = diagn$sd_lp,
                max_unique = diagn$unique_survive,
                est_k = c(diagn$est_k, NA)
            )

            tbl_print <- as.data.frame(run_dfs[[i]])
            min_n <- max(0.05 * n_samp, min(0.4 * n_samp, 100))
            bottlenecks <- dplyr::coalesce(
                with(tbl_print, pmin(max_unique, n_eff) < min_n),
                FALSE
            )
            warn_bottlenecks <- warn_bottlenecks || any(bottlenecks)
            tbl_print$bottleneck <- ifelse(bottlenecks, " * ", "")
            tbl_print$n_eff <- with(
                tbl_print,
                str_glue("{fmt_comma(n_eff)} ({sprintf('%0.1f%%', 100*eff)})")
            )
            tbl_print$eff <- NULL
            tbl_print$accept_rate <- with(tbl_print, sprintf("%0.1f%%", 100 * accept_rate))
            max_pct <- with(tbl_print, max_unique / (-n_samp * expm1(-1)))
            tbl_print$max_unique <- with(
                tbl_print,
                str_glue("{fmt_comma(max_unique)} ({sprintf('%3.0f%%', 100*max_pct)})")
            )

            names(tbl_print) <- c("Eff. samples (%)", "Acc. rate",
                                  "Log wgt. sd", " Max. unique",
                                  "Est. k", "")
            rownames(tbl_print) <- c(paste("Split", seq_len(nrow(tbl_print) - 1)), "Resample")

            if (i == 1 || isTRUE(all_runs)) {
                cat_cli("Sampling diagnostics for SMC run {i} of {n_runs} ({fmt_comma(n_samp)} samples)")
                print(tbl_print, digits = 2)
                cat("\n")
            }
        }
        out <- bind_rows(run_dfs)

        cat_cli_fmt(cli::cli_li(cli::col_grey(
            "
                Watch out for low effective samples, very low acceptance rates (less than 1%),
                large std. devs. of the log weights (more than 3 or so),
                and low numbers of unique plans.
                R-hat values for summary statistics should be between 1 and 1.05."
        )))

        if (div_bad) {
            cat_cli_fmt(cli::cli_li(
                "{.strong Low diversity:} Check for potential bottlenecks.
                            Increase the number of samples.
                            Examine the diversity plot with
                            `hist(plans_diversity({name}), breaks=24)`.
                            Consider weakening or removing constraints, or increasing
                            the population tolerance. If the acceptance rate drops
                            quickly in the final splits, try increasing
                            {.arg pop_temper} by 0.01."
            ))
        }
        if (warn_converge) {
            cat_cli_fmt(cli::cli_li(
                "{.strong SMC convergence:} Increase the number of samples.
                            If you are experiencing low plan diversity or bottlenecks as well,
                            address those issues first."
            ))
        }
        if (warn_bottlenecks) {
            cat_cli_fmt(cli::cli_li(
                "(*) {.strong Bottlenecks found:} Consider weakening or removing
                            constraints, or increasing the population tolerance.
                            If the acceptance rate drops quickly in the final splits,
                            try increasing {.arg pop_temper} by 0.01.
                            If the weight variance (Log wgt. sd) increases steadily
                            or is particularly large for the \"Resample\" step,
                            consider increasing {.arg seq_alpha}.
                            To visualize what geographic areas may be causing problems,
                            try running the following code. Highlighted areas are
                            those that may be causing the bottleneck.\n\n"
            ))
            code <- str_glue(
                "plot(<map object>, rowMeans(as.matrix({name}) == <bottleneck iteration>))"
            )
            cli::cat_line("    ", cli::code_highlight(code, "Material"))
        }
    } else if (algo %in% c("mergesplit", "flip")) {
        accept_rate <- sprintf("%0.1f%%", 100 * attr(object, "mh_acceptance"))
        cat_cli("Chain acceptance rate{?s}: {accept_rate}")

        cat_cli("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
        if (div_bad) {
            cat_cli_fmt(cli::cli_alert_danger("{.strong WARNING:} Low plan diversity"))
        }
        cat("\n")

        out <- tibble(
            accept_rate = attr(object, "mh_acceptance"),
            div_q10 = div_rg[1],
            div_q90 = div_rg[2]
        )

        cat_cli_fmt(cli::cli_li(cli::col_grey(
            "
                Watch out for low acceptance rates (less than 10%).
                R-hat values for summary statistics should be between 1 and 1.05.
                For district-level statistics (like district partisan leans), you
                should call `match_numbers()` or `number_by()` before examining
                the R-hat values."
        )))

        if (div_bad) {
            cat_cli_fmt(cli::cli_li(
                "{.strong Low diversity:} Increase the number of samples.
                            Examine the diversity plot with
                            `hist(plans_diversity({name}), breaks=24)`.
                            Consider weakening or removing constraints, or increasing
                            the population tolerance."
            ))
        }
        if (warn_converge) {
            cat_cli_fmt(cli::cli_li(
                "{.strong Chain convergence:} Increase the number of samples.
                            If you are experiencing low plan diversity, address that issue first."
            ))
        }
    }

    out
}
