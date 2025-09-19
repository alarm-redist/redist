
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
#' * **Effective samples**: the effective sample size at each iteration, computed
#' using the SMC weights. Larger is better. The percentage in parentheses is the
#' ratio of the effective samples to the total samples.
#' * **Acceptance rate**: the fraction of drawn spanning trees which yield a valid
#' redistricting plan within the population tolerance. Very small values (< 1%)
#' can indicate a bottleneck and may lead to a lack of diversity.
#' * **Standard deviation of the log weights**: More variable weights (larger s.d.)
#' indicate less efficient sampling. Values greater than 3 are likely problematic.
#' * **Maximum unique plans:** The number of unique redistricting
#' plans that survive each stage. The percentage in parentheses is the ratio of
#' this number to the total number of samples. Small values (< 100) indicate a
#' bottleneck, which leads to a loss of sample diversity and a higher variance.
#' * **Estimated `k` parameter**: For graph space plans displays how many
#' spanning tree edges were considered for cutting at each split. Mostly
#' informational, though large jumps may indicate a need to increase
#' `adapt_k_thresh`.
#' * **Bottleneck**: An asterisk will appear in the right column if a bottleneck
#' appears likely, based on the values of the other statistics.
#'
#' In the event of problematic diagnostics, the function will provide
#' suggestions for improvement.
#'
#' @param object a [redist_plans] object
#' @param district For R-hat values, which district to use for district-level
#' summary statistics. We strongly recommend calling `match_numbers()` or
#' `number_by()` before examining these district-level statistics. A value of
#' `FALSE` indicates to compute the rhats for all districts
#' @param all_runs When there are multiple SMC runs, show detailed summary
#' statistics for all runs (the default), or only the first run?
#' @param vi_max The maximum number of plans to sample in computing the pairwise
#' variation of information distance (sample diversity).
#' @param order_stats Whether or not to compute rhats on the ordered district
#' statistics.
#' @param rhat_thresh What values to use when checking convergence. It should
#' consist of two named values: `q99` and `max`.
#' Strong convergence is when all rhats are less than or equal to the `max` value
#' and weak convergence is when the 99th quantile is less than or equal to
#' `q99` and all rhats are less than or equal to `max`.
#'
#' @param \dots additional arguments (ignored)
#'
#' @return A data frame containing diagnostic information, invisibly.
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
        object, district = FALSE, all_runs = TRUE, vi_max = 100,
        order_stats = TRUE,
        rhat_thresh = getOption("redist.rhat_thresh", c(q99=1.05, max=1.1)),
        ...) {
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
    if (is.null(n_distr)) n_distr <- max(plans_m[, 1])


    if (n_distr == 1 || nrow(plans_m) == 1) {
        cli::cli_text("{fmt_comma(n_samp)}{cli::qty(n_samp)} sampled plan{?s} of
                 {n_distr} district{?s} on
                 {fmt_comma(nrow(plans_m))}{cli::qty(nrow(plans_m))} unit{?s}")
        return(invisible(1))
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

    summary_supported_algs <- c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE, MCMC_ALG_TYPE, "flip")

    # ignore if not a supported algorithm
    if(!algo %in% summary_supported_algs){
        cli::cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
        return(invisible(1))
    }

    # bool for not flip
    revamped_alg <- algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE, MCMC_ALG_TYPE)

    # checks if plans are earlier than version 5.0
    pre_v5_plans <- is.null(attr(object, "version")) || attr(object, "version") < '5'


    # get the display name
    alg_display_name <- dplyr::case_when(
        algo == SMC_ALG_TYPE ~ "SMC",
        algo == MS_SMC_ALG_TYPE ~ "SMC with Merge-Split MCMC Steps",
        algo == MCMC_ALG_TYPE ~ "Merge-Split MCMC",
        algo == "flip" ~ "Flip MCMC"
    )
    cli::cli_text("{.strong {alg_display_name}:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")

    if(pre_v5_plans){
        display_sampling_space <- GRAPH_PLAN_SPACE_SAMPLING
        display_splitting_method <- NAIVE_K_SPLITTING
        cli::cli_text("Plans sampled on {display_sampling_space} using the {display_splitting_method} forward kernel.")
        cli::cli_text("Forward kernel parameters: {.arg adapt_k_thresh}={format(all_diagn[[1]]$adapt_k_thresh, digits=3)}")
    }else if(revamped_alg && !pre_v5_plans){
        # if revamped alg check that sampling space and splitting methods all the same
        # check same sampling space
        all_sampling_spaces <- sapply(all_run_info, function(x) x$sampling_space)
        if(length(unique(all_sampling_spaces)) != 1){
            cli::cli_abort("{.fn summary} is not supported for plans sampled using different sampling spaces.")
            return(invisible(1))
        }
        sampling_space <- all_sampling_spaces[1]
        display_sampling_space <- dplyr::case_when(
            sampling_space == GRAPH_PLAN_SPACE_SAMPLING ~ "Graph Space",
            sampling_space == FOREST_SPACE_SAMPLING ~ "Spanning Forest Space",
            sampling_space == LINKING_EDGE_SPACE_SAMPLING ~ "Linking Edge Forest Space"
        )
        # check same splitting method
        all_splitting_methods <- sapply(all_run_info, function(x) x$split_method)
        if(length(unique(all_splitting_methods)) != 1){
            cli::cli_abort("{.fn summary} is not supported for plans sampled using different splitting methods")
            return(invisible(1))
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
            if(!identical(all_forward_kernel_params[[1]], all_forward_kernel_params[[i]])){
                cli::cli_abort("{.fn summary} is not supported for plans sampled using different splitting parameters")
                return(invisible(1))
            }
        }
        forward_kernel_params <- all_forward_kernel_params[[1]]

        cli::cli_text("Plans sampled on {display_sampling_space} using the {display_splitting_method} forward kernel.")

        if(split_method == NAIVE_K_SPLITTING){
            # only display adapt k threshold if k values estimated
            if(forward_kernel_params$estimate_cut_k){
                cli::cli_text("Forward kernel parameters: {.arg adapt_k_thresh}={format(forward_kernel_params$adapt_k_thresh, digits=3)}")
            }else{
                cli::cli_text("Forward kernel parameters: {.arg manual_k_params}={forward_kernel_params$manual_k_params}")
            }
        }else if(split_method == UNIF_VALID_EDGE_SPLITTING){

        }else if(split_method == EXP_BIGGER_ABS_DEV_SPLITTING){
            cli::cli_text("Forward Kernel Parameters: {.arg splitting_alpha}={format(forward_kernel_params$splitting_alpha, digits=2)}")
        }
    }


    # print algorithm specific parameters
    print_algo_specific_params(algo, all_diagn, all_run_info, pre_v5_plans)

    cli::cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
    if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")

    # now compute rhats if more than 1 chain
    cols <- names(object)
    multiple_chains <- "chain" %in% cols && dplyr::n_distinct(object[["chain"]]) > 1
    if(multiple_chains){
        one_district_only <- 1 <= district && district <= n_distr

        addl_cols <- setdiff(cols, c("chain", "draw", "district", "total_pop", "seats", "mcmc_accept"))
        if(one_district_only){
            idx <- seq_len(n_samp)
            if ("district" %in% cols) idx <- as.integer(district) + (idx - 1)*n_distr

            const_cols <- vapply(addl_cols, function(col) {
                x <- object[[col]][idx]
                all(is.na(x)) || all(x == x[1]) ||
                    any(tapply(x, object[['chain']][idx], FUN = function(z) length(unique(z))) == 1)
            }, numeric(1))
        }else{
            const_cols <- vapply(addl_cols, function(col) {
                x <- object[[col]]
                all(is.na(x)) || all(x == x[1]) ||
                    any(tapply(x, object[['chain']], FUN = function(z) length(unique(z))) == 1)
            }, numeric(1))
        }
        addl_cols <- addl_cols[!const_cols]
    }else{
        addl_cols <- c()
    }

    warn_converge <- FALSE
    # do nothing if no additional columns or no chain column


    if(multiple_chains && length(addl_cols) > 0){
        # check district input
        if(!isFALSE(district)){
            # check integer
            if(!rlang::is_integerish(district)){
                cli::cli_abort("{.arg district} must be an integer!")
            }else{
                district <- as.integer(district)
            }
            # check between 1 and ndists
            if(!all(1 <= district && district <= n_distr)){
                cli::cli_abort("{.arg district} must be between 1 and {.arg ndists}!")
            }
        }
        rhats_computed <- TRUE
        split_rhat = algo %in% c(MCMC_ALG_TYPE, "flip")

        # get rhats
        rhats_df <- compute_all_rhats(
            # drop everything but the columns to save size
            as.data.frame(object)[c("chain", "district", addl_cols)],
            addl_cols,
            order_stats,
            district,
            n_distr,
            split_rhat)

        # remove NA rhats
        if(any(is.na(rhats_df$rhat))){
            num_na_rhats <- sum(is.na(rhats_df$rhat))
            cli::cli_inform("{num_na_rhats} rhat values are `NA`")
            rhats_df <- rhats_df[!is.na(rhats_df$rhat),]
        }

        # get thresholds
        q99_rhat_thresh <- ifelse("q99" %in% rhat_thresh, rhat_thresh[["q99"]], 1.05)
        rhat_max_thresh <- ifelse("max" %in% rhat_thresh, rhat_thresh[["max"]], 1.1)


        ordered_str <- ifelse(order_stats, "ordered ", "")
        cli::cli_text("Largest R-hat values for {ordered_str}summary statistics:\n")
        # get maximum rhats for each statistic
        max_rhats <- tapply(rhats_df$rhat, rhats_df$stat_name, max, na.rm = TRUE)

        rhats_p <- vapply(max_rhats, function(x){
            ifelse(x <= q99_rhat_thresh, sprintf('%.3f', x), paste0('\U274C', round(x, 3)))
        }, FUN.VALUE = character(1))
        print(noquote(rhats_p))

        # print counts
        rhat_vals <- rhats_df$rhat



        cli::cli_ul()
        cli::cli_li("R-hat ≤ {format(q99_rhat_thresh, digits=3)}: {sum(rhat_vals <= q99_rhat_thresh)}")
        cli::cli_li("{format(q99_rhat_thresh, digits=3)} < R-hat ≤ {format(rhat_max_thresh, digits=3)}:
                    {sum(rhat_vals > q99_rhat_thresh & rhat_vals <= rhat_max_thresh)}")
        cli::cli_li("R-hat > {format(rhat_max_thresh, digits=3)}: {sum(rhat_vals > rhat_max_thresh)}")
        cli::cli_li("Total R-hats: {length(rhat_vals)}")
        cli::cli_end()


        # cat("Rhat Breakdown:\n")
        # cat("R-hat ≤ 1.05:      ", sum(rhat_vals <= 1.05), "\n")
        # cat("1.05 < R-hat ≤ 1.1:", sum(rhat_vals > 1.05 & rhat_vals <= 1.1), "\n")
        # cat("R-hat > 1.1:       ", sum(rhat_vals > 1.1), "\n")


        # get 99th quantile
        q99_rhat <- quantile(x = rhats_df$rhat, probs = .99) |>
            unname()

        # check converge
        # - weak convergence: all rhats <= 1.1 and 99th quantile <= 1.05
        # - strong covergence: all rhats <= 1.05
        if(all(rhat_vals <= 1.05)){
            convergence_status <- "strong"
        }else if(q99_rhat <= 1.05 && all(rhat_vals <= 1.1)){
            convergence_status <- "weak"
            cli::cli_alert_info("{.strong ALERT:} Chains have weakly converged.")
        }else{
            convergence_status <- "not_converged"
            warn_converge <- TRUE
            cli::cli_alert_danger("{.strong WARNING:} Chains have not converged.")
        }
    }else{
        rhats_computed <- FALSE
    }


    # Now print algorithm specific diagnostics
    if (algo == "smc" && pre_v5_plans){
        out <- legacy_print_smc_information(
            name, all_runs, object, algo, div_rg,
            all_diagn, warn_converge, div_bad)
    }else if (algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE)) {
        smc_run_dfs <- list()
        smc_ms_run_dfs <- list()
        n_runs <- length(all_diagn)
        warn_bottlenecks <- FALSE

        for (i in seq_len(n_runs)) {

            smc_summary_result_list <- get_smc_summary_df(
                all_diagn[[i]], all_run_info[[i]],
                resampled, warn_bottlenecks
            )

            smc_run_dfs[[i]] <- smc_summary_result_list$smc_summary_df
            smc_tbl_print <- smc_summary_result_list$smc_print_tbl
            warn_bottlenecks <- smc_summary_result_list$warn_bottlenecks

            if (i == 1 || isTRUE(all_runs)) {
                cli::cli_text("Sampling diagnostics for SMC run {i} of {n_runs} ({fmt_comma(n_samp)} samples)")
                print(smc_tbl_print, digits = 2)
                cat("\n")
            }

            if(algo == MS_SMC_ALG_TYPE){
                smc_ms_summary_result_list <- get_smc_ms_summary_df(
                    all_diagn[[i]], all_run_info[[i]],
                    resampled
                )
                smc_ms_run_dfs[[i]] <- smc_ms_summary_result_list$smc_ms_summary_df
                smc_ms_tbl_print <- smc_ms_summary_result_list$smc_ms_print_tbl

                if (i == 1 || isTRUE(all_runs)) {
                    cli::cli_text("Sampling diagnostics for Mergesplit Steps of SMC run {i} of {n_runs} ({fmt_comma(n_samp)} samples)")
                    print(smc_ms_tbl_print, digits = 2)
                    cat("\n")
                }
            }

        }
        out <- bind_rows(smc_ms_run_dfs)

        #step_nums <- ave(seq_along(all_run_info[[i]]$step_types), all_run_info[[i]]$step_types, FUN = seq_along)

        cli::cli_li(cli::col_grey("
            Watch out for low effective samples, very low acceptance rates (less than 1%),
            large std. devs. of the log weights (more than 3 or so),
            and low numbers of unique plans.
            R-hat values for summary statistics should be between 1 and 1.05."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Check for potential bottlenecks.
                        Increase the number of samples.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance. If the acceptance rate drops
                        quickly in the final splits, try increasing
                        {.arg pop_temper} by 0.01.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong SMC convergence:} Increase the number of samples.
                        If you are experiencing low plan diversity or bottlenecks as well,
                        address those issues first.")
        }
        if (warn_bottlenecks) {
            cli::cli_li("(*) {.strong Bottlenecks found:} Consider weakening or removing
                        constraints, or increasing the population tolerance.
                        If the acceptance rate drops quickly in the final splits,
                        try increasing {.arg pop_temper} by 0.01.
                        If the weight variance (Log wgt. sd) increases steadily
                        or is particularly large for the \"Resample\" step,
                        consider increasing {.arg seq_alpha}.
                        To visualize what geographic areas may be causing problems,
                        try running the following code. Highlighted areas are
                        those that may be causing the bottleneck.\n\n")
            code <- str_glue("plot(<map object>, rowMeans(as.matrix({name}) == <bottleneck iteration>))")
            cli::cat_line("    ", cli::code_highlight(code, "Material"))
        }
    }else if (algo %in% c("mergesplit", 'flip')) {

        accept_rate <- sprintf("%0.1f%%", 100*attr(object, "mh_acceptance"))
        cli::cli_text("Chain acceptance rate{?s}: {accept_rate}")

        if(rhats_computed){
            out <- rhats_df
        }else{
            out <- tibble(accept_rate = attr(object, "mh_acceptance"),
                          div_q10 = div_rg[1],
                          div_q90 = div_rg[2])
        }

        cli::cli_li(cli::col_grey("
            Watch out for low acceptance rates (less than 10%).
            R-hat values for summary statistics should be between 1 and 1.05.
            For district-level statistics (like district partisan leans), you
            should call `match_numbers()` or `number_by()` before examining
            the R-hat values."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Increase the number of samples.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong Chain convergence:} Increase the number of samples.
                        If you are experiencing low plan diversity, address that issue first.")
        }
    } else {
        cli::cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
    }

    invisible(out)
}


fmt_comma <- function(x){
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
print_algo_specific_params <- function(algo, all_diagn, all_run_info, pre_v5_plans){

    if (algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE)) {
        if(!pre_v5_plans){
            weight_type <- all_run_info[[1]]$weight_type
        }else{
            weight_type <- "simple"
        }
        cli::cli_bullets(c(
            "SMC Parameters:",
            "*"="{.arg weight_type} = {weight_type}",
            "*"="{.arg seq_alpha} = {format(all_diagn[[1]]$seq_alpha, digits=2)}",
            "*"="{.arg pop_temper} = {format(all_diagn[[1]]$pop_temper, digits=3)}"
        )
        )
    }else if(algo == MCMC_ALG_TYPE){
        cli::cli_bullets(c(
            "MCMC Parameters:",
            "*"="{.arg warmup} = {format(all_diagn[[1]]$warmup)}",
            "*"="{.arg thin} = {format(all_diagn[[1]]$thin)}",
            "*"="{.arg total_steps} = {format(all_diagn[[1]]$total_steps)}"
        )
        )
    }
    if(algo == MS_SMC_ALG_TYPE){
        total_ms_steps <- sum(all_run_info[[1]]$step_types == "ms")

        cli::cli_bullets(c(
            "Mergesplit Parameters:",
            "*"="{.arg total_ms_steps} = {format(total_ms_steps, digits=2)}",
            "*"="{.arg mh_accept_per_smc} = {format(all_run_info[[1]]$mh_accept_per_smc, digits=2)}",
            "*"="{.arg pair_rule} = {format(all_run_info[[1]]$pair_rule, digits=2)}"
        )
        )
    }
}



#' Get Summary and Print Dataframe for SMC steps
#'
#' @noRd
get_smc_summary_df <- function(diagn, run_info, resampled, warn_bottlenecks){
    n_samp <- run_info$nsims
    run_sampling_space <- run_info$sampling_space

    smc_accept_rate <- diagn$accept_rate[run_info$step_types == "smc"]

    smc_step_mask <- run_info$step_types == "smc"

    run_summary_df <- tibble(
        n_eff = diagn$step_n_eff,
        eff = diagn$step_n_eff/n_samp,
        accept_rate = smc_accept_rate,
        sd_log_wgt = diagn$sd_lp,
        max_unique = diagn$unique_survive[seq_len(length(smc_accept_rate))]
    )

    run_summary_names <- c("Eff. samples (%)", "Acc. rate",
                           "Log wgt. sd", " Max. unique")

    # if graph space then add the k
    if(run_sampling_space == GRAPH_PLAN_SPACE_SAMPLING){
        run_summary_df$est_k <- diagn$forward_kernel_params$cut_k_used[smc_step_mask]
        if(diagn$forward_kernel_params$estimate_cut_k){
            run_summary_names <- c(run_summary_names, "Est. k")
        }else{
            run_summary_names <- c(run_summary_names, "Manual k")
        }
    }

    # add row if resampled
    if(resampled){
        run_summary_df[nrow(run_summary_df) + 1, ] <- NA
        # add eff sample size and how many survived resample
        run_summary_df[nrow(run_summary_df), "n_eff"] <- diagn$n_eff
        run_summary_df[nrow(run_summary_df), "eff"] <- diagn$n_eff/n_samp
        run_summary_df[nrow(run_summary_df), "max_unique"] <- tail(diagn$unique_survive, 1)
    }

    # add extra column name for asterisk
    run_summary_names <- c(run_summary_names, "")


    tbl_print <- as.data.frame(run_summary_df)
    min_n <- max(0.05*n_samp, min(0.4*n_samp, 100))
    bottlenecks <- dplyr::coalesce(with(tbl_print, pmin(max_unique, n_eff) < min_n), FALSE)
    warn_bottlenecks <- warn_bottlenecks || any(bottlenecks)
    tbl_print$bottleneck <- ifelse(bottlenecks, " * ", "")
    tbl_print$n_eff <- with(tbl_print,
                            str_glue("{fmt_comma(n_eff)} ({sprintf('%0.1f%%', 100*eff)})"))
    tbl_print$eff <- NULL
    tbl_print$accept_rate <- with(tbl_print, sprintf("%0.1f%%", 100*accept_rate))
    max_pct <- with(tbl_print, max_unique/(-n_samp * expm1(-1)))
    tbl_print$max_unique <- with(tbl_print,
                                 str_glue("{fmt_comma(max_unique)} ({sprintf('%3.0f%%', 100*max_pct)})"))

    names(tbl_print) <- run_summary_names
    new_row_names <- paste("Split", seq_len(nrow(tbl_print) - 1))

    if(resampled){
        new_row_names <- c(new_row_names, "Resample")
    }

    rownames(tbl_print) <- new_row_names

    list(
        smc_summary_df = run_summary_df,
        smc_print_tbl = tbl_print,
        warn_bottlenecks = warn_bottlenecks
    )

}


#' Get Summary and Print Dataframe for Mergesplit within SMC steps
#'
#' @noRd
get_smc_ms_summary_df <- function(diagn, run_info, resampled){
    n_samp <- run_info$nsims
    run_sampling_space <- run_info$sampling_space

    ms_accept_rate <- diagn$accept_rate[run_info$step_types == "ms"]
    ms_moves_per_plan <- diagn$ms_move_counts

    run_summary_df <- tibble(
        accept_rate = ms_accept_rate,
        ms_moves = ms_moves_per_plan
    )

    run_summary_names <- c("Acc. rate", "MS Moves per Plan")

    tbl_print <- as.data.frame(run_summary_df)

    tbl_print$accept_rate <- with(tbl_print, sprintf("%0.1f%%", 100*accept_rate))

    names(tbl_print) <- run_summary_names
    new_row_names <- paste("MS Step", seq_len(nrow(tbl_print)))


    rownames(tbl_print) <- new_row_names


    list(
        smc_ms_summary_df = run_summary_df,
        smc_ms_print_tbl = tbl_print
    )

}

diag_fold <- function(x) {
    abs(x - median(x))
}

diag_ranknorm <- function(x) {
    qnorm(rank(x)/(length(x) + 1))
}

diag_calc_rhat <- function(x, grp) {
    n <- mean(table(grp))
    var_between <- n*var(tapply(x, grp, mean))
    var_within <- mean(tapply(x, grp, var))
    sqrt((var_between/var_within + n - 1)/n)
}

diag_rhat <- function(x, grp, split = FALSE) {
    if (split) {
        lengths <- rle(grp)$lengths
        grp <- grp + do.call(c, lapply(lengths, function(l) rep(c(0.0, 0.5), each = l/2)))
    }

    max(diag_calc_rhat(diag_ranknorm(x), grp),
        diag_calc_rhat(diag_ranknorm(diag_fold(x)), grp))
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
compute_all_rhats <- function(stats_df, rhat_cols, order_stats, district, ndists, split_rhat){

    # order values if needed
    if(order_stats){
        stats_df <- order_columns_by_district(
            stats_df |> filter(!is.na(chain)),
            rhat_cols, ndists
        )
    }
    # filter district if needed
    if(!isFALSE(district)){
        stats_df <- stats_df[stats_df$district %in% district,]
    }
    # now compute rhats for each column and district
    rhat_results <- lapply(rhat_cols, function(col_name) {
        # Split data by district
        # For each district, compute rhat for the column
        sapply(split(stats_df, stats_df$district), function(df) {
            diag_rhat(x = df[[col_name]], grp = df$chain, split = split_rhat)
        })
    })
    names(rhat_results) <- rhat_cols

    rhats_df <- do.call(
        rbind,
        lapply(rhat_cols,
               function(col_name) data.frame(
                   district = rhat_results[[col_name]] |> names() |> as.integer(),
                   stat_name = col_name,
                   rhat = rhat_results[[col_name]],
                   row.names = NULL
               ))
    )

    rhats_df

}

#' Legacy code to print diagnostic informaiton for old (pre Redist 5.0) plans
#'
#' @noRd
legacy_print_smc_information <- function(name, all_runs, object, algo, div_rg, all_diagn, warn_converge, div_bad){
    if (algo == "smc") {
        run_dfs <- list()
        n_runs <- length(all_diagn)
        warn_bottlenecks <- FALSE

        for (i in seq_len(n_runs)) {
            diagn <- all_diagn[[i]]
            n_samp <- nrow(diagn$ancestors)

            run_dfs[[i]] <- tibble(n_eff = c(diagn$step_n_eff, diagn$n_eff),
                                   eff = c(diagn$step_n_eff, diagn$n_eff)/n_samp,
                                   accept_rate = c(diagn$accept_rate, NA),
                                   sd_log_wgt = diagn$sd_lp,
                                   max_unique = diagn$unique_survive,
                                   est_k = c(diagn$est_k, NA))

            tbl_print <- as.data.frame(run_dfs[[i]])
            min_n <- max(0.05*n_samp, min(0.4*n_samp, 100))
            bottlenecks <- dplyr::coalesce(with(tbl_print, pmin(max_unique, n_eff) < min_n), FALSE)
            warn_bottlenecks <- warn_bottlenecks || any(bottlenecks)
            tbl_print$bottleneck <- ifelse(bottlenecks, " * ", "")
            tbl_print$n_eff <- with(tbl_print,
                                    str_glue("{fmt_comma(n_eff)} ({sprintf('%0.1f%%', 100*eff)})"))
            tbl_print$eff <- NULL
            tbl_print$accept_rate <- with(tbl_print, sprintf("%0.1f%%", 100*accept_rate))
            max_pct <- with(tbl_print, max_unique/(-n_samp * expm1(-1)))
            tbl_print$max_unique <- with(tbl_print,
                                         str_glue("{fmt_comma(max_unique)} ({sprintf('%3.0f%%', 100*max_pct)})"))

            names(tbl_print) <- c("Eff. samples (%)", "Acc. rate",
                                  "Log wgt. sd", " Max. unique",
                                  "Est. k", "")
            rownames(tbl_print) <- c(paste("Split", seq_len(nrow(tbl_print) - 1)), "Resample")

            if (i == 1 || isTRUE(all_runs)) {
                cli_text("Sampling diagnostics for SMC run {i} of {n_runs} ({fmt_comma(n_samp)} samples)")
                print(tbl_print, digits = 2)
                cat("\n")
            }
        }
        out <- bind_rows(run_dfs)

        cli::cli_li(cli::col_grey("
                Watch out for low effective samples, very low acceptance rates (less than 1%),
                large std. devs. of the log weights (more than 3 or so),
                and low numbers of unique plans.
                R-hat values for summary statistics should be between 1 and 1.05."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Check for potential bottlenecks.
                            Increase the number of samples.
                            Examine the diversity plot with
                            `hist(plans_diversity({name}), breaks=24)`.
                            Consider weakening or removing constraints, or increasing
                            the population tolerance. If the acceptance rate drops
                            quickly in the final splits, try increasing
                            {.arg pop_temper} by 0.01.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong SMC convergence:} Increase the number of samples.
                            If you are experiencing low plan diversity or bottlenecks as well,
                            address those issues first.")
        }
        if (warn_bottlenecks) {
            cli::cli_li("(*) {.strong Bottlenecks found:} Consider weakening or removing
                            constraints, or increasing the population tolerance.
                            If the acceptance rate drops quickly in the final splits,
                            try increasing {.arg pop_temper} by 0.01.
                            If the weight variance (Log wgt. sd) increases steadily
                            or is particularly large for the \"Resample\" step,
                            consider increasing {.arg seq_alpha}.
                            To visualize what geographic areas may be causing problems,
                            try running the following code. Highlighted areas are
                            those that may be causing the bottleneck.\n\n")
            code <- str_glue("plot(<map object>, rowMeans(as.matrix({name}) == <bottleneck iteration>))")
            cli::cat_line("    ", cli::code_highlight(code, "Material"))
        }
    } else if (algo %in% c("mergesplit", 'flip')) {
        accept_rate <- sprintf("%0.1f%%", 100*attr(object, "mh_acceptance"))
        cli_text("Chain acceptance rate{?s}: {accept_rate}")

        cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
        if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")
        cat("\n")

        out <- tibble(accept_rate = attr(object, "mh_acceptance"),
                      div_q10 = div_rg[1],
                      div_q90 = div_rg[2])


        cli::cli_li(cli::col_grey("
                Watch out for low acceptance rates (less than 10%).
                R-hat values for summary statistics should be between 1 and 1.05.
                For district-level statistics (like district partisan leans), you
                should call `match_numbers()` or `number_by()` before examining
                the R-hat values."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Increase the number of samples.
                            Examine the diversity plot with
                            `hist(plans_diversity({name}), breaks=24)`.
                            Consider weakening or removing constraints, or increasing
                            the population tolerance.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong Chain convergence:} Increase the number of samples.
                            If you are experiencing low plan diversity, address that issue first.")
        }
    }

    out
}


#' Get k-step Ancestors of particles
#'
#' Tells you for a given index what its original ancestor was.
#'
#' @param parent_mat Ancestor matrix where entry [i,j] equals the index of the
#' parent of particle i after step j
#' @param steps_back How many steps back we should find the ancestors of (so
#' parents are 1 step ancestors). Defaults to going all the way back to the
#' original ancestors
#' @param start_col Which column to start at. Defaults to the last column
#'
#' @returns indices of k-step ancestors for particles at iteration `start_col`
#' @noRd
get_k_step_ancestors <- function(parent_mat, steps_back = NULL, start_col = NULL){
    # check the matrix is not zero indexed
    if(0 %in% parent_mat){
        cli::cli_abort("parent_mat must be 1-indexed, not 0 indexed!")
    }


    nparent_cols <- ncol(parent_mat)

    # if null then just start on final set of plans
    if(is.null(start_col)){
        start_col <- nparent_cols
    }else{
        # else assert starting column is between 1 and number of cols
        if(!(rlang::is_scalar_integerish(start_col) &&
             start_col <= nparent_cols &&
             start_col > 1)){
            cli::cli_abort("Input start_col={start_col} is not valid.
                          Input must be a number between 1 and {nparent_cols}")
        }
    }

    # check steps back is between [1, nnparent_cols -1]
    if(is.null(steps_back)){
        steps_back <- nparent_cols-1
    }else{
        if(!(rlang::is_scalar_integerish(steps_back) &&
             steps_back <= start_col-1 &&
             steps_back >= 1)){
            cli::cli_abort("Input steps_back={steps_back} is not valid.
Input must be between 1 and the start_col value (you input {steps_back})")
        }
    }

    # vector where index i maps to the index of its ancestor
    # initialize to this
    ancestor <- seq_len(nrow(parent_mat))


    # iterate through each step back we select the successive parent indices
    for (j in start_col:(start_col - steps_back + 1)) {
        ancestor <- parent_mat[ancestor, j]
    }
    ancestor
}



#' Get Original Ancestor Matrix of particles
#'
#' Gets a matrix of the original ancestors (ie first splits) of the particles
#' at each step.
#'
#' @param parent_mat Ancestor matrix where entry [i,j] equals the index of the
#' parent of particle i after step j
#'
#' @returns indices of originals ancestors for particles at every iteration. So
#' entry [s,j] is the original ancestor of particle j after step s
#' @noRd
get_original_ancestors_mat <- function(parent_mat){

    #if the matrix has one column then its just itself
    if(ncol(parent_mat) == 1){
        return(parent_mat)
    }

    # get the original ancestors at every step from the parent matrix
    original_ancestor_mat <- sapply(
        2:ncol(parent_mat),
        function(col_num) get_k_step_ancestors(parent_mat, steps_back = col_num-1, start_col = col_num)
    )

    # add the first column where every particles ancestor is iteslf
    cbind(
        seq_len(nrow(parent_mat)),
        original_ancestor_mat
    )
}
