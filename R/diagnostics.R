
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
#' * **Maximum unique plans:** an upper bound on the number of unique redistricting
#' plans that survive each stage. The percentage in parentheses is the ratio of
#' this number to the total number of samples. Small values (< 100) indicate a
#' bottleneck, which leads to a loss of sample diversity and a higher variance.
#' * **Estimated `k` parameter**: How many spanning tree edges were considered for
#' cutting at each split. Mostly informational, though large jumps may indicate
#' a need to increase `adapt_k_thresh`.
#' * **Bottleneck**: An asterisk will appear in the right column if a bottleneck
#' appears likely, based on the values of the other statistics.
#'
#' In the event of problematic diagnostics, the function will provide
#' suggestions for improvement.
#'
#' @param object a [redist_plans] object
#' @param district For R-hat values, which district to use for district-level
#' summary statistics. We strongly recommend calling `match_numbers()` or
#' `number_by()` before examining these district-level statistics.
#' @param all_runs When there are multiple SMC runs, show detailed summary
#' statistics for all runs (the default), or only the first run?
#' @param vi_max The maximum number of plans to sample in computing the pairwise
#' variation of information distance (sample diversity).
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
summary.redist_plans <- function(object, district = 1L, all_runs = TRUE, vi_max = 100, ...) {
    cli::cli_process_done(done_class = "") # in case an earlier

    algo <- attr(object, "algorithm")
    name <- deparse(substitute(object))

    object <- subset_sampled(object)
    all_diagn <- attr(object, "diagnostics")
    plans_m <- get_plans_matrix(object)
    n_samp <- ncol(plans_m)
    n_distr <- attr(object, "ndists")
    if (is.null(n_distr)) n_distr <- max(plans_m[, 1])

    fmt_comma <- function(x) format(x, nsmall = 0, digits = 1, big.mark = ",")
    if (n_distr == 1 || nrow(plans_m) == 1) {
        cli_text("{fmt_comma(n_samp)}{cli::qty(n_samp)} sampled plan{?s} of
                 {n_distr} district{?s} on
                 {fmt_comma(nrow(plans_m))}{cli::qty(nrow(plans_m))} unit{?s}")
        return(invisible(1))
    }

    prec_pop <- attr(object, "prec_pop")
    if (is.null(prec_pop)) {
        cli_warn(c("Precinct population missing; plan diversity estimates may be misleading.",
                   ">" = 'Run `attr({name}, "prec_pop") <- <map object>$<pop column>` to fix.'))
        prec_pop <- rep(1, nrow(plans_m))
    }
    est_div <- plans_diversity(object, total_pop = prec_pop, n_max = vi_max)
    div_rg <- format(quantile(est_div, c(0.1, 0.9)), digits = 2)
    div_bad <- (mean(est_div) <= 0.35) || (mean(est_div <= 0.25) > 0.1)

    if (algo == "smc") {
        cli_text("{.strong SMC:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")
        cli_text("{.arg adapt_k_thresh}={format(all_diagn[[1]]$adapt_k_thresh, digits=3)} \u2022
                 {.arg seq_alpha}={format(all_diagn[[1]]$seq_alpha, digits=2)}")
        cli_text("{.arg pop_temper}={format(all_diagn[[1]]$pop_temper, digits=3)}")
        cat("\n")

        cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
        if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")
        cat("\n")

        cols <- names(object)
        addl_cols <- setdiff(cols, c("chain", "draw", "district", "total_pop"))
        warn_converge <- FALSE
        if ("chain" %in% cols && length(addl_cols) > 0) {
            idx <- seq_len(n_samp)
            if ("district" %in% cols) idx <- as.integer(district) + (idx - 1)*n_distr

            const_cols <- vapply(addl_cols, function(col) {
                x <- object[[col]][idx]
                all(is.na(x)) || all(x == x[1]) ||
                    any(tapply(x, object[['chain']][idx], FUN = function(z) length(unique(z))) == 1)
            }, numeric(1))
            addl_cols <- addl_cols[!const_cols]

            rhats <- vapply(addl_cols, function(col) {
                x <- object[[col]][idx]
                na_omit <- !is.na(x)
                diag_rhat(x[na_omit], object$chain[idx][na_omit])
            }, numeric(1))
            names(rhats) <- addl_cols
            cat("R-hat values for summary statistics:\n")
            rhats_p <- vapply(rhats, function(x){
                ifelse(x < 1.05, sprintf('%.3f', x), paste0('\U274C', round(x, 3)))
            }, FUN.VALUE = character(1))
            print(noquote(rhats_p))

            if (any(na.omit(rhats) >= 1.05)) {
                warn_converge <- TRUE
                cli::cli_alert_danger("{.strong WARNING:} SMC runs have not converged.")
            }
            cat("\n")

        }

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
    } else if(algo %in% c("gsmc", 'basic_smc')) {
        pop_lb <- attr(object, "pop_bounds")[1]
        pop_ub <- attr(object, "pop_bounds")[3]

        cli_text("{.strong {algo}:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units with a population between {fmt_comma(pop_lb)} and {fmt_comma(pop_ub)}")
        cli_text("{.arg pop_temper}={format(all_diagn[[1]]$pop_temper, digits=3)}")
        cat("\n")

        cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
        if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")
        cat("\n")

        cols <- names(object)
        addl_cols <- setdiff(cols, c("chain", "draw", "district", "total_pop"))
        warn_converge <- FALSE
        if ("chain" %in% cols && length(addl_cols) > 0) {
            idx <- seq_len(n_samp)
            if ("district" %in% cols) idx <- as.integer(district) + (idx - 1)*n_distr

            const_cols <- vapply(addl_cols, function(col) {
                x <- object[[col]][idx]
                all(is.na(x)) || all(x == x[1]) ||
                    any(tapply(x, object[['chain']][idx], FUN = function(z) length(unique(z))) == 1)
            }, numeric(1))
            addl_cols <- addl_cols[!const_cols]

            rhats <- vapply(addl_cols, function(col) {
                x <- object[[col]][idx]
                na_omit <- !is.na(x)
                diag_rhat(x[na_omit], object$chain[idx][na_omit])
            }, numeric(1))
            names(rhats) <- addl_cols
            cat("R-hat values for summary statistics:\n")
            rhats_p <- vapply(rhats, function(x){
                ifelse(x < 1.05, sprintf('%.3f', x), paste0('\U274C', round(x, 3)))
            }, FUN.VALUE = character(1))
            print(noquote(rhats_p))

            if (any(na.omit(rhats) >= 1.05)) {
                warn_converge <- TRUE
                cli::cli_alert_danger("{.strong WARNING:} gSMC runs have not converged.")
            }
            cat("\n")

        }

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
                                   est_k = c(diagn$est_k, NA),
                                   unique_original = c(diagn$nunique_original_ancestors, NA))

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

            #

            names(tbl_print) <- c("Eff. samples (%)", "Acc. rate",
                                  "Log wgt. sd", " Max. unique",
                                  "k", "Unique Original Ancestors", "")
            rownames(tbl_print) <- c(paste("Split", seq_len(nrow(tbl_print) - 1)), "Resample")

            if (i == 1 || isTRUE(all_runs)) {
                cli_text("Sampling diagnostics for gSMC run {i} of {n_runs} ({fmt_comma(n_samp)} samples)")
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
            cli::cli_li("{.strong gSMC convergence:} Increase the number of samples.
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

        if (algo == 'mergesplit') {
            cli_text("{.strong Merge-Split MCMC:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")
        } else {
            cli_text("{.strong Flip MCMC:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")
        }


        accept_rate <- sprintf("%0.1f%%", 100*attr(object, "mh_acceptance"))
        cli_text("Chain acceptance rate{?s}: {accept_rate}")

        cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
        if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")
        cat("\n")

        cols <- names(object)
        addl_cols <- setdiff(cols, c("chain", "draw", "district", "total_pop", "mcmc_accept"))
        warn_converge <- FALSE
        if (length(addl_cols) > 0 && "chain" %in% cols) {
            idx <- seq_len(n_samp)
            if ("district" %in% cols) {
                idx <- 1 + (idx - 1)*n_distr
            }
            chain <- object$chain

            const_cols <- vapply(addl_cols, function(col) {
                x <- object[[col]][idx]
                all(is.na(x)) || all(x == x[1]) ||
                    any(tapply(x, object[['chain']][idx], FUN = function(z) length(unique(z))) == 1)
            }, logical(1))
            addl_cols <- addl_cols[!const_cols]

            rhats <- vapply(addl_cols, function(col) {
                x <- object[[col]][idx]
                na_omit <- !is.na(x)
                diag_rhat(x[na_omit], chain[idx][na_omit], split = TRUE)
            }, numeric(1))
            names(rhats) <- addl_cols
            cat("R-hat values for summary statistics:\n")
            rhats_p <- vapply(rhats, function(x){
                ifelse(x < 1.05, sprintf('%.3f', x), paste0('\U274C', round(x, 3)))
            }, FUN.VALUE = character(1))
            print(noquote(rhats_p))

            out <- tibble(stat = addl_cols, rhat = rhats)

            if (any(rhats >= 1.05)) {
                warn_converge <- TRUE
                cli::cli_alert_danger("{.strong WARNING:} Chains have not converged.")
            }
            cat("\n")
        } else {
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
        cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
    }

    invisible(out)
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
