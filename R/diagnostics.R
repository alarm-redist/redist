
#' Diagnostic information on sampled plans
#'
#' Prints diagnostic information, which varies by algorithm. All algorithms
#' compute the [plans_diversity()] of the samples.
#'
#' For SMC and MCMC, if there are multiple runs/chains, R-hat values will be
#' computed for each summary statistic. These values should be close to 1.
#' If they are not, then there is too much between-chain variation, indicating
#' that there are not enough samples. R-hat values are calculated after
#' rank-normalization and folding. For summary statistics that vary across
#' districts, R-hat is calculated for the first district only.
#'
#' For SMC, diagnostics statistics include:
#'
#' * **Effective samples**: the effective sample size at each iteration, computed
#' using the SMC weights. Larger is better. The percentage in parentheses is the
#' ratio of the effective samples to the total samples.
#' * **Acceptance rate**: the fractino of drawn spanning trees which yield a valid
#' redistricting plan within the population tolerance. Very small values (< 1%)
#' can indicate a bottleneck and may lead to a lack of diversity.
#' * **Standard deviation of the log weights**: More variable weights (larger s.d.)
#' indicate less efficient sampling. Values greater than 3 are likely problematic.
#' * **Maximum unique plans:** an upper bound on the number of unique redistricting
#' plans that survive each stage. The percentage in parentheses is the ratio of
#' this number to expected number of unique plans under equal-probability
#' multinomial resampling. Small values (< 100) indicate a bottleneck, which
#' leads to a loss of sample diversity and a higher variance.
#' * **Estimated `k` parameter**: How many spanning tree edges were considered for
#' cutting at each split. Mostly informational, though large jumps may indicate
#' a need to increase `adapt_k_thresh`.
#' * **Bottleneck**: Will show an asterisk if a bottleneck appears likely, based on
#' the values of the other statistics.
#'
#' In the event of problematic diagnostics, the function will provide
#' suggestions for improvement.
#'
#' @param object a [redist_plans] object
#' @param \dots additional arguments (ignored)
#'
#' @method summary redist_plans
#' @return A data frame containing diagnostic information, invisibly.
#'
#' @examples
#' data(iowa)
#' iowa_map = redist_map(iowa, ndists=4, pop_tol=0.1)
#' plans = redist_smc(iowa_map, 100)
#' summary(plans)
#'
#' @md
#' @export
summary.redist_plans = function(object, ...) {
    algo = attr(object, "algorithm")
    name = deparse(substitute(object))

    object = subset_sampled(object)
    diagn = attr(object, "diagnostics")
    plans_m = get_plans_matrix(object)
    n_samp = ncol(plans_m)
    n_distr = attr(object, "ndists")
    if (is.null(n_distr)) n_distr = max(plans_m[,1])

    fmt_comma = function(x) format(x, digits=0, big.mark=",")

    prec_pop = attr(object, "prec_pop")
    if (is.null(prec_pop)) {
        cli_warn(c("Precinct population missing; plan diversity estimates may be misleading.",
                   ">"='Run `attr({name}, "prec_pop") <- <map object>$<pop column>` to fix.'))
        prec_pop = rep(1, nrow(plans_m))
    }
    est_div = plans_diversity(object, total_pop=prec_pop)
    div_rg = format(quantile(est_div, c(0.1, 0.9)), digits=2)
    div_bad = (mean(est_div) <= 0.35) || (mean(est_div <= 0.05) > 0.2)

    if (algo == "smc") {
        cli_text("{.strong SMC:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")

        cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
        if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")
        cat("\n")

        cols = names(object)
        addl_cols = setdiff(cols, c("chain", "draw", "district", "total_pop"))
        warn_converge = FALSE
        if ("chain" %in% cols && length(addl_cols) > 0) {
            idx = seq_len(n_samp)
            if ("district" %in% cols) idx = 1 + (idx - 1) * n_distr

            rhats = vapply(addl_cols, function(col) {
                x = object[[col]][idx]
                na_omit = !is.na(x)
                diag_rhat(x[na_omit], object$chain[idx][na_omit])
            }, numeric(1))
            names(rhats) = addl_cols
            cat("R-hat values for summary statistics:\n")
            print(rhats)

            if (any(rhats >= 1.05)) {
                warn_converge = TRUE
                cli::cli_alert_danger("{.strong WARNING:} SMC runs have not converged.")
            }
            cat("\n")
        }


        out = tibble(n_eff = c(diagn$step_n_eff, diagn$n_eff),
                     eff = c(diagn$step_n_eff, diagn$n_eff)/n_samp,
                     accept_rate = c(diagn$accept_rate, NA),
                     sd_log_wgt = diagn$sd_lp,
                     max_unique = diagn$unique_survive,
                     est_k = c(diagn$est_k, NA))

        tbl_print = as.data.frame(out)
        min_n = max(0.05*n_samp, min(0.4*n_samp, 100))
        bottlenecks = dplyr::coalesce(with(tbl_print, pmin(max_unique, n_eff) < min_n), FALSE)
        tbl_print$bottleneck = ifelse(bottlenecks, "     *     ", "")
        tbl_print$n_eff = with(tbl_print,
                str_glue("{fmt_comma(n_eff)} ({sprintf('%0.1f%%', 100*eff)})"))
        tbl_print$eff = NULL
        tbl_print$accept_rate = with(tbl_print, sprintf('%0.1f%%', 100*accept_rate))
        max_pct = with(tbl_print, max_unique/(-n_samp * expm1(-1)))
        tbl_print$max_unique = with(tbl_print,
                str_glue("{fmt_comma(max_unique)} ({sprintf('%3.0f%%', 100*max_pct)})"))

        colnames(tbl_print) = c("Eff. samples (%)", "Acc. rate",
                                "Log wgt. sd", " Max. unique",
                                "Est. k", "Bottleneck?")
        rownames(tbl_print) = c(paste("Split", seq_len(n_distr-1)), "Resample")

        if ("chain" %in% cols) cli::cli_alert_info("Sampling diagnostics shown for first SMC run only.")
        print(tbl_print, digits=2)
        cat("\n")

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
                        the population tolerance. If the accpetance rate drops
                        quickly in the final splits, try increasing
                        {.arg pop_temper} by 0.01.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong SMC convergence:} Increase the number of samples.
                        If you are experiencing low plan diversity or bottlenecks as well,
                        address those issues first.")
        }
        if (any(bottlenecks)) {
            cli::cli_li("{.strong Bottlenecks:} Consider weakening or removing
                        constraints, or increasing the population tolerance.
                        If the accpetance rate drops quickly in the final splits,
                        try increasing {.arg pop_temper} by 0.01.
                        To visualize what geographic areas may be causing problems,
                        try running the following code. Highlighted areas are
                        those that may be causing the bottleneck.")
            code = str_glue("plot(<map object>, colMeans(as.matrix({name}) %in% c({paste(which(bottlenecks), collapse=', ')})))")
            cli::cat_line("    ", cli::code_highlight(code, "Material"))
        }
    } else if (algo == "mergesplit") {
        cli_text("{.strong Merge-Split MCMC:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")

        accept_rate = sprintf("%0.1f%%", 100*attr(object, "mh_acceptance"))
        cli_text("Chain acceptance rate{?s}: {accept_rate}")

        cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
        if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")
        cat("\n")

        cols = names(object)
        addl_cols = setdiff(cols, c("chain", "draw", "district", "total_pop", "mcmc_accept"))
        warn_converge = FALSE
        if ("chain" %in% cols && length(addl_cols) > 0) {
            idx = seq_len(n_samp)
            if ("district" %in% cols) idx = 1 + (idx - 1) * n_distr

            rhats = vapply(addl_cols, function(col) {
                x = object[[col]][idx]
                na_omit = !is.na(x)
                diag_rhat(x[na_omit], object$chain[idx][na_omit])
            }, numeric(1))
            names(rhats) = addl_cols
            cat("R-hat values for summary statistics:\n")
            print(rhats)

            out = tibble(stat=addl_cols, rhat=rhats)

            if (any(rhats >= 1.05)) {
                warn_converge = TRUE
                cli::cli_alert_danger("{.strong WARNING:} Chains have not converged.")
            }
            cat("\n")
        } else {
            out = NULL
        }

        cli::cli_li(cli::col_grey("
            Watch out for low acceptance rates (less than 10%).
            R-hat values for summary statistics should be between 1 and 1.05."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Increase the number of samples.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong Chain convergence:} Increase the number of samples.
                        If you are experiencing low plan diversity, address that issues first.")
        }
    } else {
        cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
    }

    invisible(out)
}


diag_fold = function(x) {
    abs(x - median(x))
}

diag_ranknorm = function(x) {
    qnorm(rank(x)/(length(x)+1))
}

diag_calc_rhat = function(x, grp) {
    n = mean(table(grp))
    var_between = n * var(tapply(x, grp, mean))
    var_within = mean(tapply(x, grp, var))
    sqrt((var_between / var_within + n - 1) / n)
}

diag_rhat = function(x, grp) {
    max(diag_calc_rhat(diag_ranknorm(x), grp),
        diag_calc_rhat(diag_ranknorm(diag_fold(x)), grp))
}
