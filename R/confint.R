#' Confidence Intervals for SMC and MCMC Estimates
#'
#' Builds a confidence interval for a quantity of interest.
#' If multiple runs are available, uses the between-run variation to estimate
#' the standard error. If only one run is available, uses information on the SMC
#' particle/plan genealogy to estimate the standard error, using a variant of
#' the method of Olson & Douc (2019). The multiple-run estimator is more
#' reliable, especially for situations with many districts, and should be used
#' when parallelism is available.  All reference plans are ignored.
#'
#' @param plans a [redist_plans] object.
#' @param x the quantity to build an interval for. Tidy-evaluated within `plans`.
#' @param district for [redist_plans] objects with multiple districts, which
#' `district` to subset to. Set to `NULL` to perform no subsetting.
#' @param conf the desired confidence level.
#' @param by_chain Whether the confidence interval should indicate overall
#'   sampling uncertainty (`FALSE`) or per-chain sampling uncertainty (`TRUE`).
#'   In the latter case the intervals will be wider by a factor of `sqrt(runs)`.
#'
#' @references
#' Lee, A., & Whiteley, N. (2018). Variance estimation in the particle filter.
#' Biometrika, 105(3), 609-625.
#'
#' Olsson, J., & Douc, R. (2019). Numerically stable online estimation of
#' variance in particle filters. Bernoulli, 25(2), 1504-1535.
#'
#' H. P. Chan and T. L. Lai. A general theory of particle filters in hidden
#' Markov models and some applications. Ann. Statist., 41(6):2877–2904, 2013.
#'
#' @return A tibble with three columns: \code{X}, \code{X_lower}, and
#' \code{X_upper}, where \code{X} is the name of the vector of interest,
#' containing the mean and confidence interval. When used inside
#' \code{\link[dplyr:summarise]{summarize()}} this will create three columns in the
#' output data.
#'
#' @examples
#' library(dplyr)
#' data(iowa)
#'
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
#' plans <- redist_mergesplit(iowa_map, nsims = 200, chains = 2, silent = TRUE) |>
#'     mutate(dem = group_frac(iowa_map, dem_08, dem_08 + rep_08)) |>
#'     number_by(dem)
#' redist_smc_ci(plans, dem)
#'
#' @md
#' @concept analyze
#' @export
redist_ci <- function(plans, x, district = 1L, conf = 0.9, by_chain = FALSE) {
    algo <- attr(plans, "algorithm")
    algos_ok <- c("smc", "smc_ms", "mergesplit", "cyclewalk", "flip")

    x <- enquo(x)

    if (is.null(algo) || !algo %in% algos_ok) {
        cli::cli_abort("{.field algorithm} attribute missing from {.arg plans}.
                  Call {.fn redist_smc_ci} or {.fn redist_mcmc_ci} directly.")
    } else if (algo %in% c("smc", "smc_ms")) {
        redist_smc_ci(plans, !!x, district, conf, by_chain)
    } else {
        # MCMC
        redist_mcmc_ci(plans, !!x, district, conf, by_chain)
    }
}

#' Get the ancestor matrix used by the SMC variance estimator
#'
#' Entry `[i, j]` is the index of the ancestor of particle `i`, traced back to
#' the generation corresponding to lag `j`. Particles sharing an ancestor form
#' the groups used by the Lee-Whiteley / Olsson-Douc estimator.
#'
#' This is reconstructed from `parent_index_mat`, which records the parent of
#' every particle at every SMC step. The C++ code used to build the equivalent
#' matrix directly via
#'
#'     ancestors[i, j] <- if (dist_ctr <= lags[j]) i else ancestors[parent, j]
#'
#' so composing `parent_index_mat` lookups backwards from the final step
#' reproduces it exactly. Older `redist_plans` objects that predate
#' `parent_index_mat` fall back to the stored `ancestors` matrix.
#'
#' @param plans a `redist_plans` object from `redist_smc()`
#' @returns an integer matrix with one row per particle and one column per lag
#' @noRd
get_smc_ancestors <- function(plans) {
    parent_index_mat <- attr(plans, "internal_diagnostics")[[1]]$parent_index_mat

    if (is.null(parent_index_mat)) {
        # Pre-parent_index_mat objects; use the ancestor matrix they stored.
        ancestors <- attr(plans, "diagnostics")[[1]]$ancestors
        if (is.null(ancestors)) {
            cli::cli_abort("{.arg plans} has neither {.field parent_index_mat} nor
                            {.field ancestors} diagnostics, so standard errors
                            cannot be computed.")
        }
        return(ancestors)
    }

    n <- nrow(parent_index_mat)
    n_cols <- ncol(parent_index_mat)
    ndists <- attr(plans, "ndists")

    # When plans were resampled an extra column holding the resampling indices
    # was appended, so only the remaining columns are SMC steps.
    resampled <- isTRUE(attr(plans, "resampled"))
    n_smc_steps <- n_cols - resampled
    # each SMC step splits one region in two, so the run started here
    init_num_regions <- ndists - n_smc_steps
    lags <- 1 + unique(round((ndists - 1)^0.8 * seq(0, 0.7, length.out = 4)^0.9))

    ancestors <- matrix(0L, nrow = n, ncol = length(lags))
    for (j in seq_along(lags)) {
        # `dist_ctr` at SMC step s (0-indexed) is init_num_regions + s, so the
        # recursion resets to the identity for every step with
        # init_num_regions + s <= lags[j]. Only the last such reset matters.
        reset_step <- lags[j] - init_num_regions
        # column s + 1 holds SMC step s, so compose from the last SMC column
        # back down to the first column after the reset.
        first_col <- max(1L, reset_step + 2L)

        anc <- seq_len(n)
        # Resampling happens after every SMC step, so it always applies,
        # even for lags that reset past the end of the run.
        if (resampled) {
            anc <- parent_index_mat[anc, n_cols]
        }
        if (first_col <= n_smc_steps) {
            for (s in seq(n_smc_steps, first_col)) {
                anc <- parent_index_mat[anc, s]
            }
        }
        ancestors[, j] <- anc
    }

    ancestors
}

#' @describeIn redist_ci Compute confidence intervals for SMC output.
#' @export
redist_smc_ci <- function(plans, x, district = 1L, conf = 0.9, by_chain = FALSE) {
    plans <- subset_sampled(plans)
    x_orig <- enquo(x)
    x <- as.numeric(eval_tidy(enquo(x), plans))
    if (!"district" %in% names(plans)) {
        plans$district <- rep(1, nrow(plans))
    }
    if (!is.null(district)) {
        x <- x[plans$district == district]
    }
    N <- length(x)
    est <- mean(x)

    if ("chain" %in% names(plans)) {
        # multiple runs
        if (is.null(district)) {
            chain <- plans$chain
        } else {
            chain <- plans$chain[plans$district == district]
        }
        rhat <- diag_rhat(x, chain)
        if (is.finite(rhat) && rhat > 1.05) {
            cli::cli_warn(c("Runs have not converged for this statistic.",
                       "i" = "R-hat is {round(rhat, 3)}",
                       ">" = "Increase the number of samples."))
        }
        run_means <- tapply(x, chain, mean) |> `names<-`(NULL)

        if (isTRUE(by_chain)) {
            std_err <- sd(run_means)
        } else {
            std_err <- sd(run_means) / sqrt(max(chain) - 1) # be slightly conservative
        }
    } else {
        # See
        # OLSSON, J. and DOUC, R. (2019). Numerically stable online estimation of variance in particle filters. Bernoulli
        # 25 1504–1535.
        # LEE, A. and WHITELEY, N. (2018). Variance estimation in the particle filter. Biometrika 105 609–625.
        ancestors <- get_smc_ancestors(plans)
        std_errs <- apply(ancestors, 2, function(anc) {
            sum_inner <- tapply(x - est, anc, sum)^2
            sqrt(mean(sum_inner[as.character(anc)]) / N)
        })
        std_err <- quantile(std_errs, 0.75)
    }

    alpha <- (1 - conf) / 2
    ci <- est + qt(c(alpha, 0.5, 1 - alpha), df = N - 1) * std_err

    tibble(
        "{{ x_orig }}" := ci[2],
        "{{ x_orig }}_lower" := ci[1],
        "{{ x_orig }}_upper" := ci[3]
    )
}

#' @describeIn redist_ci Compute confidence intervals for MCMC output.
#' @param use_coda Whether or not to use the coda package to compute standard errors
#' @export
redist_mcmc_ci <- function(
    plans,
    x,
    district = 1L,
    conf = 0.9,
    by_chain = FALSE,
    use_coda = FALSE
) {
    plans <- subset_sampled(plans)
    x_orig <- enquo(x)
    x <- as.numeric(eval_tidy(enquo(x), plans))
    if (!"district" %in% names(plans)) {
        plans$district <- rep(1, nrow(plans))
    }
    if (!is.null(district)) {
        x <- x[plans$district == district]
    }
    N <- length(x)
    est <- mean(x)
    thin <- attr(plans, "diagnostics")[[1]]$thin

    if ("chain" %in% names(plans)) {
        # multiple runs
        chain <- plans$chain[plans$district == district]
    } else {
        chain <- rep(1, N)
    }

    # Fall back to coda when only a single chain is available.
    if (!use_coda && length(unique(chain)) < 2) {
        rlang::check_installed("coda", "to calculate single-chain MCMC standard errors.")
        use_coda <- TRUE
    }

    rhat <- diag_rhat(x, chain, split = TRUE)
    if (is.finite(rhat) && rhat > 1.05) {
        cli::cli_warn(c("Runs have not converged for this statistic.",
                   "i" = "R-hat is {round(rhat, 3)}",
                   ">" = "Increase the number of samples."))
    }

    if (use_coda) {
        rlang::check_installed("coda", "to calculate MCMC standard errors.")
        mcmc <- coda::mcmc.list(tapply(x, chain, coda::mcmc, thin = thin))
        std_err <- summary(mcmc)$statistics["Time-series SE"]
        if (isTRUE(by_chain)) {
            std_err <- std_err * sqrt(max(chain))
        }
    } else if ("chain" %in% names(plans)) {
        # multiple runs
        if (is.null(district)) {
            chain <- plans$chain
        } else {
            chain <- plans$chain[plans$district == district]
        }
        rhat <- diag_rhat(x, chain)
        if (is.finite(rhat) && rhat > 1.05) {
            cli::cli_warn(c("Runs have not converged for this statistic.",
                           "i" = "R-hat is {round(rhat, 3)}",
                           ">" = "Increase the number of samples."))
        }
        run_means <- tapply(x, chain, mean) |> `names<-`(NULL)

        if (isTRUE(by_chain)) {
            std_err <- sd(run_means)
        } else {
            std_err <- sd(run_means) / sqrt(max(chain) - 1) # be slightly conservative
        }
    } else {
        cli::cli_abort("Can't do non-coda std errors for single MCMC chain!")
    }

    alpha <- (1 - conf) / 2
    ci <- est + qt(c(alpha, 0.5, 1 - alpha), df = N - 1) * std_err

    tibble(
        "{{ x_orig }}" := ci[2],
        "{{ x_orig }}_lower" := ci[1],
        "{{ x_orig }}_upper" := ci[3]
    )
}


#' (Deprecated) Confidence Intervals for Importance Sampling Estimates
#'
#' Builds a confidence interval for a quantity of interest,
#' given importance sampling weights.
#'
#' @param x A numeric vector containing the quantity of interest
#' @param wgt A numeric vector containing the nonnegative importance weights.
#' Will be normalized automatically.
#' @param conf The confidence level for the interval.
#'
#' @returns A two-element vector of the form `[lower, upper]` containing
#' the importance sampling confidence interval.
#'
#' @concept post
#' @export
redist.smc_is_ci <- function(x, wgt, conf = 0.99) {
    .Deprecated("redist_smc_ci")
    wgt <- wgt / sum(wgt)
    mu <- sum(x * wgt)
    sig <- sqrt(sum((x - mu)^2 * wgt^2))
    mu + qnorm(c((1 - conf)/2, 1 - (1 - conf)/2)) * sig
}
