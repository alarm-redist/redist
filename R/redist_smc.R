#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: tidy R wrapper to run SMC redistricting code
####################################################

#' SMC Redistricting Sampler (McCartan and Imai 2020)
#'
#' `redist_smc` uses a Sequential Monte Carlo algorithm (McCartan and Imai 2020)
#' to generate nearly independent congressional or legislative redistricting
#' plans according to contiguity, population, compactness, and administrative
#' boundary constraints.
#'
#' This function draws nearly-independent samples from a specific target measure,
#' controlled by the `map`, `compactness`, and `constraints` parameters.
#'
#' Key to ensuring good performance is monitoring the efficiency of the resampling
#' process at each SMC stage.  Unless `silent=FALSE`, this function will print
#' out the effective sample size of each resampling step to allow the user to
#' monitor the efficiency.  If `verbose=TRUE` the function will also print
#' out information on the \eqn{k_i} values automatically chosen and the
#' acceptance rate (based on the population constraint) at each step.
#' Users should also check diagnostics of the sample by running
#' \code{summary.redist_plans()}.
#'
#' Higher values of `compactness` sample more compact districts;
#' setting this parameter to 1 is computationally efficient and generates nicely
#' compact districts.  Values of other than 1 may lead to highly variable
#' importance sampling weights.  In these cases, these weights are by default
#' truncated using [redist_quantile_trunc()] to stabilize the resulting
#' estimates, but if truncation is used, a specific truncation function should
#' probably be chosen by the user.
#'
#' @param map A [redist_map()] object.
#' @param nsims The number of samples to draw.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided,
#' the algorithm will only generate maps which split up to `ndists-1`
#' counties. Even there are fewer counties than `ndists - 1`, the spanning
#' trees will change the results of the simulations. There is no strength
#' parameter associated with this constraint. To adjust the number of county
#' splits further, or to constrain a second type of administrative split,
#' consider using `add_constr_splits()`, `add_constr_multisplits()`, and
#' `add_constr_total_splits()`.
#' @param compactness Controls the compactness of the generated districts, with
#' higher values preferring more compact districts. Must be nonnegative. See
#' the 'Details' section for more information, and computational
#' considerations.
#' @param constraints A [redist_constr()] object or a list containing
#' information on sampling constraints. See [constraints] for more information.
#' @param resample Whether to perform a final resampling step so that the
#' generated plans can be used immediately.  Set this to `FALSE` to
#' perform direct importance sampling estimates, or to adjust the weights
#' manually.
#' @param runs How many independent parallel runs to conduct. Each run will
#' have `nsims` simulations. Multiple runs allows for estimation of simulation
#' standard errors. Output will only be shown for the first run. For
#' compatibility with MCMC methods, runs are identified with the `chain`
#' column in the output.
#' @param ncores How many cores to use to parallelize plan generation within each
#' run. The default, 0, will use the number of available cores on the machine
#' as long as `nsims` and the number of units is large enough. If `runs>1`
#' you will need to set this manually. If more than one core is used, the
#' sampler output will not be fully reproducible with `set.seed()`. If full
#' reproducibility is desired, set `ncores=1`.
#' @param init_particles A matrix of partial plans to begin sampling from. For
#' advanced use only.  The matrix must have `nsims` columns and a row for
#' every precinct. It is important to ensure that the existing districts meet
#' contiguity and population constraints, or there may be major issues when
#' sampling.
#' @param n_steps How many steps to run the SMC algorithm for.
#' Each step splits off a new district. Defaults to all remaining districts.
#' If fewer than the number of remaining splits, reference plans are disabled.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value `k_i` for each splitting iteration. Higher values are more accurate
#' but may require more computation. Set to 1 for the most conservative
#' sampling. Must be between 0 and 1.
#' @param seq_alpha The amount to adjust the weights by at each resampling step;
#' higher values prefer exploitation, while lower values prefer exploration.
#' Must be between 0 and 1.
#' @param truncate Whether to truncate the importance sampling weights at the
#' final step by `trunc_fn`.  Recommended if `compactness` is not 1.
#' Truncation only applied if `resample=TRUE`.
#' @param trunc_fn A function which takes in a vector of weights and returns a
#' truncated vector. If the [loo][loo::loo] package is installed (strongly
#' recommended), will default to Pareto-smoothed Importance Sampling (PSIS)
#' rather than naive truncation.
#' @param pop_temper The strength of the automatic population tempering. Try
#' values of 0.01-0.05 to start if the algorithm gets stuck on the final few
#' splits.
#' @param final_infl A multiplier for the population constraint on the final
#' iteration. Used to loosen the constraint when the sampler is getting stuck
#' on the final split. `pop_temper` should be tried first, since using
#' `final_infl` will actually change the target distribution.
#' @param est_label_mult A multiplier for the number of importance samples to
#' use in estimating the number of ways to sequentially label the districts.
#' Lower values increase speed at the cost of accuracy.  Only applied when
#' there are more than 13 districts.
#' @param ref_name a name for the existing plan, which will be added as a
#' reference plan, or `FALSE` to not include the initial plan in the
#' output. Defaults to the column name of the existing plan.
#' @param verbose Whether to print out intermediate information while sampling.
#' Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @return `redist_smc` returns a [redist_plans] object containing the simulated
#' plans.
#'
#' @references
#' McCartan, C., & Imai, K. (Forthcoming). Sequential Monte Carlo for Sampling
#' Balanced and Compact Redistricting Plans. *Annals of Applied Statistics*.
#' Available at \url{https://arxiv.org/abs/2008.06131}.
#'
#' @examples \donttest{
#' data(fl25)
#'
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#'
#' sampled_basic <- redist_smc(fl_map, 5000)
#'
#' constr <- redist_constr(fl_map)
#' constr <- add_constr_incumbency(constr, strength = 100, incumbents = c(3, 6, 25))
#' sampled_constr <- redist_smc(fl_map, 5000, constraints = constr)
#'
#' # Multiple parallel independent runs
#' redist_smc(fl_map, 1000, runs = 2)
#'
#' # One run with multiple cores
#' redist_smc(fl_map, 1000, ncores = 2)
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_smc <- function(map, nsims, counties = NULL, compactness = 1, constraints = list(),
                       resample = TRUE, runs = 1L, ncores = 0L, init_particles = NULL,
                       n_steps = NULL, adapt_k_thresh = 0.985, seq_alpha = 0.5,
                       truncate = (compactness != 1), trunc_fn = redist_quantile_trunc,
                       pop_temper = 0, final_infl = 1, est_label_mult = 1,
                       ref_name = NULL, verbose = FALSE, silent = FALSE) {
    map <- validate_redist_map(map)
    V <- nrow(map)
    adj <- get_adj(map)

    if (compactness < 0)
        cli_abort("{.arg compactness} must be non-negative.")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")
    if (seq_alpha <= 0 | seq_alpha > 1)
        cli_abort("{.arg seq_alpha} must lie in (0, 1].")
    if (nsims < 1)
        cli_abort("{.arg nsims} must be positive.")

    counties <- rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(counties)) {
        counties <- rep(1, V)
    } else {
        if (any(is.na(counties)))
            cli_abort("County vector must not contain missing values.")

        # handle discontinuous counties
        component <- contiguity(adj, vctrs::vec_group_id(counties))
        counties <- dplyr::if_else(component > 1,
                                   paste0(as.character(counties), "-", component),
                                   as.character(counties)) %>%
            as.factor() %>%
            as.integer()
        if (any(component > 1)) {
            cli_warn("Counties were not contiguous; expect additional splits.")
        }
    }

    # Other constraints
    if (!inherits(constraints, "redist_constr")) {
        constraints <- new_redist_constr(eval_tidy(enquo(constraints), map))
    }
    if (any(c("edges_removed", "log_st") %in% names(constraints))) {
        cli_warn(c("{.var edges_removed} or {.var log_st} constraint found in
           {.arg constraints} and will be ignored.",
            ">" = "Adjust using {.arg compactness} instead."))
    }
    if (any(c("poslby", "fry_hold") %in% names(constraints)) && compactness == 1) {
        cli_warn("{.var polsby} or {.var fry_hold} constraint found in {.arg constraints}
                 with {.arg compactness == 1). This may disrupt efficient sampling.")
    }
    constraints <- as.list(constraints) # drop data attribute

    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0

    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]
    ndists <- attr(map, "ndists")
    if (any(pop >= pop_bounds[3])) {
        too_big <- as.character(which(pop >= pop_bounds[3]))
        cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                    population larger than the district target.",
            "x" = "Redistricting impossible."))
    }

    # handle particle inits
    if (is.null(init_particles)) {
        init_particles <- matrix(0L, nrow = V, ncol = nsims)
        n_drawn <- 0L
    } else {
        if (nrow(init_particles) != V)
            cli_abort("{.arg init_particles} must have as many rows as {.arg map} has precincts.")
        if (ncol(init_particles) != nsims)
            cli_abort("{.arg init_particles} must have {.arg nsims} columns.")
        n_drawn <- as.integer(max(init_particles[, 1]))
    }
    if (is.null(n_steps)) {
        n_steps <- ndists - n_drawn - 1L
    }
    final_dists <- n_drawn + n_steps + 1L
    if (final_dists > ndists) {
        cli_abort("Too many districts already drawn to take {n_steps} steps.")
    }

    # set up parallel
    ncores_max <- parallel::detectCores()
    ncores_runs <- min(ncores_max, runs)
    ncores_per <- as.integer(ncores)
    if (ncores_per == 0) {
        if (nsims/100*length(adj)/200 < 20) {
            ncores_per <- 1L
        } else {
            ncores_per <- floor(ncores_max/ncores_runs)
        }
    }

    lags <- 1 + unique(round((ndists - 1)^0.8*seq(0, 0.7, length.out = 4)^0.9))
    control <- list(adapt_k_thresh = adapt_k_thresh,
                    seq_alpha = seq_alpha,
                    est_label_mult = est_label_mult,
                    adjust_labels = isTRUE(getOption("redist.adjust_labels", TRUE)),
                    pop_temper = pop_temper,
                    final_infl = final_infl,
                    lags = lags,
                    cores = as.integer(ncores_per))


    if (ncores_runs > 1) {
        `%oper%` <- `%dorng%`
        of <- if (Sys.info()[["sysname"]] == "Windows") {
            tempfile(pattern = paste0("smc_", substr(Sys.time(), 1, 10)), fileext = ".txt")
        } else {
            ""
        }

        if (!silent)
            cl <- makeCluster(ncores_runs, outfile = of, methods = FALSE,
                useXDR = .Platform$endian != "little")
        else
            cl <- makeCluster(ncores_runs, methods = FALSE,
                useXDR = .Platform$endian != "little")
        doParallel::registerDoParallel(cl, cores = ncores_runs)
        on.exit(stopCluster(cl))
    } else {
        `%oper%` <- `%do%`
    }

    t1 <- Sys.time()
    all_out <- foreach(chain = seq_len(runs), .inorder = FALSE, .packages="redist") %oper% {
        run_verbosity <- if (chain == 1) verbosity else 0
        t1_run <- Sys.time()

        algout <- smc_plans(nsims, adj, counties, pop, ndists,
                            pop_bounds[2], pop_bounds[1], pop_bounds[3],
                            compactness, init_particles, n_drawn, n_steps,
                            constraints, control, run_verbosity)
        # handle interrupt
        if (length(algout) == 0) {
            cli::cli_process_done()
            cli::cli_process_done()
        }

        lr <- -algout$lp
        wgt <- exp(lr - mean(lr))
        n_eff <- length(wgt)*mean(wgt)^2/mean(wgt^2)
        if (any(is.na(lr))) {
            cli_abort(c("Sampling probabilities have been corrupted.",
                "*" = "Check that none of your constraint weights are too large.
                             The output of constraint functions multiplied by the weight
                             should generally fall in the -5 to 5 range.",
                "*" = "If you are using custom constraints, make sure that your
                             constraint function handles all edge cases and never returns
                             {.val {NA}} or {.val {Inf}}",
                "*" = "If you are not using any constraints, please call
                             {.code rlang::trace_back()} and file an issue at
                             {.url https://github.com/alarm-redist/redist/issues/new}"))
        }

        n_unique <- NA
        if (resample) {
            if (!truncate) {
                mod_wgt <- wgt
            } else if (requireNamespace("loo", quietly = TRUE) && is.null(trunc_fn)) {
                mod_wgt <- wgt/sum(wgt)
                mod_wgt <- loo::weights.importance_sampling(
                    loo::psis(log(mod_wgt), r_eff = NA), log = FALSE)
            } else {
                mod_wgt <- trunc_fn(wgt)
            }
            mod_wgt <- wgt/sum(wgt)
            n_eff <- 1/sum(mod_wgt^2)

            # rs_idx <- sample(nsims, nsims, replace = TRUE, prob = mod_wgt)
            rs_idx <- resample_lowvar(mod_wgt)
            n_unique <- dplyr::n_distinct(rs_idx)
            algout$plans <- algout$plans[, rs_idx, drop = FALSE]
            # algout$log_labels = algout$log_labels[rs_idx]
            algout$ancestors <- algout$ancestors[rs_idx, , drop = FALSE]
            storage.mode(algout$ancestors) <- "integer"
        }
        storage.mode(algout$plans) <- "integer"
        t2_run <- Sys.time()

        if (!is.nan(n_eff) && n_eff/nsims <= 0.05)
            cli_warn(c("Less than 5% resampling efficiency.",
                       "*" = "Increase the number of samples.",
                       "*" = "Consider weakening or removing constraints.",
                       "i" = "If sampling efficiency drops precipitously in the final
                            iterations, population balance is likely causing a bottleneck.
                            Try increasing {.arg pop_temper} by 0.01.",
                       "i" = "If sampling efficiency declines steadily across iterations,
                            adjusting {.arg seq_alpha} upward may help a bit."))

        algout$wgt <- wgt

        algout$l_diag <- list(
            n_eff = n_eff,
            step_n_eff = algout$step_n_eff,
            adapt_k_thresh = adapt_k_thresh,
            est_label_mult = est_label_mult,
            est_k = algout$est_k,
            accept_rate = algout$accept_rate,
            sd_labels = algout$sd_labels,
            sd_lp = c(algout$sd_lp, sd(lr)),
            sd_temper = algout$sd_temper,
            cor_labels = algout$cor_labels,
            # log_labels = algout$log_labels,
            unique_survive = c(algout$unique_survive, n_unique),
            ancestors = algout$ancestors,
            seq_alpha = seq_alpha,
            pop_temper = pop_temper,
            runtime = as.numeric(t2_run - t1_run, units = "secs")
        )

        algout
    }
    if (verbosity >= 2) {
        t2 <- Sys.time()
        cli_text("{format(nsims*runs, big.mark=',')} plans sampled in
                 {format(t2-t1, digits=2)}")
    }

    plans <- do.call(cbind, lapply(all_out, function(x) x$plans))
    wgt <- do.call(c, lapply(all_out, function(x) x$wgt))
    l_diag <- lapply(all_out, function(x) x$l_diag)
    n_dist_act <- dplyr::n_distinct(plans[, 1]) # actual number (for partial plans)

    # tempering warning
    temp_ratio = do.call(c, lapply(l_diag, function(x) x$sd_temper / head(x$sd_lp, -1)))
    if (any(temp_ratio > 0.5, na.rm=TRUE)) {
        cli_warn(c("Population tempering is increasing the variance of the
                   resampling weights by over 50% at some steps.",
                   "*" = "Consider lowering {.arg pop_temper}."))
    }

    out <- new_redist_plans(plans, map, "smc", wgt, resample,
                            ndists = final_dists,
                            n_eff = all_out[[1]]$n_eff,
                            compactness = compactness,
                            constraints = constraints,
                            diagnostics = l_diag)
    if (runs > 1) {
        out <- mutate(out, chain = rep(seq_len(runs), each = n_dist_act*nsims)) %>%
            dplyr::relocate('chain', .after = "draw")
    }

    exist_name <- attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name) && ndists == final_dists) {
        ref_name <- if (!is.null(ref_name)) ref_name else exist_name
        out <- add_reference(out, map[[exist_name]], ref_name)
    }

    out
}


#' Helper function to truncate importance weights
#'
#' Defined as \code{pmin(x, quantile(x, 1 - length(x)^(-0.5)))}
#'
#' @param x the weights
#'
#' @return numeric vector
#'
#' @export
#'
#' @examples
#' redist_quantile_trunc(c(1, 2, 3, 4))
#'
redist_quantile_trunc <- function(x) pmin(x, quantile(x, 1 - length(x)^(-0.5)))
