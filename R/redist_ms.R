#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: tidy R wrapper to run Merge-Split/Recom redistricting code
####################################################

#' Merge-Split/Recombination MCMC Redistricting Sampler (Carter et al. 2019)
#'
#' \code{redist_mergesplit} uses a Markov Chain Monte Carlo algorithm (Carter et
#' al. 2019; based on DeFord et. al 2019) to generate congressional or legislative redistricting plans
#' according to contiguity, population, compactness, and administrative boundary
#' constraints. The MCMC proposal is the same as is used in the SMC sampler
#' (McCartan and Imai 2023); it is similar but not identical to those used in
#' the references.  1-level hierarchical Merge-split is supported through the
#' \code{counties} parameter; unlike in the SMC algorithm, this does not
#' guarantee a maximum number of county splits.
#'
#' This function draws samples from a specific target measure, controlled by the
#' \code{map}, \code{compactness}, and \code{constraints} parameters.
#'
#' Key to ensuring good performance is monitoring the acceptance rate, which
#' is reported at the sample level in the output.
#' Users should also check diagnostics of the sample by running
#' \code{summary.redist_plans()}.
#'
#' Higher values of \code{compactness} sample more compact districts;
#' setting this parameter to 1 is computationally efficient and generates nicely
#' compact districts.
#'
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to draw, including warmup.
#' @param warmup The number of warmup samples to discard. Recommended to be at
#' least the first 20% of samples, and in any case no less than around 100
#' samples, unless initializing from a random plan.
#' @param thin Save every `thin`-th sample. Defaults to no thinning (1).
#' @param init_plan The initial state of the map, provided as a single vector
#' to be shared across all chains, or a matrix with `chains` columns.
#' If not provided, will default to the reference map of the \code{map} object,
#' or if none exists, will sample a random initial state using
#' \code{\link{redist_smc}}. You can also request a random initial state by
#' setting \code{init_plan="sample"}.
#' @param chains The number of chains to run. Each chain will have `nsims`
#' draws. If `init_plan` is sampled, each chain will be initialized with its
#' own sampled plan. Defaults to 1.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided,
#' the algorithm will generate maps tend to follow county lines. There is no
#' strength parameter associated with this constraint. To adjust the number of
#' county splits further, or to constrain a second type of administrative
#' split, consider using `add_constr_splits()`, `add_constr_multisplits()`,
#' and `add_constr_total_splits()`.
#' @param compactness Controls the compactness of the generated districts, with
#' higher values preferring more compact districts. Must be nonnegative. See the
#' 'Details' section for more information, and computational considerations.
#' @param constraints A list containing information on constraints to implement.
#' See the 'Details' section for more information.
#' @param constraint_fn A function which takes in a matrix where each column is
#' a redistricting plan and outputs a vector of log-weights, which will be
#' added the the final weights.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
#' the algorithm does not appear to be sampling from the target distribution.
#' Must be between 0 and 1.
#' @param k The number of edges to consider cutting after drawing a spanning
#' tree. Should be selected automatically in nearly all cases.
#' @param ncores The number of parallel processes to run. Defaults to the
#' number of available cores, capped at the number of chains. Only used when
#' `chains > 1`.
#' @param cl_type The cluster type (see [parallel::makeCluster()]). Safest is
#' `"PSOCK"`, but `"FORK"` may be appropriate in some settings.
#' @param return_all If `TRUE` return all sampled plans; otherwise, just return
#' the final plan from each chain.
#' @param init_name a name for the initial plan, or \code{FALSE} to not include
#' the initial plan in the output.  Defaults to the column name of the
#' existing plan, or "\code{<init>}" if the initial plan is sampled.
#' @param silly_adj_fix Heuristic for fixing weird inputs.
#' @param verbose Whether to print out intermediate information while sampling.
#' Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @return \code{redist_mergesplit} returns an object of class
#' \code{\link{redist_plans}} containing the simulated plans. When
#' `chains > 1`, an additional `chain` column is included indicating which
#' chain each plan was drawn from.
#'
#' @references
#' Carter, D., Herschlag, G., Hunter, Z., and Mattingly, J. (2019). A
#' merge-split proposal for reversible Monte Carlo Markov chain sampling of
#' redistricting plans. arXiv preprint arXiv:1911.01503.
#'
#' McCartan, C., & Imai, K. (2023). Sequential Monte Carlo for Sampling
#' Balanced and Compact Redistricting Plans. *Annals of Applied Statistics* 17(4).
#' Available at \doi{10.1214/23-AOAS1763}.
#'
#' DeFord, D., Duchin, M., and Solomon, J. (2019). Recombination: A family of
#' Markov chains for redistricting. arXiv preprint arXiv:1911.05725.
#'
#' @examples \donttest{
#' data(fl25)
#'
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#'
#' sampled_basic <- redist_mergesplit(fl_map, 10000)
#'
#' sampled_constr <- redist_mergesplit(fl_map, 10000, constraints = list(
#'     incumbency = list(strength = 1000, incumbents = c(3, 6, 25))
#' ))
#'
#' sampled_parallel <- redist_mergesplit(fl_map, nsims = 100, chains = 2)
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_mergesplit <- function(map, nsims,
                              warmup = if (is.null(init_plan)) 10 else max(100, nsims %/% 5),
                              thin = 1L, init_plan = NULL, chains = 1L,
                              counties = NULL, compactness = 1,
                              constraints = list(), constraint_fn = function(m) rep(0, ncol(m)),
                              adapt_k_thresh = 0.99, k = NULL,
                              ncores = NULL, cl_type = "PSOCK", return_all = TRUE,
                              init_name = NULL, silly_adj_fix = FALSE,
                              verbose = FALSE, silent = FALSE) {
    if (!missing(constraint_fn)) cli::cli_warn("{.arg constraint_fn} is deprecated.")

    map <- validate_redist_map(map)
    V <- nrow(map)
    adj <- get_adj(map)
    ndists <- attr(map, "ndists")
    warmup <- max(warmup, 0L)
    thin <- as.integer(thin)
    chains <- as.integer(chains)

    if (compactness < 0)
        cli::cli_abort("{.arg compactness} must be non-negative.")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        cli::cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")
    if (nsims <= warmup)
        cli::cli_abort("{.arg nsims} must be greater than {.arg warmup}.")
    if (thin < 1 || thin > nsims - warmup)
        cli::cli_abort("{.arg thin} must be a positive integer, and no larger than {.arg nsims - warmup}.")
    if (nsims < 1)
        cli::cli_abort("{.arg nsims} must be positive.")
    if (chains < 1)
        cli::cli_abort("{.arg chains} must be positive.")

    exist_name <- attr(map, "existing_col")
    counties <- rlang::eval_tidy(rlang::enquo(counties), map)

    # Build init_plans matrix (V x chains)
    if (is.null(init_plan)) {
        if (!is.null(exist_name)) {
            init_plans <- matrix(rep(vctrs::vec_group_id(get_existing(map)), chains), ncol = chains)
            if (is.null(init_name)) {
                init_names <- rep(exist_name, chains)
            } else {
                init_names <- rep(init_name, chains)
            }
        } else {
            init_plan <- "sample"
        }
    }

    if (isTRUE(init_plan == "sample")) {
        if (!silent) cli::cli_inform("Sampling initial plans with SMC")
        init_plans <- get_plans_matrix(
            redist_smc(map, chains, counties, compactness, constraints,
                       resample = TRUE, adapt_k_thresh = adapt_k_thresh,
                       ref_name = FALSE, verbose = verbose, silent = silent, ncores = 1))
        if (is.null(init_name)) {
            init_names <- if (chains == 1) "<init>" else paste0("<init> ", seq_len(chains))
        } else {
            init_names <- paste(init_name, seq_len(chains))
        }
    } else if (!is.null(init_plan) && !is.character(init_plan)) {
        if (is.matrix(init_plan)) {
            if (ncol(init_plan) < chains)
                cli::cli_abort("{.arg init_plan} matrix must have at least {chains} column{?s}.")
            init_plans <- init_plan[, seq_len(chains), drop = FALSE]
        } else {
            init_plans <- matrix(rep(as.integer(init_plan), chains), ncol = chains)
        }
        if (is.null(init_name)) {
            init_names <- if (chains == 1) "<init>" else paste0("<init> ", seq_len(chains))
        } else {
            init_names <- rep(init_name, chains)
        }
    }

    # check init
    if (nrow(init_plans) != V)
        cli::cli_abort("{.arg init_plan} must be as long as the number of units as `map`.")
    if (max(init_plans) != ndists)
        cli::cli_abort("{.arg init_plan} must have the same number of districts as `map`.")
    if (any(apply(init_plans, 2, function(x) contiguity(adj, x)) != 1))
        cli::cli_warn("{.arg init_plan} should have contiguous districts.")

    if (is.null(counties)) {
        counties <- rep(1, V)
    } else {
        if (any(is.na(counties)))
            cli::cli_abort("County vector must not contain missing values.")

        if (silly_adj_fix) {
            for (j in seq_len(ndists)) {
                idx_distr <- which(init_plans[, 1] == j)
                adj_distr <- redist.reduce.adjacency(adj, idx_distr)
                component <- contiguity(adj_distr, vctrs::vec_group_id(counties[idx_distr]))
                counties[idx_distr] <- paste0(
                    j, ":",
                    dplyr::if_else(component > 1,
                                   paste0(as.character(counties[idx_distr]), "-", component),
                                   as.character(counties[idx_distr]))
                )
            }
            counties <- vctrs::vec_group_id(counties)
        } else {
            component <- contiguity(adj, vctrs::vec_group_id(counties))
            counties <- dplyr::if_else(component > 1,
                                       paste0(as.character(counties), "-", component),
                                       as.character(counties)) |>
                vctrs::vec_group_id()
        }
    }

    # Other constraints
    if (!inherits(constraints, "redist_constr")) {
        constraints <- new_redist_constr(eval_tidy(enquo(constraints), map))
    }
    if (any(c("edges_removed", "log_st") %in% names(constraints))) {
        cli::cli_warn(c("{.var edges_removed} or {.var log_st} constraint found in
           {.arg constraints} and will be ignored.",
            ">" = "Adjust using {.arg compactness} instead."))
    }
    if (any(c("poslby", "fry_hold") %in% names(constraints)) && compactness == 1) {
        cli::cli_warn("{.var polsby} or {.var fry_hold} constraint found in {.arg constraints}
                 with {.arg compactness != 1). This may disrupt efficient sampling.")
    }
    constraints <- as.list(constraints) # drop data attribute

    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0

    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]
    init_pop <- pop_tally(init_plans, pop, ndists)
    if (any(init_pop < pop_bounds[1]) | any(init_pop > pop_bounds[3]))
        cli::cli_abort("Provided initialization does not meet population bounds.")
    if (any(pop >= pop_bounds[3])) {
        too_big <- as.character(which(pop >= pop_bounds[3]))
        cli::cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                    population larger than the maximum district size.",
            "x" = "Redistricting impossible."))
    }

    control <- list(adapt_k_thresh = adapt_k_thresh, do_mh = TRUE)

    # Pre-estimate k (shared across all chains)
    if (is.null(k) || k <= 0) {
        k <- ms_plans(1, adj, init_plans[, 1], counties, pop, ndists, pop_bounds[2],
                      pop_bounds[1], pop_bounds[3], compactness,
                      list(), control, 0L, 1L, verbosity = 0)$est_k
    }
    k <- as.integer(k)

    if (chains == 1L) {
        # Single-chain path
        t1_run <- Sys.time()
        algout <- ms_plans(nsims, adj, init_plans[, 1], counties, pop, ndists,
                           pop_bounds[2], pop_bounds[1], pop_bounds[3], compactness,
                           constraints, control, k, thin, verbosity)
        t2_run <- Sys.time()

        storage.mode(algout$plans) <- "integer"
        acceptances <- as.logical(algout$mhdecisions)

        warmup_idx <- c(seq_len(1 + warmup %/% thin), ncol(algout$plans))
        l_diag <- list(runtime = as.numeric(t2_run - t1_run, units = "secs"))

        plans <- if (return_all) {
            algout$plans[, -warmup_idx, drop = FALSE]
        } else {
            algout$plans[, ncol(algout$plans) - 1L, drop = FALSE]
        }

        warmup_idx_acc <- c(seq_len(warmup %/% thin), length(acceptances))
        if (!return_all) {
            acceptances <- acceptances[length(acceptances) - 1L]
        } else {
            acceptances <- acceptances[-warmup_idx_acc]
        }

        out <- new_redist_plans(plans, map, "mergesplit", NULL, FALSE,
                                ndists = ndists,
                                compactness = compactness,
                                constraints = constraints,
                                adapt_k_thresh = adapt_k_thresh,
                                version = packageVersion("redist"),
                                diagnostics = list(l_diag),
                                mh_acceptance = mean(acceptances))
        out <- out %>% mutate(mcmc_accept = rep(acceptances, each = ndists))

        if (!is.null(init_names) && !isFALSE(init_name)) {
            out <- add_reference(out, init_plans[, 1], init_names[1])
        }
    } else {
        # Multi-chain parallel path
        if (is.null(ncores)) ncores <- parallel::detectCores()
        ncores <- min(ncores, chains)

        if (ncores > 1) {
            `%oper%` <- `%dorng%`
            of <- ifelse(Sys.info()[["sysname"]] == "Windows",
                         tempfile(pattern = paste0("ms_", substr(Sys.time(), 1, 10)), fileext = ".txt"),
                         "")
            if (!silent) {
                cl <- parallel::makeCluster(ncores, type = cl_type, outfile = of,
                                            methods = FALSE,
                                            useXDR = .Platform$endian != "little")
            } else {
                cl <- parallel::makeCluster(ncores, type = cl_type,
                                            methods = FALSE,
                                            useXDR = .Platform$endian != "little")
            }
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl))
        } else {
            `%oper%` <- `%do%`
        }

        out_par <- foreach::foreach(chain = seq_len(chains), .inorder = FALSE,
                                    .packages = "redist") %oper% {
            if (!silent) cat("Starting chain ", chain, "\n", sep = "")
            run_verbosity <- if (chain == 1 || verbosity == 3) verbosity else 0
            t1_run <- Sys.time()
            algout <- ms_plans(nsims, adj, init_plans[, chain], counties, pop,
                               ndists, pop_bounds[2], pop_bounds[1], pop_bounds[3],
                               compactness, constraints, control, k, thin, run_verbosity)
            t2_run <- Sys.time()

            algout$l_diag <- list(runtime = as.numeric(t2_run - t1_run, units = "secs"))
            algout$mh <- mean(as.logical(algout$mhdecisions))

            warmup_idx <- c(seq_len(1 + warmup %/% thin), nsims %/% thin + 2L)
            if (return_all) {
                algout$plans <- algout$plans[, -warmup_idx, drop = FALSE]
            } else {
                algout$plans <- algout$plans[, nsims + 1L, drop = FALSE]
            }
            storage.mode(algout$plans) <- "integer"

            warmup_idx_acc <- c(seq_len(warmup %/% thin), nsims %/% thin + 1L)
            if (!return_all) {
                algout$mhdecisions <- as.logical(algout$mhdecisions[nsims])
            } else {
                algout$mhdecisions <- as.logical(algout$mhdecisions[-warmup_idx_acc])
            }

            algout
        }

        plans <- lapply(out_par, function(algout) algout$plans)
        each_len <- ncol(plans[[1]])
        plans <- do.call(cbind, plans)
        storage.mode(plans) <- "integer"

        mh <- sapply(out_par, function(algout) algout$mh)
        l_diag <- lapply(out_par, function(algout) algout$l_diag)
        acceptances <- sapply(out_par, function(algout) algout$mhdecisions)

        out <- new_redist_plans(plans = plans, map = map, algorithm = "mergesplit",
                                wgt = NULL, resampled = FALSE,
                                compactness = compactness,
                                constraints = constraints,
                                ndists = ndists,
                                adapt_k_thresh = adapt_k_thresh,
                                mh_acceptance = mh,
                                version = packageVersion("redist"),
                                diagnostics = l_diag) %>%
            mutate(chain = rep(seq_len(chains), each = each_len * ndists),
                   mcmc_accept = rep(acceptances, each = ndists))

        if (!is.null(init_names) && !isFALSE(init_name)) {
            if (all(init_names[1] == init_names)) {
                out <- add_reference(out, init_plans[, 1], init_names[1])
            } else {
                out <- Reduce(function(cur, idx) {
                    add_reference(cur, init_plans[, idx], init_names[idx]) %>%
                        mutate(chain = dplyr::coalesce(chain, idx))
                }, rev(seq_len(chains)), init = out)
            }
        }

        out <- dplyr::relocate(out, chain, .after = "draw")
    }

    out
}
