#' Marked Edge Walk MCMC Redistricting Sampler (McWhorter and DeFord 2024)
#'
#' `redist_mew` uses a Markov Chain Monte Carlo algorithm based on marked
#' edge walks on spanning trees (McWhorter and DeFord 2024) to generate
#' congressional or legislative redistricting plans according to contiguity,
#' population, compactness, and other custom constraints.
#'
#' @param map A [redist_map] object.
#' @param nsims The number of samples to draw, including warmup.
#' @param warmup The number of warmup samples to discard. Recommended to be at
#'   least the first 20% of samples, and in any case no less than around 100
#'   samples, unless initializing from a random plan.
#' @param thin Save every `thin`-th sample. Defaults to no thinning (1).
#' @param init_plan The initial state of the map. If not provided, will default
#'   to the reference map of the `map` object, or if none exists, will
#'   sample a random initial state using [redist_smc()]. You can also
#'   request a random initial state by setting `init_plan="sample"`.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative. See
#'   the 'Details' section for more information.
#' @param constraints A [redist_constr] object or a list containing
#'   information on sampling constraints. See [redist_constr] for more
#'   information.
#' @param adapt_k_thresh The threshold value used in the heuristic to adapt the
#'   number of marked edges. Higher values are more conservative. Set to 0.9999
#'   or 1 if the algorithm does not appear to be sampling from the target
#'   distribution. Must be between 0 and 1.
#' @param warmstart If `TRUE`, the algorithm will initialize the spanning
#'   tree from the structure of `init_plan`. If `FALSE` (the
#'   default), uses Wilson's algorithm to sample a uniform random spanning tree,
#'   then finds marked edges that match `init_plan`.
#' @param init_name A name for the initial plan, or `FALSE` to not include
#'   the initial plan in the output. Defaults to the column name of the existing
#'   plan, or `"<init>"` if the initial plan is sampled.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @returns `redist_mew` returns a [redist_plans] object containing the
#'   simulated plans.
#'
#' @references
#' McWhorter, A., & DeFord, D. (2024). The Marked Edge Walk: A Novel MCMC
#' Algorithm for Sampling of Graph Partitions.
#'
#' @examples
#' data(fl25)
#'
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#'
#' # Basic sampling
#' sampled_basic <- redist_mew(fl_map, 5000)
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_mew <- function(map, nsims,
                       warmup = 0L,
                       thin = 1L,
                       init_plan = NULL,
                       compactness = 1,
                       constraints = list(),
                       adapt_k_thresh = 0.99,
                       warmstart = FALSE,
                       init_name = NULL,
                       verbose = FALSE,
                       silent = FALSE) {

    map <- validate_redist_map(map)
    V <- nrow(map)
    adj <- get_adj(map)
    ndists <- attr(map, "ndists")
    warmup <- max(warmup, 0L)
    thin <- as.integer(thin)

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

    # Handle init_plan
    exist_name <- attr(map, "existing_col")
    orig_lookup <- seq_len(ndists)

    if (is.null(init_plan) && !is.null(exist_name)) {
        init_plan <- vctrs::vec_group_id(get_existing(map))
        orig_lookup <- unique(get_existing(map))
        if (is.null(init_name)) init_name <- exist_name
    } else if (!is.null(init_plan) && is.null(init_name)) {
        init_name <- "<init>"
    }

    if (length(init_plan) == 0L || isTRUE(init_plan == "sample")) {
        if (!silent) cli::cli_alert_info("Sampling initial plan with SMC...")
        init_plan <- tryCatch({
            as.integer(get_plans_matrix(
                redist_smc(map, 10, resample = FALSE, ref_name = FALSE,
                          silent = TRUE, ncores = 1))[, 1])
        }, error = function(e) {
            cli::cli_abort(c(
                "Failed to generate initial plan using SMC.",
                "i" = "This can happen on very small or difficult maps.",
                "i" = "Try providing an {.arg init_plan} manually.",
                "x" = conditionMessage(e)
            ))
        })
        if (is.null(init_name)) init_name <- "<init>"
    }

    # Check init_plan validity
    if (length(init_plan) != V)
        cli::cli_abort("{.arg init_plan} must be as long as the number of units as `map`.")
    if (max(init_plan) != ndists)
        cli::cli_abort("{.arg init_plan} must have the same number of districts as `map`.")
    if (any(contiguity(adj, init_plan) != 1))
        cli::cli_abort("{.arg init_plan} must be contiguous.")

    # Convert init_plan to 1-indexed integers
    init_plan <- match(init_plan, sort(unique(init_plan)))

    # Process constraints
    if (inherits(constraints, "redist_constr")) {
        constraints <- as.list(constraints)
    }

    # Population parameters
    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]


    # Setup control parameters
    control <- list(
        adapt_k_thresh = adapt_k_thresh,
        warmstart = warmstart,
        verbosity = ifelse(silent, 0L, ifelse(verbose, 2L, 1L))
    )

    if (!silent) {
        cli::cli_rule(left = "Sampling using Marked Edge Walk")
        cli::cli_alert_info("{ndists} districts on {V} units")
        if (!is.null(pop_tol)) {
            cli::cli_alert_info("Target population: {round(target_pop)} ± {scales::percent(pop_tol)}")
        } else {
            cli::cli_alert_info("Target population: {round(target_pop)}")
        }
        if (!is.null(init_name) && init_name != FALSE) {
            cli::cli_alert_info("Initial plan: {init_name}")
        }
        if (warmstart) {
            cli::cli_alert_info("Using warmstart initialization")
        } else {
            cli::cli_alert_info("Initializing with Wilson's algorithm")
        }
    }

    # Call C++ function
    time_start <- proc.time()
    plans_raw <- mew_plans(
        nsims = nsims,
        adj = adj,
        init = init_plan,
        pop = pop,
        n_distr = ndists,
        target = pop_bounds[2],
        lower = pop_bounds[1],
        upper = pop_bounds[3],
        rho = compactness,
        constraints = constraints,
        control = control,
        thin = thin,
        verbosity = control$verbosity
    )
    time_elapsed <- (proc.time() - time_start)[["elapsed"]]

    # Convert plans to original district labels
    plans_m <- plans_raw$plans
    for (i in seq_len(ncol(plans_m))) {
        plans_m[, i] <- orig_lookup[plans_m[, i]]
    }

    # Remove warmup samples
    # After thinning in C++, we have nsims/thin total samples
    # We want to remove the first warmup/thin samples
    if (warmup > 0) {
        n_warmup_thinned <- ceiling(warmup / thin)
        if (n_warmup_thinned < ncol(plans_m)) {
            plans_m_final <- plans_m[, -(1:n_warmup_thinned), drop = FALSE]
        } else {
            # All samples are warmup - keep at least the last one
            plans_m_final <- plans_m[, ncol(plans_m), drop = FALSE]
        }
    } else {
        plans_m_final <- plans_m
    }

    # Build redist_plans object
    n_out <- ncol(plans_m_final)
    plans <- new_redist_plans(
        plans = plans_m_final,
        map = map,
        algorithm = "mew",
        wgt = rep(1, n_out)
    )

    # Add diagnostics (but trim to match number of non-warmup samples)
    # Diagnostics are computed every 100 iterations, so we need to subset correctly
    diag_every <- 100
    n_diag <- length(plans_raw$accept_rate)
    warmup_diag_idx <- seq_len(pmax(1, warmup %/% diag_every))

    # Subset diagnostics to remove warmup, but keep at least one value
    if (length(warmup_diag_idx) < n_diag) {
        diag_idx <- setdiff(seq_len(n_diag), warmup_diag_idx)
    } else {
        diag_idx <- n_diag # Keep only last diagnostic
    }

    attr(plans, "diagnostics") <- tibble::tibble(
        accept_rate = plans_raw$accept_rate[diag_idx],
        cycle_intersect_rate = plans_raw$cycle_intersect_rate[diag_idx],
        avg_proposal_tries = plans_raw$avg_proposal_tries[diag_idx]
    )

    if (!silent) {
        cli::cli_alert_success("Sampling complete in {round(time_elapsed, 1)}s")
        cli::cli_alert_info("Average acceptance rate: {scales::percent(mean(attr(plans, 'diagnostics')$accept_rate))}")
    }

    # Add reference plan if requested
    # NOTE: Skipping for now due to district renumbering issues
    # if (!is.null(init_name) && init_name != FALSE) {
    #     init_out <- matrix(orig_lookup[init_plan], ncol = 1)
    #     colnames(init_out) <- init_name
    #     plans <- add_reference(plans, init_out, init_name)
    # }

    plans
}
