#' CycleWalk MCMC Redistricting Sampler
#'
#' `redist_cyclewalk` uses a Markov Chain Monte Carlo algorithm based on
#' the CycleWalk method to generate redistricting plans. The algorithm maintains
#' a spanning forest representation of districts and proposes changes by finding
#' cycles between adjacent districts.
#'
#' This function draws samples from a specific target measure, controlled by the
#' `map`, `compactness`, and `constraints` parameters.
#'
#' @param map A [redist_map] object.
#' @param nsims The number of samples to draw, not including warmup.
#' @param warmup The number of warmup samples to discard.
#' @param thin Save every `thin`-th sample. Defaults to no thinning (1).
#' @param init_plan The initial state of the map. If not provided, will default to
#'   the reference map of the `map` object, or if none exists, will sample
#'   a random initial state using [redist_smc]. You can also request
#'   a random initial state by setting `init_plan="sample"`.
#' @param counties A vector containing county (or other administrative or
#'   geographic unit) labels for each unit, which may be integers ranging from 1
#'   to the number of counties, or a factor or character vector. If provided,
#'   the algorithm will generate maps that tend to follow county lines.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative.
#' @param constraints A [redist_constr] object or list of constraints.
#' @param edge_weights Optional list of edge weights for the graph. Each element
#'   should be a list with two fields: \code{edge} (a length-2 numeric vector of
#'   vertex indices) and \code{weight} (a positive number). Edges not specified
#'   default to weight 1.0. Higher weights make edges less likely to be cut or
#'   linked in proposals (cost interpretation). Example:
#'   \code{list(list(edge = c(1, 2), weight = 2.0))}.
#' @param init_name A name for the initial plan, or `FALSE` to not include
#'   the initial plan in the output. Defaults to the column name of the
#'   existing plan, or `<init>` if the initial plan is sampled.
#' @param verbose Whether to print out intermediate information while sampling.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @returns A [redist_plans] object containing the simulated plans.
#'
#' @examples
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#' sampled <- redist_cyclewalk(fl_map, 100)
#'
#' @concept simulate
#' @md
#' @export
redist_cyclewalk <- function(map, nsims,
                             warmup = 0,
                             thin = 1L, init_plan = NULL,
                             counties = NULL, compactness = 1,
                             constraints = list(),
                             edge_weights = NULL,
                             init_name = NULL,
                             verbose = FALSE, silent = FALSE) {

    map <- validate_redist_map(map)
    V <- nrow(map)
    adj <- get_adj(map)
    ndists <- attr(map, "ndists")
    warmup <- max(warmup, 0L)
    thin <- as.integer(thin)

    # Input validation
    if (compactness < 0) {
        cli::cli_abort("{.arg compactness} must be non-negative.")
    }
    if (nsims <= warmup) {
        cli::cli_abort("{.arg nsims} must be greater than {.arg warmup}.")
    }
    if (thin < 1 || thin > nsims - warmup) {
        cli::cli_abort("{.arg thin} must be a positive integer, and no larger than {.arg nsims - warmup}.")
    }
    if (nsims < 1) {
        cli::cli_abort("{.arg nsims} must be positive.")
    }

    # Handle initial plan
    exist_name <- attr(map, "existing_col")
    counties <- rlang::eval_tidy(rlang::enquo(counties), map)
    orig_lookup <- seq_len(ndists)

    if (is.null(init_plan) && !is.null(exist_name)) {
        init_plan <- vctrs::vec_group_id(get_existing(map))
        orig_lookup <- unique(get_existing(map))
        if (is.null(init_name)) init_name <- exist_name
    } else if (!is.null(init_plan) && is.null(init_name)) {
        init_name <- "<init>"
    }

    if (length(init_plan) == 0L || isTRUE(init_plan == "sample")) {
        init_plan <- as.integer(get_plans_matrix(
            redist_smc(map, 10, counties, resample = FALSE, ref_name = FALSE, silent = TRUE, ncores = 1))[, 1])
        if (is.null(init_name)) init_name <- "<init>"
    }

    # Validate init_plan
    if (length(init_plan) != V) {
        cli::cli_abort("{.arg init_plan} must be as long as the number of units as `map`.")
    }
    if (max(init_plan) != ndists) {
        cli::cli_abort("{.arg init_plan} must have the same number of districts as `map`.")
    }
    if (any(contiguity(adj, init_plan) != 1)) {
        cli::cli_warn("{.arg init_plan} should have contiguous districts.")
    }

    # Handle counties
    if (is.null(counties)) {
        counties <- rep(1L, V)
    } else {
        if (any(is.na(counties)))
            cli::cli_abort("County vector must not contain missing values.")

        # Handle discontinuous counties
        component <- contiguity(adj, vctrs::vec_group_id(counties))
        counties <- dplyr::if_else(component > 1,
                                   paste0(as.character(counties), "-", component),
                                   as.character(counties)) |>
            vctrs::vec_group_id()
    }

    # Handle constraints
    if (!inherits(constraints, "redist_constr")) {
        constraints <- new_redist_constr(rlang::eval_tidy(rlang::enquo(constraints), map))
    }
    constraints <- as.list(constraints)

    # Handle edge weights
    if (!is.null(edge_weights)) {
        edge_weights <- validate_edge_weights(edge_weights, adj, V)
    } else {
        edge_weights <- list()  # Empty list for C++
    }

    # Set verbosity
    verbosity <- 1
    if (verbose) {
        verbosity <- 3
    }
    if (silent) {
        verbosity <- 0
    }

    # Population bounds
    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]

    # Validate population
    init_pop <- pop_tally(matrix(init_plan, ncol = 1), pop, ndists)
    if (any(init_pop < pop_bounds[1]) | any(init_pop > pop_bounds[3])) {
        cli::cli_abort("Provided initialization does not meet population bounds.")
    }
    if (any(pop >= pop_bounds[3])) {
        too_big <- as.character(which(pop >= pop_bounds[3]))
        cli::cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                    population larger than the maximum district size.",
            "x" = "Redistricting impossible."))
    }

    # Call C++ implementation
    t1_run <- Sys.time()

    control <- list()  # For future control parameters

    algout <- cyclewalk_plans(
        N = nsims,
        l = adj,
        init = init_plan,
        counties = counties,
        pop = pop,
        n_distr = ndists,
        target = pop_bounds[2],
        lower = pop_bounds[1],
        upper = pop_bounds[3],
        compactness = compactness,
        constraints = constraints,
        control = control,
        edge_weights = edge_weights,
        thin = thin,
        verbosity = verbosity
    )

    t2_run <- Sys.time()

    # Process output
    storage.mode(algout$plans) <- "integer"

    # Handle MH decisions if present
    acceptances <- if (!is.null(algout$mhdecisions)) {
        as.logical(algout$mhdecisions)
    } else {
        rep(TRUE, ncol(algout$plans))
    }

    # Calculate warmup indices to remove
    warmup_idx <- c(seq_len(1 + warmup %/% thin), ncol(algout$plans))

    l_diag <- list(
        runtime = as.numeric(t2_run - t1_run, units = "secs")
    )

    # Create redist_plans object
    out <- new_redist_plans(
        algout$plans[, -warmup_idx, drop = FALSE],
        map, "cyclewalk", NULL, FALSE,
        ndists = ndists,
        compactness = compactness,
        constraints = constraints,
        version = packageVersion("redist"),
        diagnostics = l_diag,
        mh_acceptance = mean(acceptances, na.rm = TRUE)
    )

    # Add acceptance info
    warmup_idx_acc <- c(seq_len(warmup %/% thin), length(acceptances))
    out <- out |> dplyr::mutate(mcmc_accept = rep(acceptances[-warmup_idx_acc], each = ndists))

    # Add reference plan if requested
    if (!is.null(init_name) && !isFALSE(init_name)) {
        out <- add_reference(out, init_plan, init_name)
    }

    out
}

# Validate edge weights format and values
validate_edge_weights <- function(edge_weights, adj, V) {
    if (!is.list(edge_weights)) {
        cli::cli_abort("{.arg edge_weights} must be a list.")
    }

    if (length(edge_weights) == 0) {
        return(list())  # Empty list is valid (no custom weights)
    }

    # Each element should be a list with $edge and $weight
    for (i in seq_along(edge_weights)) {
        entry <- edge_weights[[i]]

        if (!is.list(entry)) {
            cli::cli_abort("Entry {i} of {.arg edge_weights} must be a list.")
        }

        # Check for edge field
        if (!"edge" %in% names(entry)) {
            cli::cli_abort("Entry {i} of {.arg edge_weights} missing {.field edge} field.")
        }

        # Check for weight field
        if (!"weight" %in% names(entry)) {
            cli::cli_abort("Entry {i} of {.arg edge_weights} missing {.field weight} field.")
        }

        edge <- entry$edge
        weight <- entry$weight

        # Validate edge format
        if (!is.numeric(edge) || length(edge) != 2) {
            cli::cli_abort("Entry {i}: {.field edge} must be a numeric vector of length 2.")
        }

        u <- as.integer(edge[1])
        v <- as.integer(edge[2])

        # Validate vertices in range
        if (u < 1 || u > V || v < 1 || v > V) {
            cli::cli_abort("Entry {i}: vertices {u} and {v} out of range [1, {V}].")
        }

        # Check edge exists in adjacency
        # Note: adj is 0-indexed (neighbors stored as 0, 1, 2, ...)
        # but user input is 1-indexed, so check (v-1) in adj[[u]]
        if (!((v - 1) %in% adj[[u]])) {
            cli::cli_abort("Entry {i}: edge ({u}, {v}) not in adjacency graph.")
        }

        # Validate weight
        if (!is.numeric(weight) || length(weight) != 1) {
            cli::cli_abort("Entry {i}: {.field weight} must be a single number.")
        }

        if (weight <= 0) {
            cli::cli_abort("Entry {i}: {.field weight} must be positive, got {weight}.")
        }
    }

    edge_weights
}

