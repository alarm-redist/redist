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
#' @param nsims The number of samples to draw per chain, not including warmup.
#' @param chains The number of parallel chains to run. Each chain will have
#'   `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#'   with its own sampled plan. Defaults to 1.
#' @param warmup The number of warmup samples to discard.
#' @param thin Save every `thin`-th sample. Defaults to no thinning (1).
#' @param instep Number of MCMC iterations per recorded sample (default 10).
#'   Higher values improve mixing but increase runtime. Matches Julia's `instep`.
#' @param cycle_walk_frac Fraction of proposals that are cycle walks vs internal
#'   forest walks (default 0.1). Cycle walks change the partition; forest walks
#'   only change the spanning tree representation.
#' @param init_plan The initial state of the map, provided as a single vector
#'   to be shared across all chains, or a matrix with `chains` columns.
#'   If not provided, will default to the reference map of the `map` object, or if
#'   none exists, will sample a random initial state using [redist_smc]. You can
#'   also request a random initial state for each chain by setting
#'   `init_plan="sample"`.
#' @param counties A vector containing county (or other administrative or
#'   geographic unit) labels for each unit, which may be integers ranging from 1
#'   to the number of counties, or a factor or character vector. If provided and
#'   `edge_weights` is NULL, automatically creates edge weights that upweight
#'   intra-county edges by 10x, encouraging plans that follow county lines.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative. A
#'   value of 0 samples uniformly over partitions. Default is `0`.
#' @param constraints A [redist_constr] object or list of constraints.
#' @param edge_weights Edge weights for the graph. Can be:
#'   - A single number: used as the intra-county weight multiplier (requires
#'     `counties`). For example, `edge_weights = 5` with `counties` upweights
#'     intra-county edges by 5x.
#'   - A list of edge weight specifications, each with `edge` (length-2 vertex
#'     indices) and `weight` (positive number). Unspecified edges default to 1.0.
#'   - NULL (default): if `counties` is provided, auto-generates 10x weights for
#'     intra-county edges; otherwise no weighting.
#' @param ncores The number of parallel processes to run. Defaults to the
#'   number of available cores, capped at the number of chains.
#' @param cl_type The cluster type (see [parallel::makeCluster()]). Safest is `"PSOCK"`,
#'   but `"FORK"` may be appropriate in some settings.
#' @param return_all If `TRUE` return all sampled plans; otherwise, just return
#'   the final plan from each chain.
#' @param init_name A name for the initial plan, or `FALSE` to not include
#'   the initial plan in the output. Defaults to the column name of the
#'   existing plan, or `<init>` if the initial plan is sampled.
#' @param verbose Whether to print out intermediate information while sampling.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @returns A [redist_plans] object containing the simulated plans. If `chains > 1`,
#'   the output will include a `chain` column indicating which chain each plan
#'   came from.
#'
#' @examples
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#' sampled <- redist_cyclewalk(fl_map, 100)
#'
#' # Multiple chains for convergence diagnostics
#' sampled_chains <- redist_cyclewalk(fl_map, 200, chains = 2, ncores = 2)
#'
#' @concept simulate
#' @md
#' @export
redist_cyclewalk <- function(map, nsims,
                             chains = 1,
                             warmup = 0,
                             thin = 1L,
                             instep = 10L,
                             cycle_walk_frac = 0.1,
                             init_plan = NULL,
                             counties = NULL, compactness = 1,
                             constraints = list(),
                             edge_weights = NULL,
                             ncores = NULL,
                             cl_type = 'PSOCK',
                             return_all = TRUE,
                             init_name = NULL,
                             verbose = FALSE, silent = FALSE) {
  map <- validate_redist_map(map)
  V <- nrow(map)
  adj <- get_adj(map)
  ndists <- attr(map, 'ndists')
  warmup <- max(warmup, 0L)
  thin <- as.integer(thin)
  instep <- as.integer(instep)
  chains <- as.integer(chains)

  # Input validation
  if (compactness < 0) {
    cli::cli_abort('{.arg compactness} must be non-negative.')
  }
  if (nsims <= warmup) {
    cli::cli_abort('{.arg nsims} must be greater than {.arg warmup}.')
  }
  if (thin < 1 || thin > nsims - warmup) {
    cli::cli_abort('{.arg thin} must be a positive integer, and no larger than {.arg nsims - warmup}.')
  }
  if (instep < 1) {
    cli::cli_abort('{.arg instep} must be a positive integer.')
  }
  if (cycle_walk_frac < 0 || cycle_walk_frac > 1) {
    cli::cli_abort('{.arg cycle_walk_frac} must be between 0 and 1.')
  }
  if (nsims < 1) {
    cli::cli_abort('{.arg nsims} must be positive.')
  }
  if (chains < 1) {
    cli::cli_abort('{.arg chains} must be positive.')
  }

  # Set up initial plans for chains
  exist_name <- attr(map, 'existing_col')
  counties <- rlang::eval_tidy(rlang::enquo(counties), map)

  # Handle different init_plan scenarios
  if (is.null(init_plan)) {
    if (!is.null(exist_name)) {
      init_plans <- matrix(rep(vctrs::vec_group_id(get_existing(map)), chains), ncol = chains)
      if (is.null(init_name)) {
        init_names <- rep(exist_name, chains)
      } else {
        init_names <- rep(init_name, chains)
      }
    } else {
      init_plan <- 'sample'
    }
  }

  if (is.null(init_plan) || (is.character(init_plan) && init_plan == 'sample')) {
    if (verbose) {
      cli::cli_inform('Sampling initial plans with SMC\n')
    }
    init_plans <- get_plans_matrix(
      redist_smc(map, chains, counties, compactness, constraints,
        resample = TRUE, ref_name = FALSE, verbose = verbose,
        silent = TRUE, ncores = 1
      )
    )
    if (is.null(init_name)) {
      init_names <- paste0('<sample ', seq_len(chains), '>')
    } else {
      init_names <- paste(init_name, seq_len(chains))
    }
  } else if (!is.null(init_plan)) {
    if (is.matrix(init_plan)) {
      if (ncol(init_plan) != chains) {
        cli::cli_abort('{.arg init_plan} matrix must have {chains} column{?s}.')
      }
      init_plans <- init_plan
    } else {
      init_plans <- matrix(rep(as.integer(init_plan), chains), ncol = chains)
    }

    if (is.null(init_name)) {
      init_names <- paste0('<init ', seq_len(chains), '>')
    } else if (is.matrix(init_plan) && chains > 1) {
      # Matrix with unique inits per chain - add suffixes
      init_names <- paste(init_name, seq_len(chains))
    } else {
      # Single vector replicated across chains
      init_names <- rep(init_name, chains)
    }
  }

  # Validate init_plans
  if (nrow(init_plans) != V) {
    cli::cli_abort('{.arg init_plan} must be as long as the number of units as `map`.')
  }
  if (max(init_plans) != ndists) {
    cli::cli_abort('{.arg init_plan} must have the same number of districts as `map`.')
  }
  if (any(apply(init_plans, 2, function(x) contiguity(adj, x)) != 1)) {
    cli::cli_warn('{.arg init_plan} should have contiguous districts.')
  }

  # Handle counties
  if (is.null(counties)) {
    counties <- rep(1L, V)
  } else {
    if (any(is.na(counties))) {
      cli::cli_abort('County vector must not contain missing values.')
    }

    # Handle discontinuous counties
    component <- contiguity(adj, vctrs::vec_group_id(counties))
    counties <- dplyr::if_else(component > 1,
      paste0(as.character(counties), '-', component),
      as.character(counties)
    ) |>
      vctrs::vec_group_id()
  }

  # Handle constraints
  if (!inherits(constraints, 'redist_constr')) {
    if (length(constraints) == 0) {
      # Empty list or NULL - use map data for constraint evaluation
      constraints <- redist_constr(map)
    } else {
      constraints <- new_redist_constr(rlang::eval_tidy(rlang::enquo(constraints), map))
    }
  }

  # Warn if user manually added log_st or edges_removed
  if (any(c('edges_removed', 'log_st') %in% names(constraints))) {
    cli::cli_warn(c('{.var edges_removed} or {.var log_st} constraint found in
           {.arg constraints} and will be ignored.',
      '>' = 'Adjust using {.arg compactness} instead.'
    ))
  }

  constraints <- as.list(constraints)

  # Handle edge weights
  # If edge_weights is a single number and counties provided, use as intra-county weight
  # If edge_weights is a list, validate and use directly
  # If NULL and counties provided, default to 10x for intra-county edges
  if (is.numeric(edge_weights) && length(edge_weights) == 1) {
    if (all(counties == 1L)) {
      cli::cli_abort('{.arg edge_weights} as a number requires {.arg counties} to be specified.')
    }
    edge_weights <- build_county_edge_weights(adj, counties, weight = edge_weights)
    if (!silent) {
      cli::cli_inform('Using county-based edge weights ({edge_weights[[1]]$weight}x for intra-county edges).')
    }
  } else if (!is.null(edge_weights)) {
    edge_weights <- validate_edge_weights(edge_weights, adj, V)
  } else if (!all(counties == 1L)) {
    # Build edge weights from counties: upweight intra-county edges by 10x
    edge_weights <- build_county_edge_weights(adj, counties, weight = 10.0)
    if (!silent) {
      cli::cli_inform('Using county-based edge weights (10x for intra-county edges).')
    }
  } else {
    edge_weights <- list()
  }

  # Set verbosity
  verbosity <- 1
  if (verbose) verbosity <- 3
  if (silent) verbosity <- 0

  # Population bounds
  pop_bounds <- attr(map, 'pop_bounds')
  pop <- map[[attr(map, 'pop_col')]]

  # Validate population for all init plans
  init_pop <- pop_tally(init_plans, pop, ndists)
  if (any(init_pop < pop_bounds[1]) | any(init_pop > pop_bounds[3])) {
    cli::cli_abort('Provided initialization does not meet population bounds.')
  }
  if (any(pop >= pop_bounds[3])) {
    too_big <- as.character(which(pop >= pop_bounds[3]))
    cli::cli_abort(c('Unit{?s} {too_big} ha{?ve/s/ve}
                    population larger than the maximum district size.',
      'x' = 'Redistricting impossible.'
    ))
  }

  # Set up parallel cluster if needed
  if (is.null(ncores)) ncores <- parallel::detectCores()
  ncores <- min(ncores, chains)

  if (ncores > 1 && chains > 1) {
    `%oper%` <- `%dorng%`
    if (!silent) {
      of <- ifelse(Sys.info()[['sysname']] == 'Windows',
        tempfile(pattern = paste0('cw_', substr(Sys.time(), 1, 10)), fileext = '.txt'),
        ''
      )
      cl <- parallel::makeCluster(ncores,
        type = cl_type, outfile = of, methods = FALSE,
        useXDR = .Platform$endian != 'little'
      )
    } else {
      cl <- parallel::makeCluster(ncores,
        type = cl_type, methods = FALSE,
        useXDR = .Platform$endian != 'little'
      )
    }

    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    `%oper%` <- `%do%`
  }

  # Run chains (in parallel or sequentially)
  out_par <- foreach::foreach(
    chain = seq_len(chains), .inorder = FALSE,
    .packages = 'redist'
  ) %oper% {
    if (!silent) cat('Starting chain ', chain, '\n', sep = '')
    run_verbosity <- if (chain == 1 || verbosity == 3) verbosity else 0

    t1_run <- Sys.time()
    control <- list()

    algout <- cyclewalk_plans(
      N = nsims,
      l = adj,
      init = init_plans[, chain],
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
      instep = instep,
      cycle_walk_frac = cycle_walk_frac,
      verbosity = run_verbosity
    )

    t2_run <- Sys.time()

    # Process output for this chain
    storage.mode(algout$plans) <- 'integer'

    acceptances <- if (!is.null(algout$mhdecisions)) {
      as.logical(algout$mhdecisions)
    } else {
      rep(NA, ncol(algout$plans))
    }

    # Extract diagnostics
    l_diag <- list(
      runtime = as.numeric(t2_run - t1_run, units = 'secs')
    )

    if (!is.null(algout$diagnostics)) {
      l_diag$accept_prob <- algout$diagnostics$accept_prob
      l_diag$cycle_length <- algout$diagnostics$cycle_length
      l_diag$n_valid_cuts <- algout$diagnostics$n_valid_cuts
      l_diag$failure_modes <- algout$diagnostics$failure_modes
    }

    # Calculate warmup indices
    warmup_idx <- c(seq_len(1 + warmup %/% thin), ncol(algout$plans))
    if (return_all) {
      algout$plans <- algout$plans[, -warmup_idx, drop = FALSE]
      warmup_idx_acc <- c(seq_len(warmup %/% thin), length(acceptances))
      algout$mhdecisions <- acceptances[-warmup_idx_acc]
    } else {
      algout$plans <- algout$plans[, ncol(algout$plans) - 1L, drop = FALSE]
      algout$mhdecisions <- as.logical(acceptances[length(acceptances) - 1L])
    }

    algout$l_diag <- l_diag
    algout$mh <- mean(as.logical(algout$mhdecisions), na.rm = TRUE)
    algout
  }

  # Combine results from all chains
  plans <- lapply(out_par, function(algout) algout$plans)
  each_len <- ncol(plans[[1]])
  plans <- do.call(cbind, plans)
  storage.mode(plans) <- 'integer'

  mh <- sapply(out_par, function(algout) algout$mh)
  l_diag <- lapply(out_par, function(algout) algout$l_diag)
  acceptances <- sapply(out_par, function(algout) algout$mhdecisions)

  # Create output
  out <- new_redist_plans(
    plans = plans,
    map = map,
    algorithm = 'cyclewalk',
    wgt = NULL,
    resampled = FALSE,
    compactness = compactness,
    constraints = constraints,
    ndists = ndists,
    mh_acceptance = mh,
    version = packageVersion('redist'),
    diagnostics = l_diag
  )

  # Add chain column and acceptance info
  if (chains > 1) {
    out <- out |>
      dplyr::mutate(
        chain = rep(seq_len(chains), each = each_len * ndists),
        mcmc_accept = rep(acceptances, each = ndists)
      )
  } else {
    out <- out |>
      dplyr::mutate(mcmc_accept = rep(acceptances, each = ndists))
  }

  # Add reference plans if requested
  if (!is.null(init_names) && !isFALSE(init_name)) {
    if (chains == 1) {
      out <- add_reference(out, init_plans[, 1], init_names[1])
    } else if (all(init_names[1] == init_names)) {
      # Shared init across chains - add once without chain assignment
      out <- add_reference(out, init_plans[, 1], init_names[1])
    } else {
      # Unique init per chain
      out <- Reduce(function(cur, idx) {
        add_reference(cur, init_plans[, idx], init_names[idx]) |>
          dplyr::mutate(chain = dplyr::coalesce(chain, idx))
      }, rev(seq_len(chains)), init = out)
    }
  }

  if (chains > 1) {
    out <- dplyr::relocate(out, chain, .after = 'draw')
  }

  out
}

# Validate edge weights format and values
validate_edge_weights <- function(edge_weights, adj, V) {
  if (!is.list(edge_weights)) {
    cli::cli_abort('{.arg edge_weights} must be a list.')
  }

  if (length(edge_weights) == 0) {
    return(list()) # Empty list is valid (no custom weights)
  }

  # Each element should be a list with $edge and $weight
  for (i in seq_along(edge_weights)) {
    entry <- edge_weights[[i]]

    if (!is.list(entry)) {
      cli::cli_abort('Entry {i} of {.arg edge_weights} must be a list.')
    }

    # Check for edge field
    if (!'edge' %in% names(entry)) {
      cli::cli_abort('Entry {i} of {.arg edge_weights} missing {.field edge} field.')
    }

    # Check for weight field
    if (!'weight' %in% names(entry)) {
      cli::cli_abort('Entry {i} of {.arg edge_weights} missing {.field weight} field.')
    }

    edge <- entry$edge
    weight <- entry$weight

    # Validate edge format
    if (!is.numeric(edge) || length(edge) != 2) {
      cli::cli_abort('Entry {i}: {.field edge} must be a numeric vector of length 2.')
    }

    u <- as.integer(edge[1])
    v <- as.integer(edge[2])

    # Validate vertices in range
    if (u < 1 || u > V || v < 1 || v > V) {
      cli::cli_abort('Entry {i}: vertices {u} and {v} out of range [1, {V}].')
    }

    # Check edge exists in adjacency
    # Note: adj is 0-indexed (neighbors stored as 0, 1, 2, ...)
    # but user input is 1-indexed, so check (v-1) in adj[[u]]
    if (!((v - 1) %in% adj[[u]])) {
      cli::cli_abort('Entry {i}: edge ({u}, {v}) not in adjacency graph.')
    }

    # Validate weight
    if (!is.numeric(weight) || length(weight) != 1) {
      cli::cli_abort('Entry {i}: {.field weight} must be a single number.')
    }

    if (weight <= 0) {
      cli::cli_abort('Entry {i}: {.field weight} must be positive, got {weight}.')
    }
  }

  edge_weights
}

# Build edge weights from county assignments
# Edges within the same county get the specified weight; others default to 1.0
build_county_edge_weights <- function(adj, counties, weight = 10.0) {
  edge_weights <- list()
  V <- length(adj)

  for (u in seq_len(V)) {
    # adj is 0-indexed, so neighbors are stored as 0, 1, 2, ...
    for (v_idx in adj[[u]]) {
      v <- v_idx + 1L # convert to 1-indexed
      # Only process each edge once (u < v)
      if (u < v && counties[u] == counties[v]) {
        edge_weights <- c(edge_weights, list(list(
          edge = c(u, v),
          weight = weight
        )))
      }
    }
  }

  edge_weights
}
