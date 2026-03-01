#' Balanced Up-Down Walk MCMC Redistricting Sampler
#'
#' `redist_bud` uses a Markov Chain Monte Carlo algorithm based on
#' the Balanced Up-Down Walk method to generate redistricting plans. The algorithm
#' maintains a spanning forest representation of districts with a district-level
#' tree connected by marked edges, and proposes changes by forming multi-district
#' cycles and resampling balanced cuts.
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
#' @param thin Save every `thin`-th sample. Defaults to 10.
#' @param init_plan The initial state of the map, provided as a single vector
#'   to be shared across all chains, or a matrix with `chains` columns.
#'   If not provided, will default to the reference map of the `map` object, or if
#'   none exists, will sample a random initial state using [redist_smc]. You can
#'   also request a random initial state for each chain by setting
#'   `init_plan="sample"`.
#' @param counties A vector containing county (or other administrative or
#'   geographic unit) labels for each unit, which may be integers ranging from 1
#'   to the number of counties, or a factor or character vector.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative. A
#'   value of 0 samples uniformly over partitions. Default is `1`.
#' @param constraints A [redist_constr] object or list of constraints.
#' @param ncores The number of parallel processes to run. Defaults to the
#'   number of available cores, capped at the number of chains.
#' @param cl_type The cluster type (see [parallel::makeCluster()]). Safest is `"PSOCK"`,
#'   but `"FORK"` may be appropriate in some settings.
#' @param return_all If `TRUE` return all sampled plans; otherwise, just return
#'   the final plan from each chain.
#' @param init_name A name for the initial plan, or `FALSE` to not include
#'   the initial plan in the output.  Defaults to the column name of the
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
#' sampled <- redist_bud(fl_map, 100)
#'
#' @concept simulate
#' @md
#' @export
redist_bud <- function(map, nsims,
                       chains = 1,
                       warmup = 0,
                       thin = 10L,
                       init_plan = NULL,
                       counties = NULL, compactness = 1,
                       constraints = list(),
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
      init_names <- paste(init_name, seq_len(chains))
    } else {
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
      constraints <- redist_constr(map)
    } else {
      constraints <- new_redist_constr(rlang::eval_tidy(rlang::enquo(constraints), map))
    }
  }
  if (any(c('edges_removed', 'log_st') %in% names(constraints))) {
    cli::cli_warn(c('{.var edges_removed} or {.var log_st} constraint found in
           {.arg constraints} and will be ignored.',
      '>' = 'Adjust using {.arg compactness} instead.'
    ))
  }
  constraints <- as.list(constraints)

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
        tempfile(pattern = paste0('bud_', substr(Sys.time(), 1, 10)), fileext = '.txt'),
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

    algout <- bud_plans(
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
      edge_weights = list(),
      thin = thin,
      instep = 1L,
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

    l_diag <- list(
      runtime = as.numeric(t2_run - t1_run, units = 'secs')
    )

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
    algorithm = 'bud',
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
      out <- add_reference(out, init_plans[, 1], init_names[1])
    } else {
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
