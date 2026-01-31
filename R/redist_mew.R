#' Marked Edge Walk MCMC Redistricting Sampler (McWhorter and DeFord 2024)
#'
#' `redist_mew` uses a Markov Chain Monte Carlo algorithm based on marked
#' edge walks on spanning trees (McWhorter and DeFord 2024) to generate
#' congressional or legislative redistricting plans according to contiguity,
#' population, compactness, and other custom constraints.
#'
#' @param map A [redist_map] object.
#' @param nsims The number of samples to draw, including warmup.
#' @param warmup The number of warmup samples to discard.
#' @param thin Save every `thin`-th sample. Defaults to no thinning (1).
#' @param init_plan The initial state of the map. If not provided, will default
#'   to the reference map of the `map` object, or if none exists, will
#'   sample a random initial state using [redist_smc()]. You can also
#'   request a random initial state by setting `init_plan="sample"`.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative.
#' @param constraints A [redist_constr] object or a list containing
#'   information on sampling constraints. See [redist_constr] for more
#'   information.
#' @param init_name A name for the initial plan, or `FALSE` to not include
#'   the initial plan in the output. Defaults to the column name of the existing
#'   plan, or `"<init>"` if the initial plan is sampled.
#' @param chains The number of parallel chains to run. Each chain will have
#'   `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#'   with its own sampled plan. Defaults to 1 (no parallelization).
#' @param ncores The number of parallel processes to run. Defaults to the
#'   number of available cores. Only used if `chains > 1`.
#' @param cl_type The cluster type (see [parallel::makeCluster()]). Safest is `"PSOCK"`,
#'   but `"FORK"` may be appropriate in some settings. Only used if `chains > 1`.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @returns `redist_mew` returns a [redist_plans] object containing the
#'   simulated plans.
#'
#' @references
#' McWhorter, A., & DeFord, D. (2025). The Marked Edge Walk: A Novel MCMC
#' Algorithm for Sampling of Graph Partitions.
#'
#' @examples
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#'
#' redist_mew(fl_map, 200)
#'
#' redist_mew(fl_map, 200, chains = 2, ncores = 2)
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
                       init_name = NULL,
                       chains = 1L,
                       ncores = NULL,
                       cl_type = 'PSOCK',
                       verbose = FALSE,
                       silent = FALSE) {
  map <- validate_redist_map(map)
  V <- nrow(map)
  adj <- get_adj(map)
  ndists <- attr(map, 'ndists')
  warmup <- max(warmup, 0L)
  thin <- as.integer(thin)
  chains <- as.integer(chains)

  # Validate parameters
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

  # Process initial plan
  exist_name <- attr(map, 'existing_col')
  orig_lookup <- seq_len(ndists)

  # Handle init_plan for potentially multiple chains
  if (is.null(init_plan)) {
    if (!is.null(exist_name)) {
      init_plan <- vctrs::vec_group_id(get_existing(map))
      orig_lookup <- unique(get_existing(map))
      if (is.null(init_name)) init_name <- exist_name
    } else {
      init_plan <- 'sample'
    }
  } else if (!is.null(init_plan) && is.null(init_name)) {
    init_name <- '<init>'
  }

  # Convert to matrix format for chains
  if (is.matrix(init_plan)) {
    if (ncol(init_plan) != chains) {
      cli::cli_abort('{.arg init_plan} matrix must have {chains} column{?s}.')
    }
    init_plans <- init_plan
    if (is.null(init_name)) {
      init_names <- paste0('<init> ', seq_len(chains))
    } else {
      init_names <- rep(init_name, chains)
    }
  } else if (length(init_plan) == 0L || isTRUE(init_plan == 'sample')) {
    if (!silent && chains > 1) cli::cli_alert_info('Sampling initial plans with SMC...')
    if (!silent && chains == 1) cli::cli_alert_info('Sampling initial plan with SMC...')
    init_plans <- tryCatch(
      {
        as.matrix(get_plans_matrix(
          redist_smc(map, chains,
            resample = chains > 1, ref_name = FALSE,
            silent = TRUE, ncores = 1
          )
        ))
      },
      error = function(e) {
        cli::cli_abort(c(
          'Failed to generate initial plan{?s} using SMC.',
          'i' = 'This can happen on very small or difficult maps.',
          'i' = 'Try providing an {.arg init_plan} manually.',
          'x' = conditionMessage(e)
        ))
      }
    )
    if (is.null(init_name)) {
      init_names <- if (chains == 1) '<init>' else paste0('<init> ', seq_len(chains))
    } else {
      init_names <- if (chains == 1) init_name else paste(init_name, seq_len(chains))
    }
  } else {
    # Single vector - replicate across chains
    init_plans <- matrix(rep(as.integer(init_plan), chains), ncol = chains)
    if (is.null(init_name)) {
      init_names <- if (chains == 1) '<init>' else paste0('<init> ', seq_len(chains))
    } else {
      init_names <- rep(init_name, chains)
    }
  }

  # Validate initial plans
  if (nrow(init_plans) != V) {
    cli::cli_abort('{.arg init_plan} must be as long as the number of units as `map`.')
  }
  if (max(init_plans) != ndists) {
    cli::cli_abort('{.arg init_plan} must have the same number of districts as `map`.')
  }
  if (any(apply(init_plans, 2, function(x) contiguity(adj, x)) != 1)) {
    cli::cli_abort('{.arg init_plan} must be contiguous.')
  }

  # Convert init_plans to 1-indexed integers
  for (i in seq_len(chains)) {
    init_plans[, i] <- match(init_plans[, i], sort(unique(init_plans[, i])))
  }

  # Process constraints and population parameters
  if (inherits(constraints, 'redist_constr')) {
    constraints <- as.list(constraints)
  }

  pop_bounds <- attr(map, 'pop_bounds')
  pop <- map[[attr(map, 'pop_col')]]

  # Setup control parameters
  control <- list(
    verbosity = ifelse(silent, 0L, ifelse(verbose, 2L, 1L))
  )

  # Print sampling information
  if (!silent) {
    cli::cli_rule(left = 'Sampling using Marked Edge Walk')
    cli::cli_alert_info('{ndists} districts on {V} units')
    if (chains > 1) {
      cli::cli_alert_info('Running {chains} chains')
    }
  }

  # Run sampling algorithm
  if (chains == 1) {
    # Single chain - no parallelization
    time_start <- proc.time()
    if (!silent && !is.null(init_names) && init_names[1] != FALSE) {
      cli::cli_alert_info('Initial plan: {init_names[1]}')
    }

    plans_raw <- mew_plans(
      nsims = nsims,
      adj = adj,
      init = init_plans[, 1],
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
    time_elapsed <- (proc.time() - time_start)[['elapsed']]

    out_list <- list(list(plans_raw = plans_raw, time_elapsed = time_elapsed))
  } else {
    # Multiple chains - use parallelization
    if (is.null(ncores)) ncores <- parallel::detectCores()
    ncores <- min(ncores, chains)

    of <- ifelse(Sys.info()[['sysname']] == 'Windows',
      tempfile(pattern = paste0('mew_', substr(Sys.time(), 1, 10)), fileext = '.txt'),
      ''
    )
    if (!silent) {
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

    out_list <- foreach::foreach(chain = seq_len(chains), .inorder = FALSE, .packages = 'redist') %dorng% {
      if (!silent) cat('Starting chain ', chain, '\n', sep = '')
      run_verbosity <- if (chain == 1 || control$verbosity == 2L) control$verbosity else 0L
      run_control <- list(verbosity = run_verbosity)

      time_start <- proc.time()
      plans_raw <- mew_plans(
        nsims = nsims,
        adj = adj,
        init = init_plans[, chain],
        pop = pop,
        n_distr = ndists,
        target = pop_bounds[2],
        lower = pop_bounds[1],
        upper = pop_bounds[3],
        rho = compactness,
        constraints = constraints,
        control = run_control,
        thin = thin,
        verbosity = run_verbosity
      )
      time_elapsed <- (proc.time() - time_start)[['elapsed']]

      list(plans_raw = plans_raw, time_elapsed = time_elapsed)
    }
  }

  # Post-process results
  # Process each chain's output
  plans_list <- lapply(out_list, function(out) {
    plans_raw <- out$plans_raw
    time_elapsed <- out$time_elapsed

    # Convert plans to original district labels
    plans_m <- plans_raw$plans
    for (i in seq_len(ncol(plans_m))) {
      plans_m[, i] <- orig_lookup[plans_m[, i]]
    }

    # Remove warmup samples
    if (warmup > 0) {
      n_warmup_thinned <- ceiling(warmup / thin)
      if (n_warmup_thinned < ncol(plans_m)) {
        plans_m_final <- plans_m[, -(1:n_warmup_thinned), drop = FALSE]
      } else {
        plans_m_final <- plans_m[, ncol(plans_m), drop = FALSE]
      }
    } else {
      plans_m_final <- plans_m
    }

    # Process diagnostics
    diag_every <- 100
    n_diag <- length(plans_raw$accept_rate)
    warmup_diag_idx <- seq_len(pmax(1, warmup %/% diag_every))

    if (length(warmup_diag_idx) < n_diag) {
      diag_idx <- setdiff(seq_len(n_diag), warmup_diag_idx)
    } else {
      diag_idx <- n_diag
    }

    l_diag <- list(
      accept_rate = plans_raw$accept_rate[diag_idx],
      cycle_intersect_rate = plans_raw$cycle_intersect_rate[diag_idx],
      avg_proposal_tries = plans_raw$avg_proposal_tries[diag_idx],
      runtime = time_elapsed
    )

    list(plans = plans_m_final, diagnostics = l_diag)
  })

  # Combine plans from all chains
  plans_m_all <- do.call(cbind, lapply(plans_list, function(x) x$plans))
  each_len <- ncol(plans_list[[1]]$plans)

  # Handle diagnostics - if single chain, flatten to match expected format
  if (chains == 1) {
    l_diag_final <- plans_list[[1]]$diagnostics
  } else {
    l_diag_final <- lapply(plans_list, function(x) x$diagnostics)
  }

  # Construct output
  plans <- new_redist_plans(
    plans = plans_m_all,
    map = map,
    algorithm = 'mew',
    wgt = rep(1, ncol(plans_m_all)),
    ndists = ndists,
    compactness = compactness,
    constraints = constraints,
    version = packageVersion('redist'),
    diagnostics = l_diag_final
  )

  # Add chain column if multiple chains
  if (chains > 1) {
    plans <- plans %>%
      dplyr::mutate(chain = rep(seq_len(chains), each = each_len * ndists)) %>%
      dplyr::relocate(chain, .after = 'draw')
  }

  # Add reference plan if requested
  if (!is.null(init_names) && !isFALSE(init_name)) {
    if (chains == 1 || all(init_names[1] == init_names)) {
      plans <- add_reference(plans, init_plans[, 1], init_names[1])
    } else {
      plans <- Reduce(function(cur, idx) {
        add_reference(cur, init_plans[, idx], init_names[idx]) %>%
          dplyr::mutate(chain = dplyr::coalesce(chain, idx))
      }, rev(seq_len(chains)), init = plans)
    }
  }

  plans
}

utils::globalVariables('chain')
