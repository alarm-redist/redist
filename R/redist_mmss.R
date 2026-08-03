#' Multi-Merge-Split MCMC Redistricting Sampler
#'
#' `redist_mmss` uses a Markov Chain Monte Carlo algorithm that generalizes
#' the merge-split sampler (Carter et al. 2019) to merge and re-split `l`
#' districts at each step. When `l=2`, this reduces to the standard merge-split
#' algorithm. When `l=ndists`, this functions similar to a much less efficient SMC.
#' Higher values of `l` allow larger rearrangements of the plan at each step,
#' improving mixing at the cost of a lower acceptance rate.
#'
#' This function draws samples from a specific target measure, controlled by the
#' `map`, `compactness`, and `constraints` parameters.
#'
#' Key to ensuring good performance is monitoring the acceptance rate, which
#' is reported at the sample level in the output.
#' Users should also check diagnostics of the sample by running
#' `summary.redist_plans()`.
#'
#' Higher values of `compactness` sample more compact districts;
#' setting this parameter to 1 is computationally efficient and generates nicely
#' compact districts.
#'
#' @param map A [redist_map] object.
#' @param nsims The number of samples to draw, including warmup.
#' @param warmup The number of warmup samples to discard. Recommended to be at
#'   least the first 20% of samples, and in any case no less than around 100
#'   samples, unless initializing from a random plan.
#' @param thin Save every `thin`-th sample. Defaults to no thinning (1).
#' @param l The number of districts to merge and re-split at each step.
#'   Must be at least 2. Defaults to 3. Higher values allow larger moves but may
#'   have lower acceptance rates.
#' @param init_plan The initial state of the map, provided as a single vector
#'   to be shared across all chains, or a matrix with `chains` columns.
#'   If not provided, will default to the reference map of the `map` object, or
#'   if none exists, will sample a random initial state using [redist_smc].
#'   You can also request a random initial state for each chain by setting
#'   `init_plan="sample"`.
#' @param chains The number of chains to run. Each chain will have
#'   `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#'   with its own sampled plan. Defaults to 1.
#' @param counties A vector containing county (or other administrative or
#'   geographic unit) labels for each unit, which may be integers ranging from 1
#'   to the number of counties, or a factor or character vector.  If provided,
#'   the algorithm will generate maps that tend to follow county lines. Multiple
#'   nested levels may be supplied as a concatenated vector, ordered coarsest to
#'   finest, e.g. `counties = c(county, mcd)`.
#' @param project_plan_hier Whether to split administrative units into connected
#'   pseudo-units within the initial plan before sampling. This changes the
#'   administrative units used by the target. Defaults to `FALSE`.
#' @param enforce_hierarchy Whether to enforce the global hierarchical nesting
#'   condition when `counties` is supplied. If `FALSE`, administrative
#'   units guide the proposals and local connectivity only. Defaults to `FALSE`.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative. See
#'   the 'Details' section for more information, and computational
#'   considerations.
#' @param constraints A [redist_constr] object or list of constraints.
#' @param max_retries The maximum number of spanning trees to draw per proposal
#'   sequence before giving up. Higher values reduce proposal failures at the cost
#'   of more computation per step. The default of 200 is sufficient for most
#'   problems; increase for very tight population tolerances or large districts.
#' @param k_seq An integer vector giving the fixed top-eqn{k} window used at
#'   each sequential cut. A single value is reused at every cut; `NULL` uses
#'   the default sequence `l, l - 1, ..., 2`.
#' @param ncores The number of parallel processes to run. Defaults to the
#'   number of available cores, capped at the number of chains. Only used when
#'   `chains > 1`.
#' @param cl_type The cluster type (see [parallel::makeCluster()]). Safest is
#'   `"PSOCK"`, but `"FORK"` may be appropriate in some settings.
#' @param return_all If `TRUE` return all sampled plans; otherwise, just return
#'   the final plan from each chain.
#' @param init_name A name for the initial plan, or `FALSE` to not include
#'   the initial plan in the output.
#' @param verbose Whether to print out intermediate information while sampling.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @returns A [redist_plans] object containing the simulated plans.
#' @export
#'
#' @concept simulate
#' @references
#' Kenny, C. (2026). A Multi-District Markov Chain Monte Carlo Sampler for
#' Redistricting. osf.io/preprints/socarxiv/un8sk_v1.
#'
#' @examples
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#' sampled <- redist_mmss(fl_map, 100, l = 3)
redist_mmss <- function(
  map,
  nsims,
  warmup = 0,
  thin = 1L,
  l = 3L,
  init_plan = NULL,
  chains = 1L,
  counties = NULL,
  project_plan_hier = FALSE,
  enforce_hierarchy = FALSE,
  compactness = 1,
  constraints = list(),
  max_retries = 200L,
  k_seq = NULL,
  ncores = NULL,
  cl_type = "PSOCK",
  return_all = TRUE,
  init_name = NULL,
  verbose = FALSE,
  silent = FALSE
) {
  map <- validate_redist_map(map)
  V <- nrow(map)
  adj <- get_adj(map)
  ndists <- attr(map, "ndists")
  warmup <- max(warmup, 0L)
  thin <- as.integer(thin)
  chains <- as.integer(chains)
  l <- as.integer(l)

  if (
    !is.logical(enforce_hierarchy) ||
      length(enforce_hierarchy) != 1L ||
      is.na(enforce_hierarchy)
  ) {
    cli::cli_abort('{.arg enforce_hierarchy} must be a single TRUE or FALSE.')
  }
  if (
    !is.logical(project_plan_hier) ||
      length(project_plan_hier) != 1L ||
      is.na(project_plan_hier)
  ) {
    cli::cli_abort('{.arg project_plan_hier} must be a single TRUE or FALSE.')
  }

  hierarchy_mode <- 'speedup'
  if (enforce_hierarchy) {
    hierarchy_mode <- 'strict'
  }

  # Input validation
  if (compactness < 0) {
    cli::cli_abort("{.arg compactness} must be non-negative.")
  }
  if (nsims <= warmup) {
    cli::cli_abort("{.arg nsims} must be greater than {.arg warmup}.")
  }
  if (thin < 1 || thin > nsims - warmup) {
    cli::cli_abort(
      "{.arg thin} must be a positive integer, and no larger than {.arg nsims - warmup}."
    )
  }
  if (nsims < 1) {
    cli::cli_abort("{.arg nsims} must be positive.")
  }
  if (chains < 1) {
    cli::cli_abort("{.arg chains} must be positive.")
  }
  if (l < 2) {
    cli::cli_abort("{.arg l} must be at least 2.")
  }
  if (l > ndists) {
    cli::cli_abort(
      "{.arg l} must be at most the number of districts ({ndists})."
    )
  }

  # Set up initial plans for chains
  exist_name <- attr(map, "existing_col")
  counties <- rlang::eval_tidy(rlang::enquo(counties), map)
  raw_levels <- NULL
  if (!is.null(counties) && length(counties) == V) {
    raw_levels <- list(counties)
  } else if (
    !is.null(counties) && length(counties) > V && length(counties) %% V == 0L
  ) {
    raw_levels <- split(counties, rep(seq_len(length(counties) / V), each = V))
  } else if (!is.null(counties)) {
    cli::cli_abort(
      "{.arg counties} must have one value per unit in {.arg map}, or be a concatenated set of nested levels."
    )
  }
  hierarchy_levels <- normalize_hierarchy_levels(adj, raw_levels)
  if (enforce_hierarchy && is.null(hierarchy_levels)) {
    cli::cli_abort('{.arg enforce_hierarchy} = TRUE requires {.arg counties}.')
  }
  if (project_plan_hier && is.null(hierarchy_levels)) {
    cli::cli_abort('{.arg project_plan_hier} = TRUE requires {.arg counties}.')
  }
  init_counties <- NULL
  if (!is.null(hierarchy_levels)) {
    init_counties <- hierarchy_levels[, 1]
  }

  # Handle different init_plan scenarios
  if (is.null(init_plan)) {
    if (!is.null(exist_name)) {
      init_plans <- matrix(
        rep(vctrs::vec_group_id(get_existing(map)), chains),
        ncol = chains
      )
      if (is.null(init_name)) {
        init_names <- rep(exist_name, chains)
      } else {
        init_names <- rep(init_name, chains)
      }
    } else {
      init_plan <- "sample"
    }
  }

  if (
    is.null(init_plan) || (is.character(init_plan) && init_plan == "sample")
  ) {
    if (verbose) {
      cli::cli_inform("Sampling initial plans with SMC\n")
    }
    init_plans <- get_plans_matrix(
      redist_smc(
        map,
        chains,
        init_counties,
        compactness,
        constraints,
        resample = TRUE,
        ref_name = FALSE,
        verbose = verbose,
        silent = TRUE,
        ncores = 1
      )
    )
    if (is.null(init_name)) {
      init_names <- paste0("<sample ", seq_len(chains), ">")
    } else {
      init_names <- paste(init_name, seq_len(chains))
    }
  } else if (!is.null(init_plan) && !is.character(init_plan)) {
    if (is.matrix(init_plan)) {
      if (ncol(init_plan) < chains) {
        cli::cli_abort(
          "{.arg init_plan} matrix must have at least {chains} column{?s}."
        )
      }
      init_plans <- init_plan[, seq_len(chains), drop = FALSE]
    } else {
      init_plans <- matrix(rep(as.integer(init_plan), chains), ncol = chains)
    }

    if (is.null(init_name)) {
      init_names <- paste0("<init ", seq_len(chains), ">")
    } else if (is.matrix(init_plan) && chains > 1) {
      init_names <- paste(init_name, seq_len(chains))
    } else {
      init_names <- rep(init_name, chains)
    }
  }

  # Validate init_plans
  if (nrow(init_plans) != V) {
    cli::cli_abort(
      "{.arg init_plan} must be as long as the number of units as `map`."
    )
  }
  if (max(init_plans) != ndists) {
    cli::cli_abort(
      "{.arg init_plan} must have the same number of districts as `map`."
    )
  }
  if (any(apply(init_plans, 2, function(x) contiguity(adj, x)) != 1)) {
    cli::cli_warn("{.arg init_plan} should have contiguous districts.")
  }

  if (project_plan_hier) {
    hierarchy_levels <- project_plan_hierarchy(
      adj,
      init_plans,
      hierarchy_levels
    )
  }

  # Handle hierarchy labels
  if (is.null(hierarchy_levels)) {
    counties <- rep(1L, V)
    hierarchy_levels <- NULL
  } else {
    counties <- hierarchy_levels[, 1]
  }

  # Handle constraints
  if (!inherits(constraints, "redist_constr")) {
    constraints <- new_redist_constr(rlang::eval_tidy(
      rlang::enquo(constraints),
      map
    ))
  }
  if (any(c("edges_removed", "log_st") %in% names(constraints))) {
    cli::cli_warn(c(
      "{.var edges_removed} or {.var log_st} constraint found in
           {.arg constraints} and will be ignored.",
      ">" = "Adjust using {.arg compactness} instead."
    ))
  }
  if (
    any(c("polsby", "fry_hold") %in% names(constraints)) && compactness == 1
  ) {
    cli::cli_warn(
      "{.var polsby} or {.var fry_hold} constraint found in {.arg constraints}
                 with {.arg compactness == 1}. This may disrupt efficient sampling."
    )
  }
  constraints <- as.list(constraints)

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

  init_pop <- pop_tally(init_plans, pop, ndists)
  if (any(init_pop < pop_bounds[1]) || any(init_pop > pop_bounds[3])) {
    cli::cli_abort("Provided initialization does not meet population bounds.")
  }
  if (any(pop >= pop_bounds[3])) {
    too_big <- as.character(which(pop >= pop_bounds[3]))
    cli::cli_abort(c(
      "Unit{?s} {too_big} ha{?ve/s/ve}
                    population larger than the maximum district size.",
      "x" = "Redistricting impossible."
    ))
  }

  control_hierarchy_mode <- hierarchy_mode
  if (is.null(hierarchy_levels)) {
    control_hierarchy_mode <- "none"
  }
  control <- list(
    max_retries = as.integer(max_retries),
    hierarchy_mode = match(
      control_hierarchy_mode,
      c("none", "speedup", "strict")
    ) -
      1L
  )
  if (!is.null(hierarchy_levels)) {
    control$hierarchy_levels <- hierarchy_levels
  }
  if (!is.null(k_seq)) {
    control$k_seq <- as.integer(k_seq)
  }

  # Set up parallel cluster if needed
  if (is.null(ncores)) {
    ncores <- parallel::detectCores()
  }
  ncores <- min(ncores, chains)

  if (ncores > 1 && chains > 1) {
    `%oper%` <- `%dorng%`
    if (!silent) {
      of <- ifelse(
        Sys.info()[["sysname"]] == "Windows",
        tempfile(
          pattern = paste0("mmss_", substr(Sys.time(), 1, 10)),
          fileext = ".txt"
        ),
        ""
      )
      cl <- parallel::makeCluster(
        ncores,
        type = cl_type,
        outfile = of,
        methods = FALSE,
        useXDR = .Platform$endian != "little"
      )
    } else {
      cl <- parallel::makeCluster(
        ncores,
        type = cl_type,
        methods = FALSE,
        useXDR = .Platform$endian != "little"
      )
    }
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    `%oper%` <- `%do%`
  }

  # Run chains (in parallel or sequentially)
  out_par <- foreach::foreach(
    chain = seq_len(chains),
    .inorder = FALSE,
    .packages = "redist"
  ) %oper%
    {
      if (!silent) {
        cat("Starting chain ", chain, "\n", sep = "")
      }
      run_verbosity <- if (chain == 1 || verbosity == 3) verbosity else 0

      t1_run <- Sys.time()
      algout <- mmss_plans(
        nsims,
        adj,
        init_plans[, chain],
        counties,
        pop,
        ndists,
        pop_bounds[2],
        pop_bounds[1],
        pop_bounds[3],
        compactness,
        constraints,
        control,
        thin,
        l,
        run_verbosity
      )
      t2_run <- Sys.time()

      # Process output for this chain
      storage.mode(algout$plans) <- "integer"
      acceptances <- as.logical(algout$mhdecisions)

      l_diag <- list(
        runtime = as.numeric(t2_run - t1_run, units = "secs"),
        n_m_hit = algout$n_m_hit,
        n_attempts = algout$n_attempts,
        n_ust_draws = algout$n_ust_draws,
        n_ust_fail = algout$n_ust_fail,
        n_cut_fail = algout$n_cut_fail,
        n_hierarchical_merge_reject = algout$n_hierarchical_merge_reject,
        n_hierarchical_plan_reject = algout$n_hierarchical_plan_reject,
        n_zero_boundary = algout$n_zero_boundary,
        n_proposal_success = algout$n_proposal_success,
        k_seq = algout$k_seq,
        valid_cuts_by_step = algout$valid_cuts_by_step
      )

      warmup_idx <- c(seq_len(1 + warmup %/% thin), ncol(algout$plans))
      if (return_all) {
        algout$plans <- algout$plans[, -warmup_idx, drop = FALSE]
      } else {
        algout$plans <- algout$plans[, ncol(algout$plans) - 1L, drop = FALSE]
      }

      warmup_idx_acc <- c(seq_len(warmup %/% thin), length(acceptances))
      if (!return_all) {
        algout$mhdecisions <- as.logical(acceptances[length(acceptances) - 1L])
      } else {
        algout$mhdecisions <- acceptances[-warmup_idx_acc]
      }

      algout$l_diag <- l_diag
      algout$mh <- mean(as.logical(algout$mhdecisions))
      algout
    }

  # Combine results from all chains
  plans <- lapply(out_par, function(algout) algout$plans)
  each_len <- ncol(plans[[1]])
  plans <- do.call(cbind, plans)
  storage.mode(plans) <- "integer"

  mh <- sapply(out_par, function(algout) algout$mh)
  l_diag <- lapply(out_par, function(algout) algout$l_diag)
  acceptances <- sapply(out_par, function(algout) algout$mhdecisions)

  out <- new_redist_plans(
    plans = plans,
    map = map,
    algorithm = "mmss",
    wgt = NULL,
    resampled = FALSE,
    compactness = compactness,
    constraints = constraints,
    ndists = ndists,
    mh_acceptance = mh,
    version = packageVersion("redist"),
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

  # Add reference plans
  if (!is.null(init_names) && !isFALSE(init_name)) {
    if (chains == 1) {
      out <- add_reference(out, init_plans[, 1], init_names[1])
    } else if (all(init_names[1] == init_names)) {
      out <- add_reference(out, init_plans[, 1], init_names[1])
    } else {
      out <- Reduce(
        function(cur, idx) {
          add_reference(cur, init_plans[, idx], init_names[idx]) |>
            dplyr::mutate(chain = dplyr::coalesce(chain, idx))
        },
        rev(seq_len(chains)),
        init = out
      )
    }
  }

  if (chains > 1) {
    out <- dplyr::relocate(out, chain, .after = "draw")
  }

  out
}
