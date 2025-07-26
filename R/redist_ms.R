#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: parallel merge-split
####################################################

#' Parallel Merge-Split/Recombination MCMC Redistricting Sampler (Carter et al. 2019)
#'
#' \code{redist_mergesplit} uses a Markov Chain Monte Carlo algorithm (Carter et
#' al. 2019; based on DeFord et. al 2019) to generate congressional or legislative redistricting plans
#' according to contiguity, population, compactness, and administrative boundary
#' constraints. The MCMC proposal is the same as is used in the SMC sampler
#' (McCartan and Imai 2023); it is similar but not identical to those used in
#' the references.  1-level hierarchical Merge-split is supported through the
#' \code{counties} parameter and it has the same guarantees of a maximum number
#' of county splits as the SMCS algorithm.
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
#' @inheritParams redist_smc
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to generate. The chain will run for
#' `warmup+(nsims*thin)` steps.
#' @param chains the number of parallel chains to run. Each chain will have
#' `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#' with its own sampled plan.
#' @param warmup The number of warmup samples to discard. Recommended to be at
#' least the first 20% of samples, and in any case no less than around 100
#' samples, unless initializing from a random plan.
#' @param thin Save every `thin`-th sample after running warump.
#' Defaults to no thinning (1).
#' @param init_plan The initial state of the map. If not provided, will default to
#' the reference map of the \code{map} object, or if none exists, will sample
#' a random initial state using \code{\link{redist_smc}}. You can also request
#' a random initial state by setting \code{init_plan="sample"}.
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
#' @param init_name a name for the initial plan, or \code{FALSE} to not include
#' the initial plan in the output.  Defaults to the column name of the
#' existing plan, or "\code{<init>}" if the initial plan is sampled.
#' @param init_seats The initial number of seats for each district in `init_plan`.
#' @param verbose Whether to print out intermediate information while sampling.
#' Recommended.
#' @param silent Whether to suppress all diagnostic information.
#' @param ncores The number of clusters to spawn Defaults to the
#' maximum available detected by `parallel::detectCores()`.
#' @param cl_type the cluster type (see [makeCluster()]). Safest is `"PSOCK"`,
#' but `"FORK"` may be appropriate in some settings.
#' @param return_all if `TRUE` return all sampled plans; otherwise, just return
#' the final plan from each chain.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
#' the algorithm does not appear to be sampling from the target distribution.
#' Must be between 0 and 1.
#'
#' @returns A [`redist_plans`] object with all of the simulated plans, and an
#' additional `chain` column indicating the chain the plan was drawn from.
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
#' sampled_constr <- redist_mergesplit(fl_map, 10000, chains = 10,
#'     ncores = 2,
#'     constraints = list(
#'     incumbency = list(strength = 100, incumbents = c(3, 6, 25))
#'     )
#' )
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_mergesplit <- function(
  map,
  nsims,
  warmup = if (is.null(init_plan)) 10 else max(100, nsims %/% 5),
  thin = 1L,
  chains = 1L,
  init_plan = NULL,
  counties = NULL,
  compactness = 1,
  constraints = list(),
  constraint_fn = function(m) rep(0, ncol(m)),
  sampling_space = c("graph_plan", "spanning_forest", "linking_edge"),
  split_method = NULL,
  split_params = NULL,
  merge_prob_type = "uniform",
  init_seats = NULL,
  ncores = NULL,
  cl_type = "PSOCK",
  return_all = TRUE,
  init_name = NULL,
  verbose = FALSE,
  silent = FALSE,
  diagnostic_mode = FALSE,
  control = list(),
  adapt_k_thresh = .99
) {
  if (!missing(constraint_fn)) {
    cli::cli_warn("{.arg constraint_fn} is deprecated.")
  }

  if (!missing(adapt_k_thresh)) {
    cli::cli_warn(
      "Passing {.arg adapt_k_thresh} directly is deprecated. Pass it in as an argument
                 in {.arg split_params}"
    )
    if (is.list(split_params)) {
      split_params$adapt_k_thresh <- adapt_k_thresh
      split_params$estimate_cut_k <- TRUE
    } else {
      split_params <- list(
        adapt_k_thresh = adapt_k_thresh,
        estimate_cut_k = TRUE
      )
    }
  }

    # get parallel related params
    if (is.null(ncores)) {
        ncores <- parallel::detectCores()
        if (ncores <= 0) ncores <- 1
    }
    default_smc_ncores <- ncores
    ncores <- min(ncores, chains)


  # validate constraints
  constraints <- validate_constraints(
    map = map,
    constraints = rlang::enquo(constraints)
  )

  # get map params
  map_params <- get_map_parameters(map, !!rlang::enquo(counties))
  map <- map_params$map
  V <- map_params$V
  adj_list <- map_params$adj_list
  counties <- map_params$counties
  num_admin_units <- length(unique(counties))
  pop <- map_params$pop
  pop_bounds <- map_params$pop_bounds
  # get the total number of districts
  ndists <- map_params$ndists
  nseats <- map_params$nseats
  seats_range <- map_params$seats_range
  districting_scheme <- map_params$districting_scheme

    # check that no district is the sum of two others
    for (d_size1 in seats_range) {
      for (d_size2 in seats_range) {
          if((d_size1 + d_size2) %in% seats_range){
              cli::cli_abort("SMC does not support {.arg seats_range} where one district's seats is equal to the sum of two others")
          }
      }
    }


  thin <- as.integer(thin)
  chains <- as.integer(chains)

  if (compactness < 0) {
    cli::cli_abort("{.arg compactness} must be non-negative.")
  }
  if (thin < 1) {
    cli::cli_abort("{.arg thin} must be a positive integer.")
  }
  if (nsims < 1) {
    cli::cli_abort("{.arg nsims} must be positive.")
  }

  #validate the splitting method and params


  #validate the splitting method and params
  split_stuff_list <- validate_sample_space_and_splitting_method(
      sampling_space,
      split_method,
      split_params,
      1
  )

  sampling_space <- split_stuff_list$sampling_space
  split_method <- split_stuff_list$split_method
  forward_kernel_params <- split_stuff_list$forward_kernel_params



  exist_name <- attr(map, "existing_col")
  exist_seats <- attr(map, "existing_col_seats")
  if (is.null(init_plan)) {
    if (!is.null(exist_name)) {
      init_plan <- matrix(
        rep(vctrs::vec_group_id(get_existing(map)), chains),
        ncol = chains,
        nrow = length(get_existing(map))
      )
      if (is.null(init_name)) {
        init_names <- rep(exist_name, chains)
      } else {
        init_names <- rep(init_name, chains)
      }
    } else {
      init_plan <- "sample"
    }
  } else if (!is.null(init_plan)) {
    if (inherits(init_plan, "redist_plans")) {
      if (is.null(init_seats)) {
        init_seats <- get_seats_matrix(init_plan)
      }
      init_plan <- get_plans_matrix(init_plan)
    } else if (is.matrix(init_plan)) {
      stopifnot(ncol(init_plan) == chains)
      init_plan <- init_plan
    } else {
      init_plan <- matrix(rep(as.integer(init_plan), chains), ncol = chains)
    }
    if (is.null(init_name)) {
      init_names <- paste0("<init> ", seq_len(chains))
    } else {
      init_names <- rep(init_name, chains)
    }
  }

  if (isTRUE(init_plan == "sample")) {
    if (!silent) {
      cat("Sampling initial plans with SMC\n")
    }
    # heuristic. Do at least 50 plans to not get stuck
    n_smc_nsims <- max(chains, 50)
    # get ncores if needed
    if (is.list(control) && "init_ncores" %in% names(control)) {
      init_ncores <- control[["init_ncores"]]
    } else {
        # else default to default smc ncores from earlier
      init_ncores <- default_smc_ncores
    }

    init_plan <- redist_smc(
      map,
      n_smc_nsims,
      counties,
      compactness,
      constraints,
      resample = TRUE,
      split_params = split_params,
      sampling_space = sampling_space,
      split_method = split_method,
      ncores = init_ncores,
      ref_name = FALSE,
      verbose = verbose,
      silent = silent
    )

    sampled_inidices <- sample.int(n = n_smc_nsims, size = chains, replace = F)

    init_seats <- get_seats_matrix(init_plan)[,
      sampled_inidices,
      drop = FALSE
    ]

    init_plan <- get_plans_matrix(
      init_plan
    )[, sampled_inidices, drop = FALSE]

    if (is.null(init_name)) {
      init_names <- paste0("<init> ", seq_len(chains))
    } else {
      init_names <- paste(init_name, seq_len(chains))
    }
  } else {
    if (is.null(init_seats)) {
      if (districting_scheme == "SMD") {
        init_seats <- matrix(1L, nrow = ndists, ncol = ncol(init_plan))
      } else {
          init_seats <- replicate(chains, exist_seats)
      }
    }
  }


  # check it satsifies population bounds
  # add one because we assume 1 indexed
  init_pop <- pop_tally(init_plan, pop, ndists)


  # subtract 1 to make it 0 indexed
  init_plan <- init_plan - 1
  # validate initial plans
  validate_initial_region_id_mat(init_plan, V, chains, ndists)


  bad_pops <- sapply(seq_len(chains),
         function(i){
             any(init_pop[,i] < pop_bounds[1] * init_seats[,i]) ||
                 any(init_pop[,i] > pop_bounds[3] * init_seats[,i])
         }
         ) |>
    any()

  if (bad_pops) {
    cli::cli_abort("Provided initialization does not meet population bounds.")
  }

  verbosity <- 1
  if (verbose) {
    verbosity <- 3
  }
  if (silent) {
    verbosity <- 0
  }

  control = list(
    splitting_method = split_method,
    do_mh = TRUE
  )

  # add the splitting parameters
  # need to do it like this because its a list
  control <- c(control, forward_kernel_params)

  multiprocess <- ncores > 1

  if (multiprocess) {
    `%oper%` <- `%dorng%`

    of <- ifelse(
      Sys.info()[['sysname']] == 'Windows',
      tempfile(
        pattern = paste0('ms_', substr(Sys.time(), 1, 10)),
        fileext = '.txt'
      ),
      ''
    )
    # this makes a cluster using socket (NOT FORK) with
    if (!silent) {
      cl <- makeCluster(
        ncores,
        outfile = of,
        methods = FALSE,
        useXDR = .Platform$endian != "little"
      )
    } else {
      cl <- makeCluster(
        ncores,
        methods = FALSE,
        useXDR = .Platform$endian != "little"
      )
    }

    doParallel::registerDoParallel(cl, cores = ncores)
    on.exit(stopCluster(cl))


    # Ensures only one process will print when there's multiple processes
    # to avoid cluttering output
    parallel::clusterEvalQ(cl, {
      if (!exists("is_chain1", envir = .GlobalEnv)) {
        is_chain1 <- FALSE
      }
      NULL
    })
  } else {
    `%oper%` <- `%do%`
  }

  t1 <- Sys.time()
  out_par <- foreach(
    chain = seq_len(chains),
    .inorder = FALSE,
    .packages = "redist"
  ) %oper%
    {
      if (chain == 1) {
        is_chain1 <- T
      }
      run_verbosity <- if (is_chain1 || !multiprocess) verbosity else 0
      if (!silent && is_chain1) {
        cat("Starting chain ", chain, "\n", sep = "")
        # flush.console()
      }

      t1_run <- Sys.time()
      algout <- ms_plans(
        nsims = nsims,
        warmup = warmup,
        thin = thin,
        ndists = ndists,
        total_seats = nseats,
        district_seat_sizes = seats_range,
        adj_list = adj_list,
        counties = counties,
        pop = pop,
        target = pop_bounds[2],
        lower = pop_bounds[1],
        upper = pop_bounds[3],
        rho = compactness,
        init_plan = init_plan[, chain, drop = FALSE],
        init_seats = init_seats[, chain, drop = FALSE],
        sampling_space_str = sampling_space,
        merge_prob_type = merge_prob_type,
        control = control,
        constraints = constraints,
        verbosity = run_verbosity,
        diagnostic_mode = diagnostic_mode
      )
      t2_run <- Sys.time()

      # Internal diagnostics,
      algout$internal_diagnostics <- list(
        log_mh_ratio = algout$log_mh_ratio,
        mh_acceptance = algout$mhdecisions,
        warmup_acceptances = algout$warmup_acceptances,
        post_warump_acceptances = algout$post_warump_acceptances,
        warmup_accept_rate = algout$warmup_acceptances / warmup,
        postwarmup_accept_rate = algout$post_warump_acceptances /
          (nsims * thin),
        tree_sizes = algout$tree_sizes,
        successful_tree_sizes = algout$successful_tree_sizes,
        proposed_plans = algout$proposed_plans
      )

      # add cut k for graph space
      run_forward_kernel_params <- forward_kernel_params
      if(sampling_space == GRAPH_PLAN_SPACE_SAMPLING){
          run_forward_kernel_params$cut_k_used <- algout$cut_k
      }

      algout$l_diag <- list(
        runtime = as.numeric(t2_run - t1_run, units = "secs"),
        nsims = nsims,
        thin = thin,
        warmup = warmup,
        total_acceptances = (algout$warmup_acceptances +
          algout$post_warump_acceptances),
        accept_rate = (algout$warmup_acceptances +
          algout$post_warump_acceptances) /
          algout$total_steps,
        total_steps = algout$total_steps,
        forward_kernel_params = run_forward_kernel_params
      )

      # Information about the run
      algout$run_information <- list(
        valid_region_sizes_to_split_list = algout$valid_region_sizes_to_split_list,
        valid_split_region_sizes_list = algout$valid_split_region_sizes_list,
        sampling_space = sampling_space,
        split_method = split_method,
        merge_prob_type = merge_prob_type,
        nsims = nsims,
        alg_name = "mergesplit"
      )

      # flatten the region sizes by column
      dim(algout$seats) <- NULL
      dim(algout$region_pops) <- NULL

      algout
    }
  t2 <- Sys.time()

  plans <- lapply(out_par, function(algout) {
    algout$plans
  })
  each_len <- ncol(plans[[1]])
  plans <- do.call(cbind, plans)

  seats <- do.call(c, lapply(out_par, function(x) x$seats))
  region_pops <- do.call(c, lapply(out_par, function(x) x$region_pops))


  mh <- sapply(out_par, function(algout) {
    mean(as.logical(algout$mhdecisions))
  })
  l_diag <- lapply(out_par, function(algout) algout$l_diag)
  run_information <- lapply(out_par, function(x) x$run_information)
  internal_diagnostics <- lapply(out_par, function(x) x$internal_diagnostics)

  acceptances <- sapply(out_par, function(algout) {
    algout$mhdecisions
  })

  out <- new_redist_plans(
    plans = plans,
    map = map,
    algorithm = "mergesplit",
    wgt = NULL,
    inputs_safe = TRUE,
    resampled = FALSE,
    distr_pop = region_pops,
    seats = seats,
    compactness = compactness,
    constraints = constraints,
    ndists = ndists,
    mh_acceptance = mh,
    version = packageVersion("redist"),
    diagnostics = l_diag,
    run_information = run_information,
    internal_diagnostics = internal_diagnostics,
    pop_bounds = pop_bounds,
    total_runtime = t2 - t1
  ) %>%
    mutate(
      chain = rep(seq_len(chains), each = each_len * ndists),
      mcmc_accept = rep(acceptances, each = ndists)
    )

  if (!is.null(init_names) && !isFALSE(init_name)) {
    if (all(init_names[1] == init_names)) {
      out <- add_reference(
          plans = out,
          ref_plan = init_plan[, 1],
          name = init_names[1],
          ref_seats = init_seats[, 1]
          )
    } else {
      out <- Reduce(
        function(cur, idx) {
          add_reference(plans = cur, ref_plan = init_plan[, idx],
                        name = init_names[idx], ref_seats = init_seats[, idx]) %>%
            mutate(chain = dplyr::coalesce(chain, idx))
        },
        rev(seq_len(chains)),
        init = out
      )
    }
  }

  dplyr::relocate(out, chain, .after = "draw")
}

utils::globalVariables("chain")
