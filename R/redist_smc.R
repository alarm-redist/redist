##############################################
## Author: Philip O'Sullivan
## Institution: Harvard University
## Date Created: 2025/02/15
## Purpose: Wrapper for running smc cpp code
##############################################



#' Generalized SMCS Redistricting Sampler (O'Sullivan, McCartan and Imai forthcoming)
#'
#'
#' `redist_smc` uses a Sequential Monte Carlo Sampler algorithm
#' (O'Sullivan, McCartan and Imai forthcoming) to generate representative samples of
#' congressional or legislative redistricting plans according to
#' contiguity, population, compactness, and administrative boundary constraints.
#'
#' This function draws samples from a specific target measure controlled by
#' the `map`, `compactness`, and `constraints` parameters.
#'
#' Key to ensuring good performance is monitoring the efficiency of the resampling
#' process at each SMC stage.  Unless `silent=FALSE`, this function will print
#' out the effective sample size of each resampling step to allow the user to
#' monitor the efficiency.  If `verbose=TRUE` the function will also print
#' out information on any relevant splitting parameters chosen and the
#' acceptance rate (based on the population constraint) at each step.
#' Users should also check diagnostics of the sample by running
#' \code{summary.redist_plans()}.
#'
#' Higher values of `compactness` sample more compact districts;
#' setting this parameter to 1 is computationally efficient and generates nicely
#' compact districts.  Values of other than 1 may lead to highly variable
#' importance sampling weights.
#'
#' @param map A [redist_map()] object.
#' @param nsims The number of samples to draw.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided,
#' the algorithm will only generate maps which split up to `ndists-1`
#' counties. Even if there are fewer counties than `ndists - 1`, the spanning
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
#' @param ncores How many threads to use to parallelize plan generation within each
#' process. The default, 0, will use the number of available cores on the machine
#' as long as `nsims` and the number of units is large enough. If `runs>1`
#' you will need to set this manually. If more than one core is used, the
#' sampler output will not be fully reproducible with `set.seed()`. If full
#' reproducibility is desired, set `ncores=1` and `nproc=1`.
#' @param init_particles Either a [redist_plans] object or a matrix of partial
#' plans to begin sampling from. For advanced use only. The matrix must have
#' `nsims` columns and a row for every precinct. It is important to ensure that
#' the existing districts meet contiguity and population constraints, or there
#' may be major issues when sampling.
#' @param init_seats A matrix of the number of seats per region of
#' the partial plans to begin sampling from. For advanced use only. The matrix
#' must have `nsims` columns, as many rows as there are regions in `init_particles`
#' and each column must sum to the total number of seats in the map. Not needed
#' if `init_particles` is a [redist_plans] object. If `init_particles` are passed
#' but not `init_seats` then the number of seats will be attempted
#' to be inferred.
#' @param init_weights A vector of length `nsims` of unnormalized plan weights
#' associated with `init_particles`. The weights must all be strictly positive.
#' If no weights are passed in then they will all be set to 1. If the
#' `init_plans` were not resampled then it is recommended to pass their weights in.
#' Not needed if `init_particles` is a [redist_plans] object.
#' @param sampling_space The space to sample the plans on. This does not affect
#' the plans output by the function but the sample space used can have a large
#' impact on computational cost/runtime and convergence. Current spaces supported
#' right now are
#'  - `r GRAPH_PLAN_SPACE_SAMPLING` : graph partition space
#'  - `r FOREST_SPACE_SAMPLING` : spanning forest space
#'  - `r LINKING_EDGE_SPACE_SAMPLING` : linking edge forest space
#' @param split_method The method used for splitting spanning trees in the
#' sampling process. When sampling on the space of graph partitions it must be
#' the naive top k method but any method is allowed for forest space sampling.
#' @param split_params A list of parameters related to splitting the plans.
#' Options include
#' \itemize{
#'  \item \code{splitting_schedule} What rule to use for selecting splitting
#'  sizes. The final target distribution is the same regardless of splitting
#'  schedule but the intermediate distributions can change. Current options
#'  include
#'  \itemize{
#'      \item \code{split_district_only} At each step split off a single district.
#'      \item \code{any_valid_sizes} At each step allow for regions to be split
#'      into any two sizes (assuming the sizes can eventually be split into
#'      districts.) Currently this is only supported for single member districting.
#'  }
#' }
#' Parameters for \code{split_method} Any relevant parameters for the
#'  \code{split_method}. These include the following
#'  \itemize{
#'      \item `r NAIVE_K_SPLITTING` parameters:
#'        \itemize{
#'      \item \code{adapt_k_thresh} The threshold value used in the heuristic to
#'       select a value `k_i` for each splitting iteration for graph space
#'       sampling if estimation is desired (the `k_i` can also be manually passed in.)
#'       Higher values are more accurate but may require more
#'       computation. Set to 1 for the most conservative sampling.
#'       Must be between 0 and 1.
#'      \item \code{manual_k_params} The `k_i` values to be used for each splitting
#'      iteration for graph space sampling. Beware when specifying manual values
#'      it is crucial they are close to the true values as too small `k_i` values
#'      will cause the algorithm to fail to sample from the target distribution
#'      correctly and too large values will cause a drastic performance hit. The
#'      input must either be a single integer to use for each step or a vector
#'      of integers equal to the number of smc steps.
#'      }
#'      \item `r EXP_BIGGER_ABS_DEV_SPLITTING`
#'      \itemize{
#'          \item \code{splitting_alpha} When selecting an edge to cut in the
#'          tree a valid edge is selected with probability proportional to
#'          \code{exp(-splitting_alpha * max_dev)}. \code{splitting_alpha} can
#'          be any real number. Values closer to zero result in more stable weights
#'          and larger values result in more unstable weights.
#'      }
#'  }
#' @param ms_params A list of mergesplit parameters.
#' \itemize{
#'  \item \code{ms_frequency}: How often to run merge steps. Should either be an integer
#' (meaning run after every _ smc steps) or a vector of 1 indexed step numbers
#' indicating which smc steps to run merge split. A value of -1 means just run
#' after all smc steps have been run. A value of 1 means run after every smc step.
#' \item \code{ms_moves_multiplier} Multiplier to the baseline number of mergesplit
#' moves to be performed each step. For a mergesplit step the baseline number of
#' moves is calculated as the ceiling of 1 over the previous mergesplit steps
#' acceptance rate (or smc step if no prior mergesplit steps were done). The
#' total number of moves is `ceiling(ms_moves_multiplier * baseline_num_moves)`.
#' \item \code{merge_prob_type} What probability to use to select regions to merge
#' in the mergesplit kernel. Defaults to giving all pairs equal probability.
#' }
#' @param n_steps How many steps to run the SMC algorithm for.
#' Each step splits off a new region. Defaults to all remaining districts.
#' If fewer than the number of remaining splits, reference plans are disabled.
#' @param seq_alpha The amount to adjust the weights by at each resampling step;
#' higher values prefer exploitation, while lower values prefer exploration.
#' Must be between 0 and 1. If any mergesplit steps are applied then it must be
#' set to 1.
#' @param pop_temper The strength of the automatic population tempering. Try
#' values of 0.01-0.05 to start if the algorithm gets stuck on the final few
#' splits.
#' @param ref_name a name for the existing plan, which will be added as a
#' reference plan, or `FALSE` to not include the initial plan in the
#' output. Defaults to the column name of the existing plan.
#' @param verbose Whether to print out intermediate information while sampling.
#' Recommended.
#' @param silent Whether to suppress all diagnostic information.
#' @param diagnostics What amount of diagnostic information to save. Current
#' options are
#' \itemize{
#'  \item \code{basic} Saves the final plans and for each step the
#'  incremental weights and ancestors.
#'  \item \code{all} Saves everything \code{basic} does in addition to the
#'  intermediate plan and seat matrices after each splitting and mergesplit step.
#'  Use with caution as this can use a lot of memory very quickly.
#' }
#' @param control A list of optional advanced parameters.
#' \itemize{
#'  \item \code{nproc}: The number of processes (independent instances of R)
#' spawned to simulate the plans. The processes execute runs in parallel, each
#' using `ncores` threads. If more than one process is used, the
#' sampler output will not be fully reproducible with `set.seed()`. If full
#' reproducibility is desired, set `nproc=1` and
#' `ncores = 1`. If missing defaults to a single process.
#'  \item \code{weight_type} The type of SMC weights to use. Optimal weights typically
#' have lower variance and lead to faster convergence but can be more
#' computationally expensive, especially for computationally complex constraints
#' or when `compactness` is not set to 1.
#' }
#' @param adapt_k_thresh Deprecated. Pass in through the `split_params` arg.
#' The threshold value used in the heuristic to select a
#' value `k_i` for each splitting iteration. Higher values are more accurate
#' but may require more computation. Set to 1 for the most conservative
#' sampling. Must be between 0 and 1.
#'
#' @return `redist_smc` returns a [redist_plans] object containing the simulated
#' plans.
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
#' # Two runs with multiple processes
#' redist_smc(fl_map, 1000, runs = 2, control = list(nproc = 2))
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_smc <- function(
  map,
  nsims,
  counties = NULL,
  compactness = 1,
  constraints = list(),
  resample = TRUE,
  runs = 1L,
  ncores = 0L,
  init_particles = NULL,
  init_seats = NULL,
  init_weights = NULL,
  sampling_space = c("graph_plan", "spanning_forest", "linking_edge"),
  split_method = NULL,
  split_params = NULL,
  ms_params = list(),
  n_steps = NULL,
  seq_alpha = 1L,
  pop_temper = 0,
  ref_name = NULL,
  verbose = FALSE,
  silent = FALSE,
  diagnostics = c("basic", "all"),
  control = list(weight_type = "optimal", nproc = 1L),
  adapt_k_thresh = .99
) {
  if (!missing(adapt_k_thresh)) {
    cli_warn(
      "Passing {.arg adapt_k_thresh} directly is deprecated. Pass it in as an argument
                 in {.arg split_params}"
    )
    if (is.list(split_params)) {
      split_params$adapt_k_thresh <- adapt_k_thresh
      split_params$estimate_cut_k <- TRUE
    } else {
      split_params = list(
        adapt_k_thresh = adapt_k_thresh,
        estimate_cut_k = TRUE
      )
    }
  }

  if (!is_scalar(compactness) || compactness < 0) {
    cli_abort("{.arg compactness} must be non-negative.")
  }
  if (seq_alpha <= 0 || seq_alpha > 1 || !is_scalar(seq_alpha)) {
    cli_abort("{.arg seq_alpha} must lie in (0, 1].")
  }
  if (nsims < 1) {
    cli_abort("{.arg nsims} must be positive.")
  }

  # check default inputs
  diagnostics <- rlang::arg_match(diagnostics)
  diagnostic_level <- dplyr::case_when(
    diagnostics == "basic" ~ 0,
    diagnostics == "all" ~ 1,
    .default = 0
  )

  # validate constraints
  constraints <- validate_constraints(map, rlang::enquo(constraints))
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

  # check that the seat sizes are a range
  if (
    any(
        seats_range !=
        seq.int(from = min(seats_range), to = max(seats_range))
    )
  ) {
    cli_abort(
      "For {.arg seats_range} only a continuous range of district seat sizes are allowed!"
    )
  }
  # check that no district is the sum of two other
  for (d_size1 in seats_range) {
      for (d_size2 in seats_range) {
        if((d_size1 + d_size2) %in% seats_range){
            cli_abort("SMC does not support {.arg seats_range} where one district's seats is equal to the sum of two others")
        }
      }
  }

  # get the splitting size regime
  splitting_size_regime <- get_splitting_schedule(split_params, districting_scheme)

  # get initial plan parameters
  initial_plan_params <- get_init_plan_params(
        nsims, nseats, pop, pop_bounds,
        init_particles, init_seats, init_weights
  )
  init_particles <- initial_plan_params$init_particles
  init_seats <- initial_plan_params$init_seats
  initial_log_weights <- initial_plan_params$initial_log_weights
  init_num_regions <- initial_plan_params$init_num_regions


  if (is.null(n_steps)) {
    n_steps <- ndists - init_num_regions
  }
  final_dists <- init_num_regions + n_steps
  if (final_dists > ndists) {
    cli_abort("Too many districts already drawn to take {n_steps} steps.")
  }

  #validate the splitting method and params
  split_stuff_list <- validate_sample_space_and_splitting_method(
      sampling_space,
      split_method,
      split_params,
      n_steps
  )

  sampling_space <- split_stuff_list$sampling_space
  split_method <- split_stuff_list$split_method
  forward_kernel_params <- split_stuff_list$forward_kernel_params


  # linking edge with counties not supported right now
  if (num_admin_units > 1 && sampling_space == LINKING_EDGE_SPACE_SAMPLING) {
      cli_abort(
          "Linking Edge Sampling with counties is not supported right now"
      )
  }



  total_smc_steps <- n_steps

  # get merge split parameter information
  ms_params_list <- extract_ms_params(ms_params, total_smc_steps)
  run_ms <- ms_params_list$run_ms
  merge_prob_type <- ms_params_list$merge_prob_type
  ms_moves_multiplier <- ms_params_list$ms_moves_multiplier
  ms_frequency <- ms_params_list$ms_frequency
  merge_split_step_vec <- ms_params_list$merge_split_step_vec


  # get the step types
  step_types <- ifelse(merge_split_step_vec, "ms", "smc")
  if(sum(!merge_split_step_vec) != total_smc_steps){
    cli_abort("In {.field step_types} the number of smc steps must be equal to {.field total_smc_steps}")
  }
  # assert first step is not smc
  if(merge_split_step_vec[1]){
      cli_abort("The first step cannot be mergesplit! An SMC step must be run first!")
  }
  total_ms_steps <- sum(merge_split_step_vec)
  # total number of steps to run
  total_steps <- total_smc_steps + total_ms_steps


  # if run_ms is true then seq_alpha must be 1
  if (run_ms && seq_alpha != 1) {
      seq_alpha <- 1L
      cli_warn(
          "{.arg seq_alpha} must be set to 1 if any mergesplit steps are being run!"
      )
  }

  # compute lags thing
  lags <- 1 + unique(round((ndists - 1)^0.8 * seq(0, 0.7, length.out = 4)^0.9))

  # verbosity stuff
  verbosity <- 1
  if (verbose) {
    verbosity <- 3
  }
  if (silent) {
    verbosity <- 0
  }

  # set up parallel processing stuff
  ncores_max <- parallel::detectCores()
  # if 0 cores just assume only single threaded machine
  if (ncores_max <= 0) {
    ncores_max <- 1
  }

  # ncores is the number of threads per process
  # to avoid confusion: ncores is to the number of threads assigned to each process.
  # So a call will use ncores * nproc total threads
  if (!is.null(ncores)) {
    if (
      !rlang::is_integerish(ncores) ||
        !is_scalar(ncores)
    ) {
      cli_abort("{.arg ncores} must be a single integer!")
    } else if (ncores == 0) {
      ncores <- ncores_max
    } else if (ncores < 0) {
      cli_abort(
        "{.arg ncores} can't be negative!"
      )
    }
  } else {
    ncores <- ncores_max
  }

  # Now extract control parameters
  control_params_list <- extract_control_params(control)
  nproc <- control_params_list[["nproc"]]
  weight_type <- control_params_list[["weight_type"]]

  multiprocess <- nproc > 1
  # make sure we're not spawning more proccesses than runs
  if (multiprocess) {
    nproc <- min(runs, nproc)
  }

  # warn if more processes than cores
  if (nproc > ncores_max) {
    cli_warn(
      "Inputted number of processes to spawn is greater than detected number of cores on machine"
    )
  }

  nproc <- as.integer(nproc)
  ncores <- as.integer(ncores)

  if (nproc > 1 && runs > 1) {
    `%oper%` <- `%dorng%`
    of <- if (Sys.info()[["sysname"]] == "Windows") {
      tempfile(
        pattern = paste0("smc_", substr(Sys.time(), 1, 10)),
        fileext = ".txt"
      )
    } else {
      ""
    }

    # this makes a cluster using socket (NOT FORK) with
    if (!silent) {
      cl <- makeCluster(
        nproc,
        outfile = of,
        methods = FALSE,
        useXDR = .Platform$endian != "little"
      )
    } else {
      cl <- makeCluster(
        nproc,
        methods = FALSE,
        useXDR = .Platform$endian != "little"
      )
    }
    # this makes it avoid printing the loading required package message each time
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(foreach))
      suppressPackageStartupMessages(library(rngtools))
      suppressPackageStartupMessages(library(redist))
    })
    # Makes it so only one process will print but if more runs then processes
    # it doesn't just print once
    parallel::clusterEvalQ(cl, {
      if (!exists("is_chain1", envir = .GlobalEnv)) {
        is_chain1 <- FALSE
      }
      NULL
    })
    doParallel::registerDoParallel(cl, cores = nproc)
    on.exit(stopCluster(cl))
    if (!silent) cat("Spawning ", nproc, " clusters \n")
  } else {
    `%oper%` <- `%do%`
  }

  cpp_control_list <- list(
    weight_type = weight_type,
    lags = lags,
    seq_alpha = seq_alpha,
    pop_temper = pop_temper,
    num_threads = as.integer(ncores),
    splitting_method = split_method,
    splitting_size_regime = splitting_size_regime,
    custom_size_split_list = list(),
    merge_split_step_vec = merge_split_step_vec,
    ms_moves_multiplier = ms_moves_multiplier,
    merge_prob_type = merge_prob_type
  )

  # add the splitting parameters
  cpp_control_list <- c(cpp_control_list, forward_kernel_params)

  t1 <- Sys.time()
  all_out <- foreach(
    chain = seq_len(runs),
    .inorder = FALSE,
    .packages = "redist"
  ) %oper%
    {
      if (chain == 1) {
        is_chain1 <- T
      }

      if (is_chain1 && !silent) {
        cat("Starting Chain ", chain, "\n", sep = "")
      }
      run_verbosity <- if (is_chain1 || !multiprocess) verbosity else 0
      t1_run <- Sys.time()

      algout <- run_redist_smc(
        nsims = nsims,
        ndists = ndists,
        total_seats = nseats,
        district_seat_sizes = seats_range,
        initial_num_regions = init_num_regions,
        adj_list = adj_list,
        counties = counties,
        pop = pop,
        step_types = step_types,
        target = pop_bounds[2],
        lower = pop_bounds[1],
        upper = pop_bounds[3],
        rho = compactness,
        sampling_space_str = sampling_space,
        control = cpp_control_list,
        constraints = constraints,
        verbosity = run_verbosity,
        diagnostic_level = diagnostic_level,
        region_id_mat = init_particles,
        region_sizes_mat = init_seats,
        log_weights = initial_log_weights
      )

      if (length(algout) == 0) {
        cli::cli_process_done()
      }

      # Make integer since arma::umat passed back to R as double
      storage.mode(algout$ancestors) <- "integer"


      diagnostic_mode = diagnostic_level == 1

      if (!diagnostic_mode) {
        # if not diagnostic mode
        # make the region_ids_mat_list input just null since there's nothing else
        algout$region_ids_mat_list <- NULL
        algout$region_seats_mat_list <- NULL
      }

      # if no merge split was run them remove those attributes
      if (!run_ms) {
        algout$merge_split_success_mat <- NULL
        algout$merge_split_attempt_counts <- NULL
      }


      # turn it into a character vector
      algout$step_split_types <- ifelse(
        algout$merge_split_steps,
        "ms",
        "smc"
      )

      num_ms_steps <- sum(
        algout$step_split_types == "ms"
      )


      # pull out the log weights
      lr <- algout$log_weights
      # get the standard deviations
      sd_lp <- algout$log_weight_stddev

      # numerically stabilize the weights
      wgt <- exp(lr - mean(lr))
      n_eff <- length(wgt) * mean(wgt)^2 / mean(wgt^2)

      if (any(is.na(lr))) {
        cli_abort(c(
          "Sampling probabilities have been corrupted.",
          "*" = "Check that none of your constraint weights are too large.
                                 The output of constraint functions multiplied by the weight
                                 should generally fall in the -5 to 5 range.",
          "*" = "If you are using custom constraints, make sure that your
                                 constraint function handles all edge cases and never returns
                                 {.val {NA}} or {.val {Inf}}",
          "*" = "If you are not using any constraints, please call
                                 {.code rlang::trace_back()} and file an issue at
                                 {.url https://github.com/alarm-redist/redist/issues/new}"
        ))
      }

      if (resample) {
        # get normalized weights for sampling
        normalized_wgts <- wgt / sum(wgt)

        # resample matrices in place using lowvar resampling
        rs_idx <- resample_plans_lowvar(
          normalized_wgts,
          algout$plans_mat,
          algout$seats,
          algout$plan_seats_saved
        )

        # rs_idx maps plan i to its new plan index
        # `rs_idx[i] = k` means you should replace plan i with plan k
        # that means if after we've resampled then the parent of plan
        # i was rs_idx[i]

        # now adjust for the resampling
        algout$ancestors <- algout$ancestors[rs_idx, , drop = FALSE]


        # add a final column for the resampling since for the resampled plans
        # plan[i] parent is rs_idx[i]
        algout$parent_index <- cbind(
          algout$parent_index,
          rs_idx[1:length(rs_idx)]
        )
        # fix storage in case converts to double for some reason
        storage.mode(algout$parent_index) <- "integer"

        # do unique parents
        nunique_parent_indices <- c(
          algout$nunique_parent_indices,
          dplyr::n_distinct(rs_idx[1:length(rs_idx)])
        )
      } else {
        nunique_parent_indices <- algout$nunique_parent_indices
      }


      t2_run <- Sys.time()


      # make sizes null if needed
      if (!algout$plan_seats_saved) {
        algout$seats <- NULL
      }


      if (!is.nan(n_eff) && n_eff / nsims <= 0.05) {
        cli_warn(c(
          "Less than 5% resampling efficiency.",
          "*" = "Increase the number of samples.",
          "*" = "Consider weakening or removing constraints.",
          "i" = "If sampling efficiency drops precipitously in the final
                                iterations, population balance is likely causing a bottleneck.
                                Try increasing {.arg pop_temper} by 0.01.",
          "i" = "If sampling efficiency declines steadily across iterations,
                                adjusting {.arg seq_alpha} upward may help a bit."
        ))
      }


      # add the numerically stable weights back
      algout$wgt <- wgt

      # flatten the region sizes pops by column into a long vector
      dim(algout$seats) <- NULL
      dim(algout$region_pops) <- NULL


      # Internal diagnostics,
      algout$internal_diagnostics <- list(
        parent_index_mat = algout$parent_index,
        log_incremental_weights_mat = algout$log_incremental_weights_mat,
        draw_tries_mat = algout$draw_tries_mat,
        tree_sizes = algout$tree_sizes,
        successful_tree_sizes = algout$successful_tree_sizes,
        parent_unsuccessful_tries_mat = algout$parent_unsuccessful_tries_mat,
        region_ids_mat_list = algout$region_ids_mat_list,
        region_seats_mat_list = algout$region_seats_mat_list,
        merge_split_success_mat = algout$merge_split_success_mat,
        forest_adjs_list = algout$forest_adjs_list,
        linking_edges_list = algout$linking_edges_list
      )

      # Information about the run
      algout$run_information <- list(
        weight_type = weight_type,
        nproc = nproc,
        ncores = ncores,
        custom_size_split_list = list(),
        valid_region_sizes_to_split_list = algout$valid_region_sizes_to_split_list,
        valid_split_region_sizes_list = algout$valid_split_region_sizes_list,
        sampling_space = sampling_space,
        split_method = split_method,
        splitting_schedule = splitting_size_regime,
        merge_split_step_vec = merge_split_step_vec,
        ms_moves_multiplier = ms_moves_multiplier,
        merge_prob_type = merge_prob_type,
        step_types = step_types,
        nsims = nsims
      )

      # add cut k for graph space
      run_forward_kernel_params <- forward_kernel_params
      if(sampling_space == GRAPH_PLAN_SPACE_SAMPLING){
          run_forward_kernel_params$cut_k_used <- algout$cut_k_vals
      }

      # add high level diagnostic stuff
      algout$l_diag <- list(
        n_eff = n_eff,
        step_n_eff = algout$step_n_eff,
        forward_kernel_params = run_forward_kernel_params,
        accept_rate = algout$acceptance_rates,
        sd_lp = sd_lp,
        unique_survive = nunique_parent_indices,
        ms_move_counts = algout$ms_move_counts,
        ancestors = algout$ancestors,
        seq_alpha = seq_alpha,
        pop_temper = pop_temper,
        runtime = as.numeric(t2_run - t1_run, units = "secs")
      )

      if (verbosity >= 1 && runs > 1) {
          cli_text(
              "Chain {chain}: {format(nsims, big.mark=',')} plans sampled in
                 {format(t2_run - t1_run, digits=2)}"
          )
      }

      algout
    }
  t2 <- Sys.time()

  if (verbosity >= 1) {
    cli_text(
      "{format(nsims*runs, big.mark=',')} plans sampled in
                 {format(t2-t1, digits=2)}"
    )
  }


  # combine if needed
  if (runs > 1) {
    plans <- do.call(cbind, lapply(all_out, function(x) x$plans))
    region_pops <- do.call(cbind, lapply(all_out, function(x) x$region_pops))
    seats <- do.call(c, lapply(all_out, function(x) x$seats))
    wgt <- do.call(c, lapply(all_out, function(x) x$wgt))
    l_diag <- lapply(all_out, function(x) x$l_diag)
    run_information <- lapply(all_out, function(x) x$run_information)
    internal_diagnostics <- lapply(all_out, function(x) x$internal_diagnostics)
  } else {
    # else if just one run extract directly
    plans <- all_out[[1]]$plans
    region_pops <- all_out[[1]]$region_pops
    seats <- all_out[[1]]$seats
    wgt <- all_out[[1]]$wgt
    l_diag <- list(all_out[[1]]$l_diag)
    run_information <- list(all_out[[1]]$run_information)
    internal_diagnostics <- list(all_out[[1]]$internal_diagnostics)
  }

  n_dist_act <- dplyr::n_distinct(plans[, 1]) # actual number (for partial plans)

  alg_type <- ifelse(run_ms, "smc_ms", "smc")


  out <- new_redist_plans(
      plans = plans,
      map = map,
      algorithm = alg_type,
    wgt =wgt,
    inputs_safe = TRUE,
    resampled = resample,
    distr_pop = region_pops,
    ndists = n_dist_act,
    seats = seats,
    n_eff = all_out[[1]]$n_eff,
    compactness = compactness,
    constraints = constraints,
    version = packageVersion("redist"),
    diagnostics = l_diag,
    run_information = run_information,
    internal_diagnostics = internal_diagnostics,
    num_admin_units = num_admin_units,
    total_runtime = t2 - t1
  )


  if (runs > 1) {
    out <- mutate(
      out,
      chain = rep(seq_len(runs), each = n_dist_act * nsims)
    ) %>%
      dplyr::relocate('chain', .after = "draw")
  }

  exist_name <- attr(map, "existing_col")
  if (!is.null(exist_name) && !isFALSE(ref_name) && ndists == final_dists) {

    ref_name <- if (!is.null(ref_name)) ref_name else exist_name
    out <- add_reference(out, map[[exist_name]], ref_name)
  }

  out
}


########
# Helper functions for `redist_smc`
########

#' Extracts splitting schedule from `split_params` parameter of `redist_smc`
#'
#'
#' @inheritParams run_redist_smc
#'
#' @returns A list with the following
#'     - `splitting_schedule`: The splitting schedule for SMC
#' @noRd
get_splitting_schedule <- function(split_params, districting_scheme){

    # setting the splitting size regime
    if ("splitting_schedule" %in% names(split_params)) {
        splitting_schedule <- split_params[["splitting_schedule"]]
        if (splitting_schedule == "split_district_only") {
            if (districting_scheme == "SMD") {
                splitting_size_regime = "split_district_only"
            } else if (districting_scheme == "MMD") {
                splitting_size_regime = "split_district_only_mmd"
            } else {
                cli_abort(
                    "Districting scheme {districting_scheme} is not supported!"
                )
            }
        } else if (splitting_schedule == "any_valid_sizes") {
            if (districting_scheme == "SMD") {
                splitting_size_regime = "any_valid_sizes"
            } else if (districting_scheme == "MMD") {
                cli_abort(
                    "Generaliezd region splits are not supported for Multi-member districting!"
                )
            } else {
                cli_abort(
                    "Districting scheme {districting_scheme} is not supported!"
                )
            }
        } else {
            cli_abort("{.arg splitting_schedule} must be either {.arg any_valid_sizes} or {.arg split_district_only}")
        }
    } else {
        # default to  district
        if (districting_scheme == "SMD") {
            splitting_size_regime = "split_district_only"
        } else if (districting_scheme == "MMD") {
            splitting_size_regime = "split_district_only_mmd"
        } else {
            cli_abort(
                "Districting scheme {districting_scheme} is not supported!"
            )
        }
    }

    return(splitting_size_regime)
}



#' Gets starting plan related parameters
#'
#'
#' @inheritParams run_redist_smc
#'
#' @returns A list with the following
#'     - `init_particles`: A 0-indexed initial plans matrix
#'     - `init_seats`: A matrix of region seat counts
#'     - `initial_log_weights`: A vector of initial log weights
#'     - `init_num_regions`: The number of initial regions
#' @noRd
get_init_plan_params <- function(
        nsims, nseats, pop, pop_bounds,
        init_particles, init_seats, init_weights
){
    # handle particle, seats, and weights inits
    if (is.null(init_particles)) {
        # if no initial plans passed in then create empty matrix
        init_particles <- matrix(0L)
        init_seats <- matrix(0L)
        init_num_regions <- 1L
    } else {
        if (inherits(init_particles, "redist_plans")) {
            if (is.null(init_seats)) {
                init_seats <- get_seats_matrix(init_particles)
            }
            if (is.null(init_weights)) {
                # get weights if not resampled, else just set all equal to 1
                init_plan_weights <- get_plans_weights(init_particles)
                if(isFALSE(attr(init_plan_weights, "resampled"))){
                    init_weights <- rep(1, nsims)
                }else{
                    init_weights <- as.vector(init_plan_weights)
                }
            }
            init_particles <- get_plans_matrix(init_particles) - 1L
        } else if (is.matrix(init_particles)) {
            if (is.null(init_seats)) {
                # else infer
                cli_warn(
                    "{.arg init_seats} was not passed in, attempting to infer number of seats per region."
                )
                init_seats <- infer_plan_seats(
                    init_particles,
                    nseats,
                    pop,
                    pop_bounds[1],
                    pop_bounds[3]
                )
            }
            # subtract 0 if needed
            if(min(init_particles[,1])){
                init_particles <- init_particles - 1L
            }
        } else {
            cli_abort(
                "{.arg init_particles} must be either a redist_plans object or a matrix"
            )
        }
        init_num_regions <- length(unique(init_particles[, 1]))
    }

    # get init weights
    if (!is.null(init_weights)) {
        # check its a vector or matrix
        if (is.matrix(init_weights) && is.numeric(init_weights)) {
            # check 1 column and nsim rows
            if (any(dim(init_weights) != c(nsims, 1))) {
                cli_abort(
                    "{.arg init_weights} must have only {nsims} elements!"
                )
            }
            # flatten
            init_weights <- as.vector(init_weights)
        } else if (is.vector(init_weights) && is.numeric(init_weights)) {
            if (length(init_weights) != nsims) {
                cli_abort(
                    "{.arg init_weights} must be of length {nsims}!"
                )
            }
        }
        # now check all positive
        if (any(init_weights <= 0)) {
            cli_abort(
                "All elements of {.arg init_weights} must be of length positive!"
            )
        }
    } else {
        init_weights <- rep(1, nsims)
    }
    initial_log_weights <- log(init_weights)

    init_plans_params <- list(
        init_particles = init_particles,
        init_seats = init_seats,
        initial_log_weights = initial_log_weights,
        init_num_regions = init_num_regions
    )

    return(init_plans_params)
}


#' Extracts arguments from `control` parameter of `redist_smc`
#'
#'
#' @inheritParams run_redist_smc
#'
#' @returns A list with the following
#'     - `nproc`: The number of parallel R processes to spawn. Defaults to 1.
#'     - `weight_type`: Must be either simple or optimal. Defaults to optimal
#' @noRd
extract_control_params <- function(control){

    control_param_names <- c("nproc", "weight_type")

    if (is.list(control) && any("nproc" %in% control_param_names)) {
        if ("nproc" %in% names(control)) {
            nproc <- control[["nproc"]]
            if (!rlang::is_integerish(nproc) || !is_scalar(nproc)) {
                cli_abort(
                    "{.arg nproc} in {.arg control} must be a single integer!"
                )
            } else if (nproc <= 0) {
                cli_abort(
                    "{.arg nproc} in {.arg control} must be a positive integer!"
                )
            }
        } else {
            # default to 1
            nproc <- 1L
        }

        if ("weight_type" %in% names(control)) {
            weight_type <- control[["weight_type"]]
            weight_type <- rlang::arg_match(
                arg = weight_type,
                values = c("optimal", "simple")
            )
        } else {
            # else default to optimal
            weight_type <- "optimal"
        }
    } else {
        nproc <- 1L
        weight_type <- "optimal"
    }

    control_params <- list(
        nproc = nproc,
        weight_type = weight_type
    )

    return(control_params)
}




#' Extracts Mergesplit paramters from `ms_params` parameter of `redist_smc`
#'
#' @inheritParams run_redist_smc
#' @param total_smc_steps How many SMC steps will be run
#'
#' @returns A list with the following
#'     - `run_ms`: Whether or not any mergesplit steps will be run
#'     - `merge_prob_type`: What probability to use when selecting pairs to merge.
#'     Defaults to "uniform".
#'     - `ms_moves_multiplier`: Multiplier to baseline number of steps
#'     - `ms_frequency`: How oftent to run mergesplit steps
#'     - `merge_split_step_vec`: vector whose length is the total number of steps
#'     and a value of true indicates that step is a mergesplit step.
#' @noRd
extract_ms_params <- function(ms_params, total_smc_steps){
    ms_param_names <- c("ms_moves_multiplier", "ms_frequency", "merge_prob_type")

    # create merge split parameter information
    if (is.list(ms_params) && any(ms_param_names %in% names(ms_params))) {
        run_ms <- TRUE
        # check if ms_moves_multiplier was passed else default to 1
        if ("ms_moves_multiplier" %in% names(ms_params)) {
            ms_moves_multiplier <- ms_params[["ms_moves_multiplier"]]
            # check that ms_moves_multiplier is positive
            if (
                !is_scalar(ms_moves_multiplier) || !ms_moves_multiplier > 0
            ) {
                cli_abort("{.arg ms_moves_multiplier} must be a positive scalar")
            }
        } else {
            ms_moves_multiplier <- 1L
        }

        # check if the frequency was passed else default to after every step
        if ("ms_frequency" %in% names(ms_params)) {
            ms_frequency <- ms_params[["ms_frequency"]]
        } else {
            # else default to after every step
            ms_frequency <- 1L
        }

        # check merge probability
        if ("merge_prob_type" %in% names(ms_params)) {
            merge_prob_type <- ms_params[["merge_prob_type"]]
            if (
                !is_scalar(merge_prob_type) || merge_prob_type != "uniform"
            ) {
                cli_abort("Only uniform merge probability is supported right now!")
            }
        } else {
            # else default to after every step
            merge_prob_type <- "uniform"
        }
    } else {
        run_ms <- FALSE
        merge_prob_type <- "ignore"
        ms_moves_multiplier <- NULL
        ms_frequency <- NULL
    }

    if (!run_ms) {
        merge_split_step_vec <- rep(FALSE, total_smc_steps)
    } else if (ms_frequency == 1) {
        # if frequency 1 then do after every step
        merge_split_step_vec <- rep(FALSE, total_smc_steps)
        # Now add merge split every `ms_frequency` steps
        # insertion trick
        # https://stackoverflow.com/questions/1493969/insert-elements-into-a-vector-at-given-indexes
        ind <- seq(from = ms_frequency, to = total_smc_steps, by = ms_frequency)
        val <- c(merge_split_step_vec, rep(TRUE, length(ind)))
        id <- c(seq_along(merge_split_step_vec), ind + 0.5)

        # number of merge split is sum of trues
        merge_split_step_vec <- val[order(id)]
    } else if (ms_frequency == -1) {
        # if negative 1 then just put at the end
        merge_split_step_vec <- rep(FALSE, total_smc_steps)
        merge_split_step_vec <- c(
            merge_split_step_vec,
            TRUE
        )
    }


    extracted_ms_params <- list(
        run_ms = run_ms,
        merge_prob_type = merge_prob_type,
        ms_moves_multiplier = ms_moves_multiplier,
        ms_frequency = ms_frequency,
        merge_split_step_vec = merge_split_step_vec
    )

    return(extracted_ms_params)
}
