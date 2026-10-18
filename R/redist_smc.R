##############################################
## Author: Philip O'Sullivan
## Institution: Harvard University
## Date Created: 2025/02/15
## Purpose: Wrapper for running smc cpp code
##############################################

#' Generalized SMCS Redistricting Sampler (O'Sullivan, Imai, and McCartan (2026) 
#' <doi:10.48550/arXiv.2603.22188>)
#'
#'
#' `redist_smc` uses a Sequential Monte Carlo Sampler algorithm
#' (O'Sullivan, Imai, and McCartan (2026) 
#' <doi:10.48550/arXiv.2603.22188>) to generate representative samples of
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
#' manually. The resampling is done after the last SMC split, so if the run
#' ends with mergesplit rounds they are run on the resampled plans, which can
#' restore some of the diversity lost in resampling.
#' @param runs How many independent parallel runs to conduct. Each run will
#' have `nsims` simulations. Multiple runs allows for estimation of simulation
#' standard errors. Output will only be shown for the first run. For
#' compatibility with MCMC methods, runs are identified with the `chain`
#' column in the output.
#' @param ncores How many threads to use to parallelize plan generation within each
#' process. The default value is 1 which is single-threading. If `ncores` is set
#' to 0 it will use the number of available cores on the machine
#' as long as `nsims` and the number of units is large enough.
#' If more than one core is used, the
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
#' @param init_weights A vector of length `nsims` of plan weights associated
#' with `init_particles`. The weights must all be strictly positive and are
#' normalized to sum to 1, so they may be passed in on any scale. If no weights
#' are passed in then the plans are weighted equally, unless `init_particles` is
#' a [redist_plans] object that was not resampled, in which case its weights
#' are used. If `init_particles` is a matrix of plans that were not resampled
#' then it is recommended to pass their weights in.
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
#' @param ms_params A list of mergesplit parameters. Mergesplit rounds are run
#' if any of the parameters below are passed. For maps with 10 or more districts
#' they are also run by default, with the defaults below, unless `ms_params` is
#' `NULL`. Pass `ms_params = NULL` to turn mergesplit off.
#' \itemize{
#'  \item \code{frequency}: How often to perform mergesplit rounds after SMC
#' steps. Should either be an integer (meaning run after every _ smc steps)
#' or a vector of 1 indexed step numbers indicating after which smc steps to
#' run a mergesplit round. A value of -1 means just run a single mergesplit round
#' after all smc steps have been run. A value of 1 means run after every smc step.
#' Defaults to 1.
#' \item \code{mh_accept_per_smc} Multiplier to the baseline number of mergesplit
#' steps to be performed each round. For a mergesplit round the baseline number of
#' steps is calculated as the ceiling of 1 over the previous mergesplit round
#' acceptance rate (or smc step if no prior mergesplit rounds were done). The
#' total number of steps is `ceiling(mh_accept_per_smc * baseline_num_steps)`.
#' This should be thought of as specifying how many successful mergesplit steps
#' are desired for each mergesplit round. Defaults to the number of districts
#' divided by 3, rounded down, for maps with 10 or more districts and 1
#' otherwise.
#' \item \code{pair_rule} What probability to use to select regions to merge
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
#'  incremental weights and parent indices.
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
#'  \item \code{max_split_tries} The maximum number of times to try splitting
#' plans before giving up. If more than `max_split_tries` attempts at splitting
#' are made without a success then an exception will be thrown. This can be
#' useful for maps and constraints where it might become impossible to split
#' plans after a certain step. Defaults to not checking.
#'  \item \code{est_norm_unbiased} Whether or not to perform extra work
#'  allowing for unbiased estimation of the normalizing constant of the target
#'  distribution. This entails at each step sampling an extra plan, saving the
#'  number of tries before sampling the extra plan, and then discarding the
#'  extra plan. The number of attempts can be used with the weights to form an
#'  unbiased estimator, as opposed to the biased normalizing
#'  constant estimators that are obtained with the weights alone.
#'  \item \code{return_augmented_samples} Whether or not to return the spanning
#'  forests and linking edges associated with the sampled plans (if applicable).
#'  Setting this value to `TRUE` can substantially increase the size of the
#'  `redist_plans` object returned.
#'  \item \code{keep_init_plan_info} If `init_particles` is a [redist_plans]
#'  object sampled with `redist_smc()`, whether to combine its diagnostics with
#'  those of this run, so the output has the diagnostics of every step as if
#'  all the steps had been run at once. Requires `runs = 1` both here and for
#'  the initial plans, the same `nsims`, and initial plans sampled with
#'  `resample = FALSE`. Defaults to `FALSE`.
#'  \item \code{md_alpha} What power to use for multidistrict selection
#'  probability. When choosing a multidistrict to split, one will be chosen with
#'  probability proportional to <size>^md_alpha.
#' }
#' @param adapt_k_thresh Deprecated. Pass in through the `split_params` arg.
#' The threshold value used in the heuristic to select a
#' value `k_i` for each splitting iteration. Higher values are more accurate
#' but may require more computation. Set to 1 for the most conservative
#' sampling. Must be between 0 and 1.
#' @param truncate Deprecated, do not use anymore.
#' @param trunc_fn Deprecated, do not use anymore.
#' @param final_infl Deprecated, do not use anymore.
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
    ncores = 1L,
    init_particles = NULL,
    init_seats = NULL,
    init_weights = NULL,
    sampling_space = c("linking_edge", "spanning_forest", "graph_plan"),
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
    adapt_k_thresh = 0.99,
    truncate = (compactness != 1),
    trunc_fn = redist_quantile_trunc,
    final_infl = 1
) {
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
    # deprecated arguments
    if (!missing(truncate)) {
        cli::cli_abort("{.arg truncate} is deprecated. Do not use this argument.")
    }
    if (!missing(trunc_fn)) {
        cli::cli_abort("{.arg truncate} is deprecated. Do not use this argument.")
    }
    if (!missing(final_infl)) {
        cli::cli_abort("{.arg truncate} is deprecated. Do not use this argument.")
    }

    if (!rlang::is_scalar_atomic(compactness) || compactness < 0) {
        cli::cli_abort("{.arg compactness} must be non-negative.")
    }
    if (seq_alpha <= 0 || seq_alpha > 1 || !rlang::is_scalar_atomic(seq_alpha)) {
        cli::cli_abort("{.arg seq_alpha} must lie in (0, 1].")
    }
    if (nsims < 1) {
        cli::cli_abort("{.arg nsims} must be positive.")
    }
    if (!rlang::is_scalar_logical(resample)) {
        cli::cli_abort("{.arg resample} must be a single {.val TRUE} or {.val FALSE}.")
    }

    # check default inputs
    diagnostics <- rlang::arg_match(diagnostics)
    diagnostic_level <- dplyr::case_when(
        diagnostics == "basic" ~ 0,
        diagnostics == "all" ~ 1,
        .default = 0
    )

    # validate constraints
    constraints <- validate_constraints(map, rlang::enquo(constraints), compactness)
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
    if (!setequal(seats_range, seq.int(from = min(seats_range), to = max(seats_range)))) {
        cli::cli_abort(
      "For {.arg seats_range} only a continuous range of district seat sizes are allowed!"
    )
    }

    # get the splitting size regime
    splitting_size_regime <- get_splitting_schedule(
        split_params,
        districting_scheme
    )

    # get initial plan parameters
    # kept so its diagnostics can be combined with this run's if requested
    init_plans_object <- init_particles
    initial_plan_params <- get_init_plan_params(
        nsims,
        nseats,
        pop,
        pop_bounds,
        init_particles,
        init_seats,
        init_weights
    )
    init_particles <- initial_plan_params$init_particles
    init_seats <- initial_plan_params$init_seats
    initial_weights <- initial_plan_params$initial_weights
    init_num_regions <- initial_plan_params$init_num_regions

    if (is.null(n_steps)) {
        n_steps <- ndists - init_num_regions
    }
    final_dists <- init_num_regions + n_steps
    if (final_dists > ndists) {
        cli::cli_abort("Too many districts already drawn to take {n_steps} steps.")
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
    # if (num_admin_units > 1 && sampling_space == LINKING_EDGE_SPACE_SAMPLING) {
    #     cli::cli_abort(
    #         "Linking Edge Sampling with counties is not supported right now"
    #     )
    # }

    total_smc_steps <- n_steps

    # get merge split parameter information
    ms_params_list <- extract_ms_params(ms_params, total_smc_steps, ndists)
    run_ms <- ms_params_list$run_ms
    pair_rule <- ms_params_list$pair_rule
    mh_accept_per_smc <- ms_params_list$mh_accept_per_smc
    ms_frequency <- ms_params_list$ms_frequency
    merge_split_step_vec <- ms_params_list$merge_split_step_vec

    # get the step types
    step_types <- ifelse(merge_split_step_vec, "ms", "smc")
    if (sum(!merge_split_step_vec) != total_smc_steps) {
        cli::cli_abort(
      "In {.field step_types} the number of smc steps must be equal to {.field total_smc_steps}"
    )
    }
    # assert first step is not smc
    if (merge_split_step_vec[1]) {
        cli::cli_abort(
      "The first step cannot be mergesplit! An SMC step must be run first!"
    )
    }
    total_ms_steps <- sum(merge_split_step_vec)
    # total number of steps to run
    total_steps <- total_smc_steps + total_ms_steps

    # if run_ms is true then seq_alpha must be 1
    if (run_ms && seq_alpha != 1) {
        seq_alpha <- 1L
        cli::cli_warn(
      "{.arg seq_alpha} must be set to 1 if any mergesplit steps are being run!"
    )
    }

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
                !rlang::is_scalar_integerish(ncores)
        ) {
            cli::cli_abort("{.arg ncores} must be a single integer!")
        } else if (ncores == 0) {
            ncores <- ncores_max
        } else if (ncores < 0) {
            cli::cli_abort(
        "{.arg ncores} can't be negative!"
      )
        }
    } else {
        ncores <- ncores_max
    }

    # Now extract control parameters
    control_params_list <- extract_control_params(control, compactness)
    nproc <- control_params_list[["nproc"]]
    weight_type <- control_params_list[["weight_type"]]
    cache_weights <- control_params_list[["cache_weights"]]
    max_split_tries <- control_params_list[["max_split_tries"]]
    est_norm_unbiased <- control_params_list[["est_norm_unbiased"]]
    return_augmented_samples <- control_params_list[["return_augmented_samples"]]

    if (return_augmented_samples && sampling_space == GRAPH_PLAN_SPACE_SAMPLING) {
        cli::cli_abort(c(
            "{.arg return_augmented_samples} is not available for graph space sampling.",
            "i" = "Only spanning forest and linking edge space plans have augmented samples."
        ))
    }

    # The diagnostics of the initial plans, combined with this run's at the end
    # so the output looks like a single run through all the steps
    init_plans_info <- NULL
    if (control_params_list[["keep_init_plan_info"]] &&
        inherits(init_plans_object, "redist_plans")) {
        if (!attr(init_plans_object, "algorithm") %in% c("smc", "smc_ms")) {
            cli::cli_abort(
                "{.arg keep_init_plan_info} requires initial plans sampled with {.fn redist_smc}."
            )
        }
        if (runs != 1 || length(attr(init_plans_object, "diagnostics")) != 1) {
            cli::cli_abort(
                "{.arg keep_init_plan_info} is only supported when both the initial plans
                and this call use {.code runs = 1}."
            )
        }
        if (!isFALSE(attr(init_plans_object, "resampled"))) {
            cli::cli_abort(c(
                "Combining diagnostics with {.arg keep_init_plan_info} is not supported when
                the partial plans have already been resampled.",
                "i" = "Sample the initial plans with {.code resample = FALSE}."
            ))
        }
        init_internal_diagnostics <- attr(init_plans_object, "internal_diagnostics")[[1]]
        if (nrow(init_internal_diagnostics$parent_index_mat) != nsims) {
            cli::cli_abort(
                "{.arg keep_init_plan_info} requires {.arg nsims} to equal the number of
                plans the initial plans were sampled with."
            )
        }
        init_plans_info <- list(
            diagnostics = attr(init_plans_object, "diagnostics")[[1]],
            run_information = attr(init_plans_object, "run_information")[[1]],
            internal_diagnostics = init_internal_diagnostics,
            total_runtime = attr(init_plans_object, "total_runtime")
        )
    }

    multiprocess <- nproc > 1
    # make sure we're not spawning more proccesses than runs
    if (multiprocess) {
        nproc <- min(runs, nproc)
    }

    # warn if more processes than cores
    if (nproc > ncores_max) {
        cli::cli_warn(
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
                outfile = nullfile(),
                methods = FALSE,
                useXDR = .Platform$endian != "little"
            )
        }
        # Makes it so only one process will print but if more runs then processes
        # it doesn't just print once
        parallel::clusterEvalQ(cl, {
            if (!exists("is_chain1", envir = .GlobalEnv)) {
                is_chain1 <- FALSE
            }
            NULL
        })
        doParallel::registerDoParallel(cl, cores = nproc)
        init_workers(cl)
        on.exit(stopCluster(cl))
        if (!silent) cat("Spawning ", nproc, " clusters \n")
    } else {
        `%oper%` <- `%do%`
    }

    cpp_control_list <- list(
    weight_type = weight_type,
    cache_weights = cache_weights,
    max_split_tries = max_split_tries,
    est_norm_unbiased = est_norm_unbiased,
    return_augmented_samples = return_augmented_samples,
    resample = resample,
    seq_alpha = seq_alpha,
    pop_temper = pop_temper,
    num_threads = as.integer(ncores),
    splitting_method = split_method,
    splitting_size_regime = splitting_size_regime,
    custom_size_split_list = list(),
    merge_split_step_vec = merge_split_step_vec,
    mh_accept_per_smc = mh_accept_per_smc,
    pair_rule = pair_rule
  )

    # add md alpha if included
    if ("md_alpha" %in% names(control_params_list)){
        md_alpha <- control_params_list[["md_alpha"]]
        cpp_control_list <- c(cpp_control_list,
                              list(md_alpha = md_alpha))
    }


    # add the splitting parameters
    cpp_control_list <- c(cpp_control_list, forward_kernel_params)

    # normalized initial plan weights, if not passed in all plans start with
    # the same weight
    if (!is.null(initial_weights)) {
        cpp_control_list$initial_weights <- initial_weights
    }

    # When continuing an earlier run, pass in what it would have carried into
    # this run's first step if it had kept going: its weights are split with
    # seq_alpha, its last estimated cut k starts the next estimate, and its last
    # mergesplit acceptance rate sets the size of the first mergesplit round
    if (!is.null(init_plans_info)) {
        cpp_control_list$continue_from_init_plans <- TRUE
        init_step_types <- init_plans_info$run_information$step_types
        init_cut_k <- init_plans_info$diagnostics$forward_kernel_params$cut_k_used[
            init_step_types == "smc"
        ]
        init_cut_k <- init_cut_k[init_cut_k > 0]
        if (length(init_cut_k) > 0) {
            cpp_control_list$initial_last_cut_k <- as.integer(utils::tail(init_cut_k, 1))
        }
        init_ms_accept_rates <- init_plans_info$diagnostics$accept_rate[
            init_step_types == "ms"
        ]
        if (length(init_ms_accept_rates) > 0) {
            cpp_control_list$initial_ms_acceptance_rate <- utils::tail(init_ms_accept_rates, 1)
        }
    }

    # Initialize the plans with their saved spanning forests (and linking
    # edges) if the sampling space can use them, else they are drawn at random.
    # Spanning forest space ignores any linking edges.
    init_augmented_samples <- initial_plan_params$init_augmented_samples
    if (!is.null(init_augmented_samples)) {
        if (sampling_space == FOREST_SPACE_SAMPLING) {
            cpp_control_list$initial_spanning_forests <- init_augmented_samples$spanning_forests
        } else if (sampling_space == LINKING_EDGE_SPACE_SAMPLING &&
                   !is.null(init_augmented_samples$linking_edges)) {
            cpp_control_list$initial_spanning_forests <- init_augmented_samples$spanning_forests
            cpp_control_list$initial_linking_edges <- init_augmented_samples$linking_edges
        }
    }

    t1 <- Sys.time()
    all_out <- foreach(
        chain = seq_len(runs),
        .inorder = FALSE,
        .errorhandling = "stop"
    ) %oper%
        {
            if (chain == 1) {
                is_chain1 <- TRUE
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
                region_sizes_mat = init_seats
            )

            if (length(algout) == 0) {
                # means user cancelled or something went wrong
                cli::cli_abort("Process Interrupted!")
            }

            diagnostic_mode <- diagnostic_level == 1

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

            # The order the steps were done in, including the final resampling
            # which is right after the last split. The sampler stores the
            # diagnostics indexed by every step with the resampling last, so
            # those are put in this order.
            all_step_types <- build_all_step_types(algout$step_split_types, resample)
            step_order <- get_step_order(all_step_types)
            if (length(algout$nunique_plans) == length(step_order)) {
                algout$nunique_plans <- algout$nunique_plans[step_order]
            }
            # the per step diagnostic lists are also named (smc_step1,
            # smc_step2, ms_step1, ..., resample) so the plans, forests, and
            # linking edges from a step can be matched by name or index
            step_names <- get_step_list_names(all_step_types)
            for (step_list in c("region_ids_mat_list", "forest_adjs_list", "linking_edges_list",
                                "region_seats_mat_list")) {
                if (length(algout[[step_list]]) == length(step_order)) {
                    algout[step_list] <- list(algout[[step_list]][step_order])
                    if (step_list != "region_seats_mat_list") {
                        names(algout[[step_list]]) <- step_names
                    }
                }
            }

            # pull out the log weights
            lr <- algout$log_weights
            # get the standard deviations
            sd_lp <- algout$log_weight_stddev

            # numerically stabilize the weights
            wgt <- exp(lr - mean(lr))
            n_eff <- length(wgt) * mean(wgt)^2 / mean(wgt^2)

            if (anyNA(lr)) {
                cli::cli_abort(c(
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

            # Resampling is done inside the sampler if needed. If performed it
            # also sizes the ancestry diagnostics to match, appending the resampling
            # index as a final parent index column and a final entry to the
            # unique parent and unique plan counts.
            nunique_parent_indices <- algout$nunique_parent_indices
            nunique_plans <- algout$nunique_plans

            t2_run <- Sys.time()

            # make sizes null if needed
            if (!algout$plan_seats_saved) {
                algout$seats <- NULL
            }

            if (!is.nan(n_eff) && n_eff / nsims <= 0.05) {
                cli::cli_warn(c(
          "Less than 5% resampling efficiency.",
          "*" = "Increase the number of samples.",
          "*" = "Increase `mh_accept_per_smc` or start adding mergesplit steps if not already using them.",
          "*" = "Consider weakening or removing constraints if sampling efficiency is low."
        ))
            }

            # add the numerically stable weights back
            algout$wgt <- wgt

            # flatten the region sizes pops by column into a long vector
            dim(algout$seats) <- NULL
            dim(algout$region_pops) <- NULL

            # add the time information
            time_breakdowns <- list(
                smc_step_parameter_estimation_times = algout$smc_step_parameter_estimation_times,
                smc_split_times = algout$smc_split_times,
                smc_weight_times = algout$smc_weight_times,
                ms_step_parameter_estimation_times = algout$ms_step_parameter_estimation_times,
                ms_step_times = algout$ms_step_times
            )

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
        # `ifelse()` is vectorized and shapes its result like the condition, so
        # `ifelse(FALSE, x, NULL)` collapses to NULL and the assignment fails.
        tries_before_extra_particle = if (est_norm_unbiased) {
            algout$tries_before_extra_particle
        } else {
            NULL
        },
        log_blank_map_target_density = ifelse(
            init_num_regions == 1,
            algout$log_blank_map_target_density,
            NA),
        time_breakdowns = time_breakdowns,
        granular_times = algout$granular_times,
        forest_adjs_list = algout$forest_adjs_list,
        linking_edges_list = algout$linking_edges_list
      )




            # Information about the run
            algout$run_information <- list(
        weight_type = weight_type,
        nproc = nproc,
        ncores = ncores,
        custom_size_split_list = list(),
        sampling_space = sampling_space,
        split_method = split_method,
        splitting_schedule = splitting_size_regime,
        merge_split_step_vec = merge_split_step_vec,
        mh_accept_per_smc = mh_accept_per_smc,
        pair_rule = pair_rule,
        step_types = step_types,
        nsims = nsims,
        # if resampled and the run ends in mergesplit rounds, the resampling
        # was done after the last split, before those rounds
        # This was added after the realization that performing MCMC moves after
        # the resample as opposed to before is never worse.
        resample_before_final_ms = TRUE,
        # the steps in the order they were done, including any resampling
        all_step_types = all_step_types
      )

            # add cut k for graph space
            run_forward_kernel_params <- forward_kernel_params
            if (sampling_space == GRAPH_PLAN_SPACE_SAMPLING) {
                run_forward_kernel_params$cut_k_used <- algout$cut_k_vals
            }
            # add multidistrict selection parameters
            if (splitting_size_regime == "any_valid_sizes"){
                run_forward_kernel_params$multidistrict_selection_alpha <- algout$multidistrict_selection_alpha
            }

            # add high level diagnostic stuff
            algout$l_diag <- list(
        n_eff = n_eff,
        step_n_eff = algout$step_n_eff,
        forward_kernel_params = run_forward_kernel_params,
        accept_rate = algout$acceptance_rates,
        sd_lp = sd_lp,
        unique_survive = nunique_parent_indices,
        nunique_plans = nunique_plans,
        ms_step_counts = algout$ms_step_counts,
        seq_alpha = seq_alpha,
        pop_temper = pop_temper,
        runtime = as.numeric(t2_run - t1_run, units = "secs")
      )

            if (run_verbosity >= 1 && runs > 1) {
                cli::cli_text(
          "Chain {chain}: {format(nsims, big.mark=',')} plans sampled in
                 {format(t2_run - t1_run, digits=2)}"
        )
            }

            algout
        }
    t2 <- Sys.time()

    if (verbosity >= 1) {
        cli::cli_text(
      "{format(nsims*runs, big.mark=',')} plans sampled in
                 {format(t2-t1, digits=2)}"
    )
    }

    # combine the initial plans' diagnostics with this run's (runs is 1 here)
    if (!is.null(init_plans_info) && !is.null(all_out[[1]])) {
        all_out[[1]] <- combine_init_plan_info(init_plans_info, all_out[[1]])
    }

    bad_runs <- sapply(all_out, is.null)

    # combine if needed
    if (runs > 1) {
        if (all(bad_runs)) {
            return(NULL)
        } else if (any(bad_runs)) {
            all_out <- all_out[!bad_runs]
            runs <- length(all_out)
        }
        plans <- do.call(cbind, lapply(all_out, function(x) x$plans))
        region_pops <- do.call(cbind, lapply(all_out, function(x) x$region_pops))
        seats <- do.call(c, lapply(all_out, function(x) x$seats))
        wgt <- do.call(c, lapply(all_out, function(x) x$wgt))
        l_diag <- lapply(all_out, function(x) x$l_diag)
        run_information <- lapply(all_out, function(x) x$run_information)
        internal_diagnostics <- lapply(all_out, function(x) x$internal_diagnostics)
    } else {
        # if null then return null
        if (all(bad_runs)) {
            return(NULL)
        }

        # else if just one run extract directly
        plans <- all_out[[1]]$plans
        region_pops <- all_out[[1]]$region_pops
        seats <- all_out[[1]]$seats
        wgt <- all_out[[1]]$wgt
        l_diag <- list(all_out[[1]]$l_diag)
        run_information <- list(all_out[[1]]$run_information)
        internal_diagnostics <- list(all_out[[1]]$internal_diagnostics)
    }

    # Spanning forests and linking edges of the final plans, indexed by the
    # sampled draws (the plans matrix columns that are not reference plans, see
    # `subset_augmented_samples`). The linking edge `draw` column is the sampled
    # draw so each run is offset by the draws of the runs before it.
    augmented_samples <- NULL
    if (return_augmented_samples) {
        linking_edges <- NULL
        if (sampling_space == LINKING_EDGE_SPACE_SAMPLING) {
            linking_edges <- do.call(rbind, lapply(seq_along(all_out), function(i) {
                run_linking_edges <- all_out[[i]]$linking_edges
                run_linking_edges$draw <- run_linking_edges$draw + (i - 1L) * nsims
                run_linking_edges
            }))
        }
        augmented_samples <- list(
            edge_list = graph_edge_list(adj_list),
            spanning_forests = do.call(cbind, lapply(all_out, function(x) x$spanning_forests)),
            linking_edges = linking_edges
        )
    }

    n_dist_act <- dplyr::n_distinct(plans[, 1]) # actual number (for partial plans)

    alg_type <- ifelse(run_ms, "smc_ms", "smc")

    out <- new_redist_plans(
        plans = plans,
        map = map,
        algorithm = alg_type,
        wgt = wgt,
        inputs_safe = TRUE,
        resampled = resample,
        distr_pop = region_pops,
        ndists = n_dist_act,
        nseats = nseats,
        seats = seats,
        n_eff = all_out[[1]]$n_eff,
        compactness = compactness,
        constraints = constraints,
        counties = counties,
        version = packageVersion("redist"),
        diagnostics = l_diag,
        run_information = run_information,
        internal_diagnostics = internal_diagnostics,
        augmented_samples = augmented_samples,
        num_admin_units = num_admin_units,
        # includes the initial plans' sampling time if their diagnostics were kept
        total_runtime = as.numeric(t2 - t1, units = "secs") +
            sum(init_plans_info$total_runtime)
    )

    if (runs > 1) {
        out <- mutate(
            out,
            chain = rep(seq_len(runs), each = n_dist_act * nsims)
        ) |>
            dplyr::relocate("chain", .after = "draw")
    }

    exist_name <- attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name) && ndists == final_dists) {
        ref_name <- if (!is.null(ref_name)) ref_name else exist_name

        ref_plan_list <- get_ref_plan_and_seats(map)

        out <- add_reference(
            plans = out,
            ref_plan = ref_plan_list$ref_plan,
            name = ref_name,
            ref_seats = ref_plan_list$ref_seats
        )
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
get_splitting_schedule <- function(split_params, districting_scheme) {
    # setting the splitting size regime
    if ("splitting_schedule" %in% names(split_params)) {
        splitting_schedule <- split_params[["splitting_schedule"]]
        if (splitting_schedule == "split_district_only") {
            if (districting_scheme == "single") {
                splitting_size_regime <- "split_district_only"
            } else if (districting_scheme == "multiple") {
                splitting_size_regime <- "split_district_only_mmd"
            } else {
                cli::cli_abort(
          "Districting scheme {districting_scheme} is not supported!"
        )
            }
        } else if (splitting_schedule == "any_valid_sizes") {
            if (districting_scheme == "single") {
                splitting_size_regime <- "any_valid_sizes"
            } else if (districting_scheme == "multiple") {
                cli::cli_abort(
          "Generaliezd region splits are not supported for Multi-member districting!"
        )
            } else {
                cli::cli_abort(
          "Districting scheme {districting_scheme} is not supported!"
        )
            }
        } else {
            cli::cli_abort(
        "{.arg splitting_schedule} must be either {.arg any_valid_sizes} or {.arg split_district_only}"
      )
        }
    } else {
        # default to  any_valid_sizes for single member
        if (districting_scheme == "single") {
            splitting_size_regime <- "any_valid_sizes"
        } else if (districting_scheme == "multiple") {
            splitting_size_regime <- "split_district_only_mmd"
        } else {
            cli::cli_abort(
        "Districting scheme {districting_scheme} is not supported!"
      )
        }
    }

    splitting_size_regime
}


#' Names of the per step diagnostic lists
#'
#' The steps are named by type and counted separately (`smc_step1`,
#' `smc_step2`, `ms_step1`, ...), with the final resampling named `resample`.
#'
#' @param all_step_types Character vector of `"smc"`, `"ms"`, and
#' `"resample"` steps in the order they were done, see `get_all_step_types()`
#'
#' @returns A character vector of names
#' @noRd
get_step_list_names <- function(all_step_types) {
    step_type_counts <- stats::ave(
        seq_along(all_step_types), all_step_types, FUN = seq_along
    )
    ifelse(
        all_step_types == "resample",
        "resample",
        paste0(all_step_types, "_step", step_type_counts)
    )
}


#' The order the steps of a new SMC run are done in
#'
#' The final resampling, if done, is right after the last split. So if the run
#' ends in mergesplit rounds they are run on the resampled plans.
#'
#' @param step_types Character vector of `"smc"` and `"ms"` step types
#' @param resampled Whether the plans were resampled
#'
#' @returns `step_types` with `"resample"` added after the last split if
#' resampled
#' @noRd
build_all_step_types <- function(step_types, resampled) {
    if (!resampled) {
        return(step_types)
    }
    append(step_types, "resample", after = max(which(step_types == "smc")))
}


#' The order the steps of an SMC run were done in
#'
#' Runs store this as `all_step_types` in their `run_information`. Older runs
#' that didn't store it resampled after every other step, so it is inferred as
#' their step types followed by `"resample"` if they were resampled.
#'
#' @param run_info The run's `run_information`
#' @param resampled Whether the plans were resampled
#'
#' @returns A character vector of `"smc"`, `"ms"`, and `"resample"` steps in
#' the order they were done
#' @noRd
get_all_step_types <- function(run_info, resampled) {
    if (!is.null(run_info$all_step_types)) {
        return(run_info$all_step_types)
    }
    c(run_info$step_types, if (resampled) "resample")
}


#' Indices putting per step values into the order the steps were done
#'
#' The sampler stores values indexed by every step with the final resampling's
#' value last. Indexing them with this puts them in the order the steps were
#' done in.
#'
#' @param all_step_types The order the steps were done in, see
#' `get_all_step_types()`
#'
#' @returns An integer vector of indices
#' @noRd
get_step_order <- function(all_step_types) {
    is_resample <- all_step_types == "resample"
    ifelse(is_resample, length(all_step_types), cumsum(!is_resample))
}


#' Combines the diagnostics of initial plans with a run started from them
#'
#' Used for `keep_init_plan_info` so that the output of a run started from
#' partial plans has the diagnostics of every step, as if it had been sampled
#' in a single run. The per step diagnostics of the initial plans' run come
#' first, followed by this run's. The initial plans must not have been
#' resampled, so their run has no final resampling step.
#'
#' @param init_info A list of the `diagnostics`, `run_information`, and
#' `internal_diagnostics` of the initial plans' run
#' @param algout The output of the run started from the initial plans after
#' its `l_diag`, `run_information`, and `internal_diagnostics` are created
#'
#' @returns `algout` with the combined diagnostics
#' @noRd
combine_init_plan_info <- function(init_info, algout) {
    init_diag <- init_info$diagnostics
    init_run_info <- init_info$run_information
    init_internal <- init_info$internal_diagnostics

    # Merges two lists of timings element by element, so per step vectors are
    # concatenated, per plan by step matrices are bound, and nested lists are
    # merged the same way. Anything else, like flags, is taken from `new`.
    combine_by_name <- function(old, new) {
        for (name in names(new)) {
            old_val <- old[[name]]
            new_val <- new[[name]]
            if (is.null(old_val)) next
            if (is.list(new_val)) {
                new[[name]] <- combine_by_name(old_val, new_val)
            } else if (is.matrix(new_val)) {
                new[[name]] <- cbind(old_val, new_val)
            } else if (is.numeric(new_val)) {
                new[[name]] <- c(old_val, new_val)
            }
        }
        new
    }

    # names for this run's per step list entries continue the initial run's.
    # The initial plans weren't resampled so their steps are their step types
    step_types <- c(init_run_info$step_types, algout$run_information$step_types)
    init_all_step_types <- get_all_step_types(init_run_info, FALSE)
    all_step_types <- c(init_all_step_types, algout$run_information$all_step_types)
    new_step_names <- get_step_list_names(all_step_types)[-seq_along(init_all_step_types)]
    combine_step_list <- function(old, new) {
        if (length(new) == length(new_step_names)) {
            names(new) <- new_step_names
        }
        c(old, new)
    }

    # high level diagnostics
    diag <- algout$l_diag
    for (name in c("step_n_eff", "accept_rate", "sd_lp", "unique_survive",
                   "nunique_plans", "ms_step_counts")) {
        diag[name] <- list(c(init_diag[[name]], diag[[name]]))
    }
    if (!is.null(diag$forward_kernel_params$cut_k_used)) {
        diag$forward_kernel_params$cut_k_used <- c(
            init_diag$forward_kernel_params$cut_k_used,
            diag$forward_kernel_params$cut_k_used
        )
    }
    diag$runtime <- init_diag$runtime + diag$runtime

    # information about the steps run
    run_info <- algout$run_information
    run_info$step_types <- step_types
    run_info$all_step_types <- all_step_types
    run_info$merge_split_step_vec <- c(
        init_run_info$merge_split_step_vec, run_info$merge_split_step_vec
    )

    # internal diagnostics, `[<-` with a list keeps NULL entries in place
    internal <- algout$internal_diagnostics
    for (name in c("parent_index_mat", "log_incremental_weights_mat", "draw_tries_mat",
                   "tree_sizes", "successful_tree_sizes", "parent_unsuccessful_tries_mat",
                   "merge_split_success_mat")) {
        internal[name] <- list(cbind(init_internal[[name]], internal[[name]]))
    }
    # only usable if tracked for every step
    internal["tries_before_extra_particle"] <- list(
        if (is.null(init_internal$tries_before_extra_particle) ||
            is.null(internal$tries_before_extra_particle)) {
            NULL
        } else {
            c(init_internal$tries_before_extra_particle, internal$tries_before_extra_particle)
        }
    )
    # only the run that started from a blank map can compute it
    internal$log_blank_map_target_density <- init_internal$log_blank_map_target_density
    internal$time_breakdowns <- combine_by_name(
        init_internal$time_breakdowns, internal$time_breakdowns
    )
    internal$granular_times <- combine_by_name(
        init_internal$granular_times, internal$granular_times
    )
    for (name in c("region_ids_mat_list", "forest_adjs_list", "linking_edges_list")) {
        internal[name] <- list(combine_step_list(init_internal[[name]], internal[[name]]))
    }
    internal["region_seats_mat_list"] <- list(
        c(init_internal$region_seats_mat_list, internal$region_seats_mat_list)
    )

    algout$l_diag <- diag
    algout$run_information <- run_info
    algout$internal_diagnostics <- internal
    algout
}


#' Gets starting plan related parameters
#'
#'
#' @inheritParams run_redist_smc
#'
#' @returns A list with the following
#'     - `init_particles`: A 0-indexed initial plans matrix
#'     - `init_seats`: A matrix of region seat counts
#'     - `initial_weights`: A vector of initial weights normalized to sum to 1,
#'     or `NULL` if all plans start with the same weight
#'     - `init_num_regions`: The number of initial regions
#'     - `init_augmented_samples`: The saved spanning forests and linking edges
#'     of `init_particles` if it is a [redist_plans] object that has them, else
#'     `NULL`
#' @noRd
get_init_plan_params <- function(
    nsims,
    nseats,
    pop,
    pop_bounds,
    init_particles,
    init_seats,
    init_weights
) {
    # saved spanning forests and linking edges of the initial plans, if any
    init_augmented_samples <- NULL

    # handle particle, seats, and weights inits
    if (is.null(init_particles)) {
        # if no initial plans passed in then create empty matrix
        init_particles <- matrix(0L)
        init_seats <- matrix(0L)
        init_num_regions <- 1L
    } else {
        if (inherits(init_particles, "redist_plans")) {
            init_particles <- subset_sampled(init_particles)
            if (is.null(init_seats)) {
                init_seats <- get_seats_matrix(init_particles)
            }
            if (is.null(init_weights)) {
                # get weights if not resampled, else leave them equal. The
                # stored weights of resampled plans are the weights from before
                # resampling so they don't correspond to the plans.
                init_plan_weights <- get_plans_weights(init_particles)
                if (isFALSE(attr(init_plan_weights, "resampled"))) {
                    init_weights <- as.vector(init_plan_weights)
                }
            }
            init_augmented_samples <- attr(init_particles, "augmented_samples")
            init_particles <- get_plans_matrix(init_particles)
        } else if (is.matrix(init_particles)) {
            if (is.null(init_seats)) {
                # else infer
                cli::cli_warn(
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
        } else {
            cli::cli_abort(
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
                cli::cli_abort(
          "{.arg init_weights} must have only {nsims} elements!"
        )
            }
            # flatten
            init_weights <- as.vector(init_weights)
        } else if (
            (is.vector(init_weights) && is.numeric(init_weights)) &&
                (length(init_weights) != nsims)
        ) {
            cli::cli_abort(
          "{.arg init_weights} must be of length {nsims}!"
        )
        }
        # now check all positive
        if (any(init_weights <= 0)) {
            cli::cli_abort(
        "All elements of {.arg init_weights} must be of length positive!"
      )
        }

        if(sum(init_weights) != 1){
            initial_weights <- init_weights / sum(init_weights)
        }else{
            initial_weights <- init_weights
        }
    } else {
        # all plans start with the same weight
        initial_weights <- NULL
    }

    init_plans_params <- list(
    init_particles = init_particles,
    init_seats = init_seats,
    initial_weights = initial_weights,
    init_num_regions = init_num_regions,
    init_augmented_samples = init_augmented_samples
  )

    init_plans_params
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
extract_control_params <- function(control, compactness) {
    control_param_names <- c("nproc", "weight_type",
                             "cache_weights", "max_split_tries",
                             "est_norm_unbiased", "return_augmented_samples",
                             "keep_init_plan_info")

    default_nproc <- 1L
    default_weight_type <- "optimal"
    default_cache_weights <- compactness != 1
    default_max_split_tries <- 0L
    default_est_norm_unbiased <- TRUE

    if ("nproc" %in% names(control)) {
        nproc <- control[["nproc"]]
        if (!rlang::is_scalar_integerish(nproc)) {
            cli::cli_abort(
          "{.arg nproc} in {.arg control} must be a single integer!"
        )
        } else if (nproc <= 0) {
            cli::cli_abort(
          "{.arg nproc} in {.arg control} must be a positive integer!"
        )
        }
    } else {
        # default value
        nproc <- default_nproc
    }

    if ("weight_type" %in% names(control)) {
        weight_type <- control[["weight_type"]]
        weight_type <- rlang::arg_match(
            arg = weight_type,
            values = c("optimal", "simple")
        )
    } else {
        # else default to optimal
        weight_type <- default_weight_type
    }

    if ("cache_weights" %in% names(control)) {
        cache_weights <- control[["cache_weights"]]
        if (!rlang::is_scalar_logical(cache_weights)) {
            cli::cli_abort("{.arg cache_weights} must be a scalar boolean")
        }
    } else {
        # else default
        cache_weights <- default_cache_weights
    }

    if ("max_split_tries" %in% names(control)) {
        max_split_tries <- control[["max_split_tries"]]
        if (!rlang::is_scalar_integerish(max_split_tries)) {
            cli::cli_abort(
              "{.arg max_split_tries} in {.arg control} must be a single integer!"
          )
        } else if (max_split_tries <= 0) {
            cli::cli_abort(
              "{.arg max_split_tries} in {.arg control} must be a positive integer!"
          )
        }
    } else {
        # default to 0 meaning don't check
        max_split_tries <- default_max_split_tries
    }

    if ("est_norm_unbiased" %in% names(control)) {
        est_norm_unbiased <- control[["est_norm_unbiased"]]
        if (!rlang::is_scalar_logical(est_norm_unbiased)) {
            cli::cli_abort("{.arg est_norm_unbiased} must be a scalar boolean")
        }
    } else {
        # else default
        est_norm_unbiased <- default_est_norm_unbiased
    }

    if ("return_augmented_samples" %in% names(control)) {
        return_augmented_samples <- control[["return_augmented_samples"]]
        if (!rlang::is_scalar_logical(return_augmented_samples) ||
            is.na(return_augmented_samples)) {
            cli::cli_abort("{.arg return_augmented_samples} must be a scalar boolean")
        }
    } else {
        return_augmented_samples <- FALSE
    }

    if ("keep_init_plan_info" %in% names(control)) {
        keep_init_plan_info <- control[["keep_init_plan_info"]]
        if (!rlang::is_scalar_logical(keep_init_plan_info) || is.na(keep_init_plan_info)) {
            cli::cli_abort("{.arg keep_init_plan_info} must be a scalar boolean")
        }
    } else {
        keep_init_plan_info <- FALSE
    }

    control_params <- list(
    nproc = nproc,
    weight_type = weight_type,
    cache_weights = cache_weights,
    max_split_tries = max_split_tries,
    est_norm_unbiased = est_norm_unbiased,
    return_augmented_samples = return_augmented_samples,
    keep_init_plan_info = keep_init_plan_info
  )

    if ("md_alpha" %in% names(control)){
        md_alpha <- control[["md_alpha"]]
        if(
            !rlang::is_scalar_atomic(md_alpha) ||
            !is.numeric(md_alpha)
        ){
            cli::cli_abort("{.arg md_alpha} must be a numeric scalar")
        }
        control_params <- c(
            control_params,
            list(md_alpha = md_alpha)
        )
    }

    control_params
}


#' Extracts Mergesplit paramters from `ms_params` parameter of `redist_smc`
#'
#' @inheritParams run_redist_smc
#' @param total_smc_steps How many SMC steps will be run
#' @param ndists The number of districts in the map, which sets the defaults
#'
#' @returns A list with the following
#'     - `run_ms`: Whether or not any mergesplit steps will be run
#'     - `pair_rule`: What probability to use when selecting pairs to merge.
#'     Defaults to "uniform".
#'     - `mh_accept_per_smc`: Multiplier to baseline number of steps
#'     - `ms_frequency`: How oftent to run mergesplit steps
#'     - `merge_split_step_vec`: vector whose length is the total number of steps
#'     and a value of true indicates that step is a mergesplit step.
#' @noRd
extract_ms_params <- function(ms_params, total_smc_steps, ndists) {
    ms_param_names <- c("mh_accept_per_smc", "frequency", "pair_rule")
    # maps with this many districts run mergesplit by default
    default_ms_min_ndists <- 10

    # Mergesplit runs if any mergesplit parameters are passed or, for maps with
    # enough districts, by default. Passing `NULL` turns it off.
    run_ms <- is.list(ms_params) && (
        any(ms_param_names %in% names(ms_params)) || ndists >= default_ms_min_ndists
    )

    # create merge split parameter information
    if (run_ms) {
        # check if mh_accept_per_smc was passed else default to a third of the
        # districts for large maps and 1 otherwise
        if ("mh_accept_per_smc" %in% names(ms_params)) {
            mh_accept_per_smc <- ms_params[["mh_accept_per_smc"]]
            # check that mh_accept_per_smc is positive
            if (
                !rlang::is_scalar_atomic(mh_accept_per_smc) ||
                    !mh_accept_per_smc > 0
            ) {
                cli::cli_abort("{.arg mh_accept_per_smc} must be a positive scalar")
            }
        } else if (ndists >= default_ms_min_ndists) {
            mh_accept_per_smc <- floor(ndists / 3)
        } else {
            mh_accept_per_smc <- 1L
        }

        # check if the frequency was passed else default to after every step
        if ("frequency" %in% names(ms_params)) {
            ms_frequency <- ms_params[["frequency"]]
            # ensure its integers
            if (!rlang::is_integerish(ms_frequency)) {
                cli::cli_abort("{.arg frequency} must be a integer valued")
            }
        } else {
            # else default to after every step
            ms_frequency <- 1L
        }

        # check merge probability
        if ("pair_rule" %in% names(ms_params)) {
            pair_rule <- ms_params[["pair_rule"]]
            if (
                !rlang::is_scalar_character(pair_rule) ||
                    pair_rule != "uniform"
            ) {
                cli::cli_abort("Only uniform merge probability is supported right now!")
            }
        } else {
            # else default to after every step
            pair_rule <- "uniform"
        }
    } else {
        run_ms <- FALSE
        pair_rule <- "ignore"
        mh_accept_per_smc <- NULL
        ms_frequency <- NULL
    }

    if (!run_ms) {
        merge_split_step_vec <- rep(FALSE, total_smc_steps)
    } else if (rlang::is_scalar_integerish(ms_frequency)) {
        # check if its a scalar
        if (ms_frequency >= 1) {
            # then we do merge split after every ms_frequency-th step.
            ms_frequency <- min(ms_frequency, total_smc_steps)

            merge_split_step_vec <- rep(FALSE, total_smc_steps)
            insert_positions <- seq(ms_frequency, total_smc_steps, by = ms_frequency)

            # Build the result
            offset <- 0

            for (pos in insert_positions) {
                idx <- pos + offset
                merge_split_step_vec <- append(merge_split_step_vec, TRUE, after = idx)
                offset <- offset + 1
            }
        } else if (ms_frequency <= -1) {
            # if its a negative number then that means the last ms_frequency steps
            # make sure its absolute value is less than the number of smc steps
            if (abs(ms_frequency) > total_smc_steps) {
                ms_frequency <- -total_smc_steps
            }

            merge_split_step_vec <- rep(FALSE, total_smc_steps)
            # This makes it so we run mergesplit after the last abs(ms_frequency)
            # SMC steps
            freq <- abs(ms_frequency)
            stopifnot(freq <= length(merge_split_step_vec))
            merge_split_step_vec <- c(
              rep(FALSE, length(merge_split_step_vec) - freq),
              rep(c(FALSE, TRUE), times = freq)
          )
        } else {
            cli::cli_abort("{.arg frequency} cannot be 0")
        }
    } else {
        # else that means its a vector of specific SMC steps to run it after
        cli::cli_abort("Specific step {.arg frequency} is not implemented yet!")
    }

    extracted_ms_params <- list(
    run_ms = run_ms,
    pair_rule = pair_rule,
    mh_accept_per_smc = mh_accept_per_smc,
    ms_frequency = ms_frequency,
    merge_split_step_vec = merge_split_step_vec
  )

    extracted_ms_params
}


#' Shared implementation of the SMC normalizing constant estimators
#'
#' Both estimators divide the summed incremental weights at each SMC step by
#' the number of spanning trees that were drawn to get there. They differ only
#' in whether the extra particle that is drawn and discarded counts toward that
#' denominator: including it makes the estimator unbiased, leaving it out gives
#' the cheaper biased version.
#'
#' @param plans A `redist_plans` object.
#' @param include_extra_tries Whether to add `tries_before_extra_particle` to
#'   the attempt counts.
#'
#' @returns One log normalizing constant estimate per run.
#' @noRd
est_log_norm_const <- function(plans, include_extra_tries) {
    if (!inherits(plans, "redist_plans")) {
        cli::cli_abort("{.arg plans} must be a {.cls redist_plans} type!")
    }
    if (attr(plans, "version") < "5.0") {
        cli::cli_abort("{.arg plans} must have been sampled with Redist 5.0 or later!")
    }
    # the estimators assume the weights are fully resampled at every step
    seq_alphas <- vapply(
        attr(plans, "diagnostics"), function(x) x$seq_alpha %||% 1, numeric(1)
    )
    if (any(seq_alphas != 1)) {
        cli::cli_abort(c(
            "Normalizing constant estimation is only supported for plans sampled
            with {.code seq_alpha = 1}.",
            "x" = "{.arg plans} were sampled with {.arg seq_alpha} =
                   {.val {unique(seq_alphas)}}."
        ))
    }

    internal_diagnostics <- attr(plans, "internal_diagnostics")

    if (include_extra_tries) {
        # Check the extra particle counts are actually there. Testing for the
        # name is not enough: sampling with `est_norm_unbiased = FALSE` leaves
        # the name in place with a NULL value, which would otherwise surface
        # much later as an opaque length mismatch.
        extra_info_added <- sapply(
            internal_diagnostics,
            function(x) !is.null(x$tries_before_extra_particle)
        ) |>
            all()
        if (!extra_info_added) {
            cli::cli_abort(c(
                "{.arg plans} must have been sampled with a version of Redist that
                 tracked {.arg tries_before_extra_particle} in
                 {.arg internal_diagnostics}",
                "i" = "Sample with {.code control = list(est_norm_unbiased = TRUE)},
                       or use {.fn est_norm_biased}, which does not need it."
            ))
        }
    }

    # get initial weight
    initial_log_weights <- sapply(
        internal_diagnostics,
        function(x) x$log_blank_map_target_density
    )
    # Warn if they are not all the same
    if (!all(abs(initial_log_weights[1] - initial_log_weights) < 1e-14)) {
        cli::cli_warn("{.arg log_blank_map_target_density} in {.arg internal_diagnostics}
                       are not equal accross runs. These plans may have been sampled on
                      different maps or with different values for {.arg counties}")
    }

    wgt_mats <- lapply(
        internal_diagnostics,
        function(x) x$log_incremental_weights_mat
    )

    # extract the estimated k value if needed
    sampling_space <- attr(plans, "run_information")[[1]]$sampling_space

    smc_steps <- lapply(
        attr(plans, "run_information"),
        function(a_run_list) a_run_list$step_types == "smc"
    )

    # log of the number of attempts made at each step, optionally including the
    # extra particle that was drawn and thrown away
    log_attempt_counts <- Map(
        function(idiag, smc_step_vec) {
            tries <- colSums(idiag$draw_tries_mat)[smc_step_vec]
            if (include_extra_tries) {
                tries <- tries + idiag$tries_before_extra_particle
            }
            log(tries)
        },
        internal_diagnostics,
        smc_steps
    )

    if (sampling_space == GRAPH_PLAN_SPACE_SAMPLING) {
        # get the k values
        param_mat <- lapply(
            seq_along(attr(plans, "diagnostics")),
            function(i) {
                attr(plans, "diagnostics")[[i]]$forward_kernel_params$cut_k_used[smc_steps[[i]]] |>
                    log()
            }
        )
    } else {
        # else just do nothing
        param_mat <- lapply(
            smc_steps,
            function(x) rep(0L, sum(x))
        )
    }

    # adds the log parameter value to the log weights
    # equivalent to weights * parameter value and
    # then exponentiates and sums and takes log
    log_of_wgt_step_sums <- Map(
        function(W, p) {
            stopifnot(ncol(W) == nrow(p))
            sweep(W, MARGIN = 2, STATS = p, FUN = "+") |>
                exp() |>
                colSums() |>
                log()
        },
        wgt_mats,
        param_mat
    )

    # now for each smc step we want to divide the sum of the weights
    # by the number of attempts and take their product
    log_norm_estimates <- mapply(
        function(W, tries, initital_log_weight) {
            stopifnot(length(W) == length(tries))
            initital_log_weight + sum(W - tries)
        },
        log_of_wgt_step_sums,
        log_attempt_counts,
        lapply(initial_log_weights, function(x) x)
    )

    return(log_norm_estimates)
}


#' Estimates the log of normalizing constant of plans sampled with `redist_smc`
#'
#' Returns a biased estimate of the log normalizing constant for plans sampled
#' using `redist_smc`. The denominator at each step is the number of spanning
#' trees drawn for the particles that were kept, which leaves out the extra
#' discarded particle that [est_norm_unbiased()] uses to remove the bias.
#'
#' Unlike [est_norm_unbiased()] this needs no extra bookkeeping from the
#' sampler, so it works on any plans object from `Redist 5.0` or later that was
#' sampled with `seq_alpha = 1`. Plans sampled with any other `seq_alpha` give
#' an error.
#'
#' @param plans A `redist_plans` object generated using `Redist 5.0` or later.
#' @inheritParams redist_smc
#'
#' @return An estimate of the normalizing constant from each run
#'
#' @md
#' @order 1
#' @export
est_norm_biased <- function(map, plans) {
    est_log_norm_const(plans, include_extra_tries = FALSE)
}


#' Unbiased Estimates the log of normalizing constant of plans sampled with `redist_smc`
#'
#' Returns an unbiased estimate of the log normalizing constant for plans
#' sampled using `redist_smc`. Requires that the plans were sampled with
#' `control = list(est_norm_unbiased = TRUE)`, which draws one extra particle
#' per step and discards it so that the attempt counts give an unbiased
#' estimate. The plans must also have been sampled with `seq_alpha = 1`; plans
#' sampled with any other `seq_alpha` give an error.
#'
#' @param plans A `redist_plans` object generated using `Redist 5.1` or later.
#' @inheritParams redist_smc
#'
#' @return An estimate of the normalizing constant from each run
#'
#' @md
#' @order 1
#' @export
est_norm_unbiased <- function(map, plans) {
    est_log_norm_const(plans, include_extra_tries = TRUE)
}


#' Deprecated Helper function to truncate importance weights
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
