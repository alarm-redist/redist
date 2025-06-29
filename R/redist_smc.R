##############################################
## Author: Philip O'Sullivan
## Institution: Harvard University
## Date Created: 2025/02/15
## Purpose: Wrapper for running smc cpp code
##############################################

DEBUG_MODE <- FALSE


#' Generalized SMCS Redistricting Sampler (O'Sullivan, McCartan and Imai ???)
#'
#'
#' `redist_smc` uses a Sequential Monte Carlo Sampler algorithm
#' (O'Sullivan, McCartan and Imai ???) to generate representative samples of
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
#' importance sampling weights.  In these cases, these weights are by default
#' truncated using [redist_quantile_trunc()] to stabilize the resulting
#' estimates, but if truncation is used, a specific truncation function should
#' probably be chosen by the user.
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
#' @param init_nseats A matrix of the number of seats per region of
#' the partial plans to begin sampling from. For advanced use only. The matrix
#' must have `nsims` columns, as many rows as there are regions in `init_particles`
#' and each column must sum to the total number of seats in the map. Not needed
#' if `init_particles` is a [redist_plans] object. If `init_particles` are passed
#' but not `init_nseats` then the number of seats will be attempted
#' to be inferred.
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
#' Must be between 0 and 1.
#' @param truncate Whether to truncate the importance sampling weights at the
#' final step by `trunc_fn`.  Recommended if `compactness` is not 1.
#' Truncation only applied if `resample=TRUE`.
#' @param trunc_fn A function which takes in a vector of weights and returns a
#' truncated vector. If the [loo][loo::loo] package is installed (strongly
#' recommended), will default to Pareto-smoothed Importance Sampling (PSIS)
#' rather than naive truncation.
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
#' computationally expensive, especially for computationally complex constraints.
#' }
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
#' # One run with multiple processes
#' redist_smc(fl_map, 1000, nproc = 2)
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_smc <- function(
        map, nsims, counties = NULL,
        compactness = 1, constraints = list(),
        resample = TRUE,
        runs = 1L, ncores = 0L,
        init_particles = NULL,
        init_nseats = NULL,
        sampling_space = c("graph_plan", "spanning_forest", "linking_edge"),
        split_method = c("naive_top_k","uniform_valid_edge", "expo_bigger_abs_dev"),
        split_params = list(adapt_k_thresh = .99),
        ms_params = list(),
        n_steps = NULL, seq_alpha = 1L,
        truncate = (compactness != 1), trunc_fn = redist_quantile_trunc,
        pop_temper = 0, ref_name = NULL,
        verbose = FALSE, silent = FALSE, diagnostics = c("basic", "all"),
        control = list(weight_type = "optimal", nproc = 1L)
)
{

    # defunct code not used right now
    custom_size_split_list = list()

    # not supported right now
    if(truncate){
        cli::cli_abort("Truncation not suppored right now!")
    }

    # check default inputs
    sampling_space <- rlang::arg_match(sampling_space)
    split_method <- rlang::arg_match(split_method)
    diagnostics <- rlang::arg_match(diagnostics)
    diagnostic_level <- dplyr::case_when(
        diagnostics == "basic" ~ 0,
        diagnostics == "all" ~ 1,
        .default = 0
    )


    # come up with a better name for diagnostic level code



    # want things to be as similar as possible to current code (dev branch)
    #   - better to add function inputs, bad to remove old ones
    #   - For mergesplit parallel put it back and then just have it call the
    #        other on

    if (!assertthat::is.scalar(compactness) || compactness < 0)
        cli_abort("{.arg compactness} must be non-negative.")
    if (seq_alpha <= 0 || seq_alpha > 1 || !assertthat::is.scalar(seq_alpha))
        cli_abort("{.arg seq_alpha} must lie in (0, 1].")

    # if graph space default to k stuff
    if(sampling_space == GRAPH_PLAN_SPACE_SAMPLING){
        if(missing(split_method)){
            split_method <- NAIVE_K_SPLITTING
        }
        if(missing(split_params)){
            split_params = list(
                adapt_k_thresh=.99
            )
        }
    }else if(sampling_space == FOREST_SPACE_SAMPLING || sampling_space == LINKING_EDGE_SPACE_SAMPLING){
        # the others default to uniform
        if(missing(split_method)){
            split_method <- UNIF_VALID_EDGE_SPLITTING
        }
    }

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
    ndists <- attr(map, "ndists")
    total_seats <- attr(map, "total_seats")
    district_seat_sizes <- attr(map, "district_seat_sizes")
    storage.mode(district_seat_sizes) <- "integer"
    districting_scheme <- attr(map, "districting_scheme")

    # check that the seat sizes are a range
    if(any(district_seat_sizes != seq.int(from = min(district_seat_sizes), to = max(district_seat_sizes) ))){
        cli::cli_abort("For {.arg district_seat_sizes} only a continuous range of district seat sizes are allowed!")
    }

    # setting the splitting size regime
    if("splitting_schedule" %in% names(split_params)){
        splitting_schedule <- split_params[["splitting_schedule"]]
        if(splitting_schedule == "split_district_only"){
            if(districting_scheme == "SMD"){
                splitting_size_regime = "split_district_only"
            }else if(districting_scheme == "MMD"){
                splitting_size_regime = "split_district_only_mmd"
            }else{
                cli::cli_abort("Districting scheme {districting_scheme} is not supported!")
            }
        }else if(splitting_schedule == "any_valid_sizes"){
            if(districting_scheme == "SMD"){
                splitting_size_regime = "any_valid_sizes"
            }else if(districting_scheme == "MMD"){
                cli::cli_abort("Generaliezd region splits are not supported for Multi-member districting!")
            }else{
                cli::cli_abort("Districting scheme {districting_scheme} is not supported!")
            }
        }else{ # else its custom
            cli::cli_abort("Custom splitting schedules are not supported right now!")
            # only support doing a single size right now
            # validate it
            validate_custom_size_split_list(ndists, n_steps, init_num_regions, custom_size_split_list)
            splitting_size_regime = "one_custom_size"
        }
    }else{
        # default to any size if SMD and district if MMD
        if(districting_scheme == "SMD"){
            splitting_size_regime = "any_valid_sizes"
        }else if(districting_scheme == "MMD"){
            splitting_size_regime = "split_district_only_mmd"
        }else{
            cli::cli_abort("Districting scheme {districting_scheme} is not supported!")
        }
    }


    # handle particle inits
    if (is.null(init_particles)) {
        # if no initial plans passed in then create empty matrix
        init_particles <- matrix(0L)
        init_nseats <- matrix(0L)
        init_num_regions <- 1L
    }else {
        if (inherits(init_particles, "redist_plans")){
            if(is.null(init_nseats)){
                init_nseats <- get_nseats_matrix(init_particles)
            }
            init_particles <- get_plans_matrix(init_particles) - 1L

        }
        if(is.null(init_nseats)){
            # else infer
            cli::cli_warn("{.arg init_nseats} was not passed in, attempting to infer number of seats per region.")
            init_nseats <- infer_plan_nseats(
                init_particles, total_seats, pop,
                pop_bounds[1], pop_bounds[3]
            )
        }
        # if user input then check its valid
        init_num_regions <- length(unique(init_particles[,1]))
        # validate_initial_region_id_mat(init_particles, V, nsims, init_num_regions)
        # validate_initial_region_sizes_mat(init_nseats, ndists, nsims, init_num_regions)
    }
    if (is.null(n_steps)) {
        n_steps <- ndists - init_num_regions
    }
    final_dists <- init_num_regions + n_steps
    if (final_dists > ndists) {
        cli_abort("Too many districts already drawn to take {n_steps} steps.")
    }


    #validate the splitting method and params
    split_params <- validate_sample_space_and_splitting_method(
        sampling_space, split_method, split_params, n_steps
    )
    total_smc_steps <- n_steps



    ms_param_names <- c("ms_moves_multiplier", "ms_frequency", "merge_prob_type")

    # create merge split parameter information
    if(is.list(ms_params) && any(ms_param_names %in% names(ms_params))){
        run_ms <- TRUE
        # check if ms_moves_multiplier was passed else default to 1
        if("ms_moves_multiplier" %in% names(ms_params)){
            ms_moves_multiplier <- ms_params[["ms_moves_multiplier"]]
            # check that ms_moves_multiplier is positive
            if(!assertthat::is.scalar(ms_moves_multiplier) || !ms_moves_multiplier > 0){
                cli::cli_abort("{.arg ms_moves_multiplier} must be a positive scalar")
            }
        }else{
            ms_moves_multiplier <- 1L
        }

        # check if the frequency was passed else default to after every step
        if("ms_frequency" %in% names(ms_params)){
            ms_frequency <- ms_params[["ms_frequency"]]
        }else{
            # else default to after every step
            ms_frequency <- 1L
        }

        # check merge probability
        if("merge_prob_type" %in% names(ms_params)){
            merge_prob_type <- ms_params[["merge_prob_type"]]
            if(!assertthat::is.scalar(merge_prob_type) || merge_prob_type != "uniform"){
                cli::cli_abort("Only uniform merge probability is supported right now!")
            }
        }else{
            # else default to after every step
            merge_prob_type <- "uniform"
        }

    }else{
        run_ms <- FALSE
        merge_prob_type <- "ignore"
        ms_moves_multiplier <- NULL
        ms_frequency <- NULL
    }


    if(!run_ms){
        merge_split_step_vec <- rep(FALSE, total_smc_steps)
    }else if(ms_frequency == 1){
        # if frequency 1 then do after every step
        merge_split_step_vec <- rep(FALSE, total_smc_steps)
        # Now add merge split every `ms_frequency` steps
        # insertion trick
        # https://stackoverflow.com/questions/1493969/insert-elements-into-a-vector-at-given-indexes
        ind <- seq(from = ms_frequency, to = total_smc_steps, by = ms_frequency)
        val <- c( merge_split_step_vec, rep(TRUE,length(ind)) )
        id  <- c( seq_along(merge_split_step_vec), ind+0.5 )

        # number of merge split is sum of trues
        merge_split_step_vec <- val[order(id)]
    }else if(ms_frequency == -1){
        # if negative 1 then just put at the end
        merge_split_step_vec <- rep(FALSE, total_smc_steps)
        merge_split_step_vec <- c(
            merge_split_step_vec, TRUE
        )
    }




    # get the step types
    step_types <- ifelse(merge_split_step_vec, "ms", "smc")
    assertthat::assert_that(sum(!merge_split_step_vec) == total_smc_steps)
    # assert first step is not smc
    assertthat::assert_that(!merge_split_step_vec[1])
    total_ms_steps <- sum(merge_split_step_vec)
    # total number of steps to run
    total_steps <- total_smc_steps + total_ms_steps
    any_ms_steps_ran <- run_ms



    # compute lags thing
    lags <- 1 + unique(round((ndists - 1)^0.8*seq(0, 0.7, length.out = 4)^0.9))

    # verbosity stuff
    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0

    # set up parallel processing stuff
    ncores_max <- parallel::detectCores()

    # legacy, ncores is essentially the number of threads per process
    num_threads_per_process <- ncores
    if(!is.null(num_threads_per_process)){
        if(!rlang::is_integerish(num_threads_per_process) || !assertthat::is.scalar(num_threads_per_process)){
            cli::cli_abort("{.arg num_threads_per_process} in {.arg control} must be a single integer!")
        }else if(num_threads_per_process == 0){
            num_threads_per_process <- ncores_max
        }else if(num_threads_per_process < 0){
            cli::cli_abort("{.arg num_threads_per_process} in {.arg control} can't be negative!")
        }
    }else{
        num_threads_per_process <- ncores_max
    }



    # Now handle control parameters
    control_param_names <- c("nproc", "weight_type")

    if(is.list(control) && any("nproc" %in% control_param_names)){
        if("nproc" %in% names(control)){
            nproc <- control[["nproc"]]
            if(!rlang::is_integerish(nproc) || !assertthat::is.scalar(nproc)){
                cli::cli_abort("{.arg nproc} in {.arg control} must be a single integer!")
            }else if(nproc <= 0){
                cli::cli_abort("{.arg nproc} in {.arg control} must be a positive integer!")
            }
        }else{
            # default to 1
            nproc <- 1L
        }

        if("weight_type" %in% names(control)){
            weight_type <- control[["weight_type"]]
            weight_type <- rlang::arg_match(
                arg = weight_type,
                values = c("optimal", "simple")
                )

        }else{
            # else default to optimal
            weight_type <- "optimal"
        }
    }else{
        nproc <- 1L
        weight_type <- "optimal"
    }

    multiprocess <- nproc > 1
    # make sure we're not spawning more proccesses than runs
    nproc <- min(runs, nproc)

    # warn if more processes than cores
    if(nproc > ncores_max){
        cli_warn("Inputted number of processes to spawn is greater than detected number of cores on machine")
    }

    nproc <- as.integer(nproc)
    num_threads_per_process <- as.integer(num_threads_per_process)


    if (nproc > 1 && runs > 1) {
        `%oper%` <- `%dorng%`
        of <- if (Sys.info()[["sysname"]] == "Windows") {
            tempfile(pattern = paste0("smc_", substr(Sys.time(), 1, 10)), fileext = ".txt")
        } else {
            ""
        }

        # this makes a cluster using socket (NOT FORK) with
        if (!silent)
            cl <- makeCluster(nproc, outfile = of, methods = FALSE,
                              useXDR = .Platform$endian != "little")
        else
            cl <- makeCluster(nproc, methods = FALSE,
                              useXDR = .Platform$endian != "little")
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
        cat("Spawning " , nproc, " clusters \n")
    } else {
        `%oper%` <- `%do%`
    }


    cpp_control_list <- list(
        weight_type=weight_type,
        lags=lags,
        seq_alpha=seq_alpha,
        pop_temper = pop_temper,
        num_threads=as.integer(num_threads_per_process),
        splitting_method = split_method,
        splitting_size_regime = splitting_size_regime,
        custom_size_split_list=custom_size_split_list,
        merge_split_step_vec = merge_split_step_vec,
        ms_moves_multiplier = ms_moves_multiplier,
        merge_prob_type = merge_prob_type
    )

    # add the splitting parameters
    cpp_control_list <- c(cpp_control_list, split_params)


    t1 <- Sys.time()
    all_out <- foreach(chain = seq_len(runs), .inorder = FALSE) %oper% {
        if(chain == 1){
            is_chain1 <- T
        }

        if(is_chain1 && !silent) cat("Starting Chain ", chain, "\n", sep = "")
        run_verbosity <- if (is_chain1 || !multiprocess) verbosity else 0
        t1_run <- Sys.time()

        algout <- redist::run_redist_gsmc(
            nsims=nsims,
            ndists=ndists, total_seats=total_seats,
            district_seat_sizes = district_seat_sizes,
            initial_num_regions=init_num_regions,
            adj_list=adj_list,
            counties=counties,
            pop=pop,
            step_types=step_types,
            target=pop_bounds[2],
            lower=pop_bounds[1],
            upper=pop_bounds[3],
            rho=compactness,
            sampling_space = sampling_space,
            control = cpp_control_list,
            constraints = constraints,
            verbosity=run_verbosity,
            diagnostic_level=diagnostic_level,
            region_id_mat = init_particles,
            region_sizes_mat = init_nseats
        )

        if (length(algout) == 0) {
            cli::cli_process_done()
        }
        if(DEBUG_MODE) print("Done with c++ Code!")


        diagnostic_mode = diagnostic_level == 1

        if(!diagnostic_mode){
            # if not diagnostic mode
            # make the region_ids_mat_list input just null since there's nothing else
            algout$region_ids_mat_list <- NULL
            algout$region_sizes_mat_list <- NULL
        }else{
            # make intermediate plans 1 indexed
            for (i in seq_len(length(algout$region_ids_mat_list))) {
                algout$region_ids_mat_list[[i]] <- algout$region_ids_mat_list[[i]] + 1L
            }
            # add plans as well

        }
        if(DEBUG_MODE) print("Checkpoint 1!")

        # if no merge split was run them remove those attributes
        if(!any_ms_steps_ran){
            algout$merge_split_success_mat <- NULL
            algout$merge_split_attempt_counts <- NULL
        }

        gc()
        if(DEBUG_MODE) print("Checkpoint 2 - gc!")

        # turn it into a character vector
        algout$step_split_types <- ifelse(
            algout$merge_split_steps, "ms", "smc"
        )

        num_ms_steps <- sum(
            algout$step_split_types == "ms"
        )

        # make parent succesful tries matrix counting the number of
        # times a parent index was successfully sampled
        # NOTE: Not storing to save space
        # parent_successful_tries_mat <- apply(
        #     algout$parent_index, 2, tabulate, nbins = nsims
        # )


        # pull out the log weights
        lr <- algout$log_weights

        wgt <- exp(lr - mean(lr))
        n_eff <- length(wgt)*mean(wgt)^2/mean(wgt^2)

        if(DEBUG_MODE) print("Checkpoint 3 - weight and lr!")
        if (any(is.na(lr))) {
            cli_abort(c("Sampling probabilities have been corrupted.",
                        "*" = "Check that none of your constraint weights are too large.
                                 The output of constraint functions multiplied by the weight
                                 should generally fall in the -5 to 5 range.",
                        "*" = "If you are using custom constraints, make sure that your
                                 constraint function handles all edge cases and never returns
                                 {.val {NA}} or {.val {Inf}}",
                        "*" = "If you are not using any constraints, please call
                                 {.code rlang::trace_back()} and file an issue at
                                 {.url https://github.com/alarm-redist/redist/issues/new}"))
        }

        if (resample) {

            if (!truncate) {
                normalized_wgts <- wgt/sum(wgt)
            } else if (requireNamespace("loo", quietly = TRUE) && is.null(trunc_fn)) {
                cli::cli_abort("loo truncation not suppored right now!")
                normalized_wgts <- wgt/sum(wgt)
                truncated_
                normalized_wgts <- loo::weights.importance_sampling(
                    loo::psis(log(mod_wgt), r_eff = NA), log = FALSE)
            } else {
                # truncate the weights
                wgt <- trunc_fn(wgt)
                # get normalized weights
                normalized_wgts <- wgt/sum(wgt)
            }

            n_eff <- 1/sum(normalized_wgts^2)

            # resample matrices in place
            rs_idx <- resample_plans_lowvar(
                normalized_wgts,
                algout$plans_mat,
                algout$region_sizes_mat, algout$plan_sizes_saved
            )
            if(DEBUG_MODE) print("Checkpoint 3.5 - did in place reordering!")

            n_unique <- dplyr::n_distinct(rs_idx)
            # now adjust for the resampling
            algout$ancestors <- algout$ancestors[rs_idx, , drop = FALSE]

            # add a final column for the resampling
            # NOTE: I THINK THIS IS WRONG, MIGHT NEED TO FLIP COLUMN
            algout$parent_index <- cbind(algout$parent_index, rs_idx[1:length(rs_idx)])

            # do unique parents
            nunique_parent_indices <- c(
                algout$nunique_parent_indices,
                dplyr::n_distinct(rs_idx[1:length(rs_idx)]))

            #TODO probably need to adjust the rest of these as well
            storage.mode(algout$ancestors) <- "integer"
        }else{
            nunique_parent_indices <- algout$nunique_parent_indices
        }
        if(DEBUG_MODE) print("Checkpoint 4 - after resample!")

        t2_run <- Sys.time()
        # get original ancestor matrix from parent index
        algout$original_ancestors_mat <- get_original_ancestors_mat(
            algout$parent_index
        )
        if(DEBUG_MODE) print("Checkpoint 5 - after og anvestor mat!")


        # now for the smc step only diagnostics make it so
        # the merge split steps are just NA
        dummy_vec <- rep(NA, length(algout$merge_split_steps))

        # do effective sample size
        dummy_vec[!algout$merge_split_steps] <- algout$step_n_eff
        algout$step_n_eff <- dummy_vec
        # do log weight sd
        dummy_vec[!algout$merge_split_steps] <- algout$log_weight_stddev
        sd_lp <- c(dummy_vec, sd(lr))

        dummy_vec <- rep(NA, length(algout$merge_split_steps) + 1)
        # do unique original ancestors
        dummy_vec[!c(algout$merge_split_steps,FALSE)] <- apply(algout$original_ancestors_mat, 2, dplyr::n_distinct)
        nunique_original_ancestors <- dummy_vec
        # do unique parents
        dummy_vec[!c(algout$merge_split_steps,FALSE)] <- nunique_parent_indices
        nunique_parent_indices <- dummy_vec
        if(DEBUG_MODE) print("Checkpoint 5.1 - got summary info!")
        # make sizes null if needed
        if(!algout$plan_sizes_saved){
            algout$region_sizes_mat <- NULL
        }
        if(DEBUG_MODE) print("Checkpoint 5.5 - got summary info!")


        if (!is.nan(n_eff) && n_eff/nsims <= 0.05)
            cli_warn(c("Less than 5% resampling efficiency.",
                       "*" = "Increase the number of samples.",
                       "*" = "Consider weakening or removing constraints.",
                       "i" = "If sampling efficiency drops precipitously in the final
                                iterations, population balance is likely causing a bottleneck.
                                Try increasing {.arg pop_temper} by 0.01.",
                       "i" = "If sampling efficiency declines steadily across iterations,
                                adjusting {.arg seq_alpha} upward may help a bit."))

        # add the numerically stable weights back
        algout$wgt <- wgt

        # flatten the region sizes by column
        dim(algout$region_sizes_mat) <- NULL

        storage.mode(algout$original_ancestors_mat) <- "integer"
        storage.mode(algout$parent_index) <- "integer"

        if(DEBUG_MODE) print("Checkpoint 6 - before various diagnostics!")
        # Internal diagnostics,
        algout$internal_diagnostics <- list(
            parent_index_mat = algout$parent_index,
            original_ancestors_mat = algout$original_ancestors_mat,
            log_incremental_weights_mat = algout$log_incremental_weights_mat,
            draw_tries_mat = algout$draw_tries_mat,
            tree_sizes = algout$tree_sizes,
            successful_tree_sizes = algout$successful_tree_sizes,
            parent_unsuccessful_tries_mat = algout$parent_unsuccessful_tries_mat,
            region_ids_mat_list = algout$region_ids_mat_list,
            region_sizes_mat_list = algout$region_sizes_mat_list,
            merge_split_success_mat = algout$merge_split_success_mat,
            merge_split_attempt_counts = algout$merge_split_attempt_counts,
            forest_adjs_list = algout$forest_adjs_list,
            linking_edges_list = algout$linking_edges_list
        )

        # Information about the run
        algout$run_information <- list(
            weight_type=weight_type,
            nproc = nproc,
            num_threads = num_threads_per_process,
            custom_size_split_list=custom_size_split_list,
            valid_region_sizes_to_split_list=algout$valid_region_sizes_to_split_list,
            valid_split_region_sizes_list=algout$valid_split_region_sizes_list,
            sampling_space=sampling_space,
            split_method = split_method,
            splitting_size_regime = splitting_size_regime,
            merge_split_step_vec = merge_split_step_vec,
            ms_moves_multiplier = ms_moves_multiplier,
            merge_prob_type = merge_prob_type,
            step_types = step_types,
            nsims = nsims
        )

        # add high level diagnostic stuff
        # DOUBLE CHECK ALL THE SAME
        # Need to standardize with merge split
        algout$l_diag <- list(
            n_eff = n_eff,
            step_n_eff = algout$step_n_eff,
            adapt_k_thresh = split_params$adapt_k_thresh, # adapt_k_thresh, NEED TO DEAL WITH
            est_k = algout$est_k,
            split_params=split_params,
            accept_rate = algout$acceptance_rates,
            sd_lp = sd_lp,
            sd_temper = rep(NA, total_steps),
            unique_survive = nunique_parent_indices,
            ancestors = algout$ancestors,
            seq_alpha = seq_alpha,
            pop_temper = pop_temper,
            runtime = as.numeric(t2_run - t1_run, units = "secs"),
            nunique_original_ancestors = nunique_original_ancestors
        )

        algout

    }
    t2 <- Sys.time()

    if (verbosity >= 2) {
        cli_text("{format(nsims*runs, big.mark=',')} plans sampled in
                 {format(t2-t1, digits=2)}")
    }

    if(DEBUG_MODE) print("Checkpoint 7 - Out of for loop!")

    # combine if needed
    if(runs > 1){
        plans <- do.call(cbind, lapply(all_out, function(x) x$plans))
        region_sizes <- do.call(c, lapply(all_out, function(x) x$region_sizes_mat))
        wgt <- do.call(c, lapply(all_out, function(x) x$wgt))
        l_diag <- lapply(all_out, function(x) x$l_diag)
        run_information <- lapply(all_out, function(x) x$run_information)
        internal_diagnostics <- lapply(all_out, function(x) x$internal_diagnostics)
    }else{
        # else if just one run extract directly
        plans <- all_out[[1]]$plans
        region_sizes <- all_out[[1]]$region_sizes_mat
        wgt <- all_out[[1]]$wgt
        l_diag <- list(all_out[[1]]$l_diag)
        run_information <- list(all_out[[1]]$run_information)
        internal_diagnostics <- list(all_out[[1]]$internal_diagnostics)
    }


    if(DEBUG_MODE) print("Checkpoint 7.2 - Past the All Combine")


    n_dist_act <- dplyr::n_distinct(plans[, 1]) # actual number (for partial plans)

    alg_type <- ifelse(any_ms_steps_ran, "smc_ms","smc")
    if(DEBUG_MODE) print("Checkpoint 7.5 -About to create new plans!")
    out <- new_redist_plans(plans, map, alg_type, wgt, resample,
                            ndists = n_dist_act,
                            region_sizes = region_sizes,
                            n_eff = all_out[[1]]$n_eff,
                            compactness = compactness,
                            constraints = constraints,
                            version = packageVersion("redist"),
                            diagnostics = l_diag,
                            run_information = run_information,
                            internal_diagnostics = internal_diagnostics,
                            num_admin_units = num_admin_units,
                            entire_runtime = t2-t1)

    if(DEBUG_MODE) print("Checkpoint 8 - Created new plans!")
    if (runs > 1) {
        out <- mutate(out, chain = rep(seq_len(runs), each = n_dist_act*nsims)) %>%
            dplyr::relocate('chain', .after = "draw")
    }

    exist_name <- attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name) && ndists == final_dists) {
        if(DEBUG_MODE) print("Checkpoint 8.1 - adding reference plan!")
        ref_name <- if (!is.null(ref_name)) ref_name else exist_name
        out <- add_reference(out, map[[exist_name]], ref_name)
    }

    out

}



#' Helper function to truncate importance weights
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
