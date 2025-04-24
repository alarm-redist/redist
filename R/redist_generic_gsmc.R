##############################################
## Author: Philip O'Sullivan
## Institution: Harvard University
## Date Created: 2025/02/15
## Purpose: Wrapper for running gsmc cpp code
##############################################



#' gSMC Redistricting Sampler (O'Sullivan, McCartan and Imai ???)
#'
#' `generic_redist_gsmc` uses a Sequential Monte Carlo algorithm (O'Sullivan, McCartan and Imai ???)
#' to generate representative samples of congressional or legislative
#' redistricting plans according to contiguity, population, compactness, and
#' administrative boundary constraints.
#'
#' Sampling can be performed either on the space of plans (graph partitions) or
#' on the space of spanning forests. This is an internal function designed to
#' support either sampling space and any arbitrary tree splitting method. The
#' functions `redist_gsmc` and `treedist_gsmc` are largely wrappers to calling
#' this function.
#'
#' This function draws samples from a specific target measure controlled by
#' the `map` parameters.
#'
#' @inheritParams redist_smc
#' @param split_district_only Whether or not to split plans by splitting off one
#' district at a time.
#' @param sampling_space The space to sample the plans on. Right now the supported
#' options are graph partitions (`r GRAPH_PLAN_SPACE_SAMPLING`) and spanning
#' forests (`r FOREST_SPACE_SAMPLING`)
#' @param splitting_method The method used for splitting spanning trees in the
#' sampling process. When sampling on the space of graph partitions it must be
#' the naive top k method but any method is allowed for forest space sampling.
#' @param splitting_params A list of parameters associated with the splitting
#' method passed in. Specific parameters depends on the splitting types
#' @param weight_type The type of SMC weights to use. Optimal weights typically
#' have lower variance but can be more computationally expensive, especially for
#' complex constraints or when sampling on the space of spanning forests.
#' @param multiprocess Whether or not to launch multiple processes (sometimes
#' better to disable to avoid using too much memory. NOTE: Non multiprocessing
#' appears to introduce validation bugs right now)
#' @param num_processes The number of processes (independent instances of R to
#' spawn) spawned to simulate the plans. The processes can execute runs in parallel
#' @param num_threads_per_process The number of threads assigned to each process.
#' This is the number of threads used when performing a specific run. If simulations
#' are memory constrained it can be better to lower the number of processes and increase
#' the threads per process.
#' @param init_region_ids_mat A matrix of partial plans to begin sampling from. For
#' advanced use only.  The matrix must have `nsims` columns and a row for
#' every precinct. It is important to ensure that the existing districts meet
#' contiguity and population constraints, or there may be major issues when
#' sampling.
#' @param init_region_sizes_mat A matrix of region sizes of the partial plans to
#' begin sampling from. For advanced use only.  The matrix must have `nsims`
#' columns and at least as many rows as the number of regions in the
#' `init_region_ids_mat` and each column must sum to the number of seats in the
#' map.
#' @param num_splitting_steps How many steps to run the SMC algorithm for.
#' Each step splits off a new district. Defaults to all remaining districts.
#' If fewer than the number of remaining splits, reference plans are disabled.
#'
#' @return `redist_gsmc` returns a [redist_plans] object containing the simulated
#' plans.
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
generic_redist_gsmc <- function(
        map, nsims, counties = NULL, compactness = 1,
        constraints = list(),
        runs = 1L,
        split_district_only = FALSE, weight_type = "optimal",
        sampling_space, splitting_method, splitting_params,
        ms_freq = 0,
        ms_steps_multiplier = 1L,
        run_ms = 0 < ms_freq && ms_freq <= attr(map, "ndists"),
        merge_prob_type = "uniform",
        resample = TRUE,
        num_processes=0L, num_threads_per_process=0L,
        multiprocess=TRUE,
        pop_temper = 0,
        init_region_ids_mat = NULL,
        init_region_sizes_mat = NULL,
        custom_size_split_list = NULL,
        num_splitting_steps = NULL,
        ref_name = NULL,
        verbose = FALSE, silent = FALSE, diagnostic_level = 0,
        counties_q = NULL, use_counties_q = F)
{

    if (compactness < 0)
        cli_abort("{.arg compactness} must be non-negative.")

    # if graph space default to k stuff
    if(sampling_space == GRAPH_PLAN_SPACE_SAMPLING){
        if(missing(splitting_method)){
            splitting_method <- NAIVE_K_SPLITTING
        }
        if(missing(splitting_params)){
            splitting_params = list(
                adapt_k_thresh=.99
            )
        }
    }else if(sampling_space == FOREST_SPACE_SAMPLING || sampling_space == LINKING_EDGE_SPACE_SAMPLING){
        # the others default to uniform
        if(missing(splitting_method)){
            splitting_method <- UNIF_VALID_EDGE_SPLITTING
        }
    }

    # validate constraints
    constraints_q <- rlang::enquo(constraints)
    constraints <- validate_constraints(map=map, constraints_q=constraints_q, use_constraints_q=TRUE)
    # get the total number of districts
    ndists <- attr(map, "ndists")

    # update quosure if quosure not passed in
    if(!use_counties_q){
        counties_q <- rlang::enquo(counties)
    }

    map_params <- get_map_parameters(map, counties_q=counties_q, use_counties_q=T)
    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    num_admin_units <- length(unique(counties))
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds


    # handle particle inits
    if (is.null(init_region_ids_mat)) {
        # if no initial plans passed in then create empty matrix
        init_region_ids_mat <- matrix(0L)
        init_num_regions <- 1L
    } else {
        # make sure size matrix passed in as well, inferring size not supported now
        if(is.null(init_region_sizes_mat)){
            cli::cli_abort("Initial Plans Matrix provided but no associated sizes. Inferring region sizes is not supported right now!")
        }
        # if user input then check its valid
        init_num_regions <- length(unique(init_region_ids_mat[,1]))
        validate_initial_region_id_mat(init_region_ids_mat, V, nsims, init_num_regions)
    }
    if (is.null(num_splitting_steps)) {
        num_splitting_steps <- ndists - init_num_regions
    }
    final_dists <- init_num_regions + num_splitting_steps
    if (final_dists > ndists) {
        cli_abort("Too many districts already drawn to take {num_splitting_steps} steps.")
    }

    # handle region size inits
    if (is.null(init_region_sizes_mat)) {
        # if no initial plans passed in then create empty matrix
        init_region_sizes_mat <- matrix(0L)
    } else {
        validate_initial_region_sizes_mat(init_region_sizes_mat, ndists, nsims, init_num_regions)
    }

    #validate the splitting method and params
    splitting_params <- validate_sample_space_and_splitting_method(
        sampling_space, splitting_method, splitting_params, num_splitting_steps
    )
    total_smc_steps <- num_splitting_steps

    # check weights are ok
    if(!weight_type %in% c("optimal", "simple")){
        cli_abort("{.arg weight_type} must be either `optimal` or `simple`!")
    }

    # create merge split parameter information

    # check that ms_steps_multiplier is an integer
    assertthat::assert_that(
        floor(ms_steps_multiplier) == ms_steps_multiplier &&
            ms_steps_multiplier > 0,
        msg = "`ms_steps_multiplier` must be a positive integer!"
    )


    # check if there will be any merge split steps
    any_ms_steps_ran <- ms_freq <= total_smc_steps && run_ms

    merge_split_step_vec <- rep(FALSE, total_smc_steps)

    if(any_ms_steps_ran){
        # Now add merge split every `ms_freq` steps
        # insertion trick
        # https://stackoverflow.com/questions/1493969/insert-elements-into-a-vector-at-given-indexes
        ind <- seq(from = ms_freq, to = total_smc_steps, by = ms_freq)
        val <- c( merge_split_step_vec, rep(TRUE,length(ind)) )
        id  <- c( seq_along(merge_split_step_vec), ind+0.5 )

        # number of merge split is sum of trues
        merge_split_step_vec <- val[order(id)]
    }


    # get the types
    step_types <- ifelse(merge_split_step_vec, "ms", "smc")
    assertthat::assert_that(sum(!merge_split_step_vec) == total_smc_steps)
    # assert first step is not smc
    assertthat::assert_that(!merge_split_step_vec[1])
    total_ms_steps <- sum(merge_split_step_vec)
    # total number of steps to run
    total_steps <- total_smc_steps + total_ms_steps


    # setting the splitting size regime
    if(split_district_only){
        splitting_size_regime = "split_district_only"
    }else if(is.null(custom_size_split_list)){
        # this is means allow any valid sizes
        splitting_size_regime = "any_valid_sizes"
    }else{ # else its custom
        # only support doing a single size right now
        # validate it
        validate_custom_size_split_list(ndists, num_splitting_steps, init_num_regions, custom_size_split_list)
        splitting_size_regime = "one_custom_size"
    }

    # compute lags thing
    lags <- 1 + unique(round((ndists - 1)^0.8*seq(0, 0.7, length.out = 4)^0.9))


    # verbosity stuff
    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0

    # set up parallel processing stuff
    ncores_max <- parallel::detectCores()


    if(num_processes > ncores_max){
        cli_warn("Inputted number of processes to spawn is greater than detected number of cores on machine")
    }else if(num_processes == 0){
        if(multiprocess){# if multiprocess then spawn min(runs, max cores) processes
            num_processes <- min(runs, ncores_max)
        }else{
            num_processes <- 1
        }
    }else{
        # make sure we're not spawning more proccesses than runs
        num_processes <- min(runs, num_processes)
    }

    if (num_threads_per_process == 0) {
        if(!multiprocess || num_processes == 1){
            # if no multiprocessing then the single process gets all threads
            num_threads_per_process <- ncores_max
        }else{
            # else each process gets ncores_max/num_processes threads
            num_threads_per_process <- floor(ncores_max/num_processes)
            num_threads_per_process <- max(1, num_threads_per_process)
        }
    }
    num_threads_per_process <- as.integer(num_threads_per_process)


    if (num_processes > 1 && multiprocess && runs > 1) {
        `%oper%` <- `%dorng%`
        of <- if (Sys.info()[["sysname"]] == "Windows") {
            tempfile(pattern = paste0("smc_", substr(Sys.time(), 1, 10)), fileext = ".txt")
        } else {
            ""
        }

        # this makes a cluster using socket (NOT FORK) with
        if (!silent)
            cl <- makeCluster(num_processes, outfile = of, methods = FALSE,
                              useXDR = .Platform$endian != "little")
        else
            cl <- makeCluster(num_processes, methods = FALSE,
                              useXDR = .Platform$endian != "little")
        # this makes it avoid printing the loading required package message each time
        parallel::clusterEvalQ(cl, {
            suppressPackageStartupMessages(library(foreach))
            suppressPackageStartupMessages(library(rngtools))
            suppressPackageStartupMessages(library(gredist))
        })
        # weird code, probably remove in production and find better way to ensure printing
        # but essentially makes it so only one process will print but if more runs then processes
        # it doesn't just print once
        parallel::clusterEvalQ(cl, {
            if (!exists("is_chain1", envir = .GlobalEnv)) {
                is_chain1 <- FALSE
            }
            NULL
        })
        doParallel::registerDoParallel(cl, cores = num_processes)
        on.exit(stopCluster(cl))
        cat("Spawning " , num_processes, " clusters \n")
    } else {
        `%oper%` <- `%do%`
    }


    control <- list(
        weight_type=weight_type,
        lags=lags,
        pop_temper = pop_temper,
        num_threads=as.integer(num_threads_per_process),
        splitting_method = splitting_method,
        splitting_size_regime = splitting_size_regime,
        custom_size_split_list=custom_size_split_list,
        merge_split_step_vec = merge_split_step_vec,
        ms_steps_multiplier = ms_steps_multiplier,
        merge_prob_type = merge_prob_type
    )

    # add the splitting parameters
    control <- c(control, splitting_params)


    t1 <- Sys.time()
    # all_out <- foreach(chain = seq_len(runs), .inorder = FALSE, .packages="gredist") %oper% {
    all_out <- foreach(chain = seq_len(runs), .inorder = FALSE) %oper% {
        if(chain == 1){
            is_chain1 <- T
        }

        if(is_chain1 && !silent) cat("Starting Chain ", chain, "\n", sep = "")
        run_verbosity <- if (is_chain1 || !multiprocess) verbosity else 0
        t1_run <- Sys.time()

        algout <- gredist::run_redist_gsmc(
            nsims=nsims,
            ndists=ndists,
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
            control = control,
            constraints = constraints,
            verbosity=run_verbosity,
            diagnostic_level=diagnostic_level,
            region_id_mat = init_region_ids_mat,
            region_sizes_mat = init_region_sizes_mat
        )

        if (length(algout) == 0) {
            cli::cli_process_done()
        }
        print("Out!")



        # convert the order added results into an actual list of arrays where
        # for each list entry n and column entry i that is a vector of length n
        # mapping the region id to its sorted order. The way to interpret is
        # the index where algout$region_order_added_list[[test_n]][,test_i] == r
        # is the new index region id r was mapped to.
        # In other words which(algout$region_order_added_list[[n]][,i] == r) is
        # the new ordered value r should be set to

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
        print("Checkpoint 1!")

        # if no merge split was run them remove those attributes
        if(!any_ms_steps_ran){
            algout$merge_split_success_mat <- NULL
            algout$merge_split_attempt_counts <- NULL
        }

        gc()
        print("Checkpoint 2 - gc!")

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
        lr <- algout$log_incremental_weights_mat[,total_smc_steps]

        wgt <- exp(lr - mean(lr))
        n_eff <- length(wgt)*mean(wgt)^2/mean(wgt^2)

        print("Checkpoint 3 - weight and lr!")
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
                                 {.url https://github.com/alarm-gredist/gredist/issues/new}"))
        }



        if (resample) {
            normalized_wgts <- wgt/sum(wgt)
            n_eff <- 1/sum(normalized_wgts^2)

            rs_idx <- resample_lowvar(normalized_wgts)

            n_unique <- dplyr::n_distinct(rs_idx)
            # makes algout$plans[i] now equal to algout$plans[rs_idx[i]]
            algout$plans_mat <- algout$plans_mat[, rs_idx, drop = FALSE]
            # now adjust for the resampling
            algout$ancestors <- algout$ancestors[rs_idx, , drop = FALSE]

            # add a final column for the resampling
            # NOTE: I THINK THIS IS WRONG, MIGHT NEED TO FLIP COLUMN
            algout$parent_index <- cbind(algout$parent_index, rs_idx[1:length(rs_idx)])

            # do unique parents
            nunique_parent_indices <- c(
                algout$nunique_parent_indices,
                dplyr::n_distinct(rs_idx[1:length(rs_idx)]))


            if(diagnostic_mode){
                # makes algout$final_region_labs[i] now equal to algout$final_region_labs[rs_idx[i]]
                # to account for resampling
                # algout$final_region_labs[,rs_idx, drop = FALSE]
            }

            #TODO probably need to adjust the rest of these as well
            storage.mode(algout$ancestors) <- "integer"
        }else{
            nunique_parent_indices <- algout$nunique_parent_indices
        }
        print("Checkpoint 4 - after resample!")

        t2_run <- Sys.time()
        # get original ancestor matrix from parent index
        algout$original_ancestors_mat <- get_original_ancestors_mat(
            algout$parent_index
        )
        print("Checkpoint 5 - after og anvestor mat!")


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




        # nunique_original_ancestors <- c(nunique_original_ancestors,
        #                          dplyr::n_distinct(algout$original_ancestors_mat[, ncol(algout$original_ancestors_mat)]))

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

        storage.mode(algout$original_ancestors_mat) <- "integer"
        storage.mode(algout$parent_index) <- "integer"

        print("Checkpoint 6 - before various diagnostics!")
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
            num_processes = num_processes,
            num_threads = num_threads_per_process,
            custom_size_split_list=custom_size_split_list,
            valid_region_sizes_to_split_list=algout$valid_region_sizes_to_split_list,
            valid_split_region_sizes_list=algout$valid_split_region_sizes_list,
            sampling_space=sampling_space,
            splitting_method = splitting_method,
            splitting_size_regime = splitting_size_regime,
            merge_split_step_vec = merge_split_step_vec,
            ms_steps_multiplier = ms_steps_multiplier,
            merge_prob_type = merge_prob_type,
            step_types = step_types,
            nsims = nsims
        )

        # add high level diagnostic stuff
        algout$l_diag <- list(
            n_eff = n_eff,
            step_n_eff = algout$step_n_eff,
            adapt_k_thresh = splitting_params$adapt_k_thresh, # adapt_k_thresh, NEED TO DEAL WITH
            est_k = algout$est_k,
            splitting_params=splitting_params,
            accept_rate = algout$acceptance_rates,
            sd_lp = sd_lp,
            sd_temper = rep(NA, total_steps),
            unique_survive = nunique_parent_indices,
            ancestors = algout$ancestors,
            seq_alpha = .99,
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
    #return(all_out)
    print("Checkpoint 7 - Out of for loop!")


    plans <- do.call(cbind, lapply(all_out, function(x) x$plans))
    plan_sizes <- do.call(cbind, lapply(all_out, function(x) x$plan_sizes))
    wgt <- do.call(c, lapply(all_out, function(x) x$wgt))
    l_diag <- lapply(all_out, function(x) x$l_diag)
    run_information <- lapply(all_out, function(x) x$run_information)
    internal_diagnostics <- lapply(all_out, function(x) x$internal_diagnostics)
    n_dist_act <- dplyr::n_distinct(plans[, 1]) # actual number (for partial plans)

    alg_type <- ifelse(any_ms_steps_ran, "smc_ms","smc")
    out <- new_redist_plans(plans, map, alg_type, wgt, resample,
                            ndists = n_dist_act,
                            n_eff = all_out[[1]]$n_eff,
                            compactness = compactness,
                            constraints = constraints,
                            version = packageVersion("gredist"),
                            diagnostics = l_diag,
                            plan_sizes = plan_sizes,
                            run_information = run_information,
                            internal_diagnostics = internal_diagnostics,
                            pop_bounds = pop_bounds,
                            num_admin_units = num_admin_units,
                            entire_runtime = t2-t1)

    if (runs > 1) {
        out <- mutate(out, chain = rep(seq_len(runs), each = n_dist_act*nsims)) %>%
            dplyr::relocate('chain', .after = "draw")
    }

    exist_name <- attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name) && ndists == final_dists) {
        ref_name <- if (!is.null(ref_name)) ref_name else exist_name
        out <- add_reference(out, map[[exist_name]], ref_name)
    }

    out

}
