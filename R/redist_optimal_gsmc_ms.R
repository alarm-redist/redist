#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2024/08/18
# Purpose: tidy R wrapper to run gSMC with merge split steps
# redistricting code
####################################################


#' OptimalgSMC with merge split Redistricting Sampler
#'
#' @param k_params Either a single value to use as the splitting parameter for
#' every round or a vector of length N-1 where each value is the one to use for
#' a split.
#' @param multiprocess Whether or not to launch multiple processes (sometimes
#' better to disable to avoid using too much memory.)
#'
#' @return `redist_smc` returns a [redist_plans] object containing the simulated
#' plans.
#'
#' @export
redist_optimal_gsmc_ms <- function(state_map, M, counties = NULL,
                       k_params = 6, split_district_only = FALSE,
                       ms_freq = 1, ms_steps_multiplier = 1L,
                       run_ms = 0 < ms_freq && ms_freq <= N,
                       resample = TRUE, runs = 1L,
                       ncores = 0L, multiprocess=TRUE,
                       pop_temper = 0,
                       init_region_ids_mat = NULL, init_dvals_mat = NULL,
                       verbose = FALSE, silent = FALSE, diagnostic_mode = FALSE){
    N <- attr(state_map, "ndists")
    ndists <- attr(state_map, "ndists")
    nsims <- M

    # figure out the alg type
    if(split_district_only){
        alg_type <- "smc_ms"
    }else{
        alg_type <- "gsmc_ms"
    }

    # no constraints
    constraints <- list()


    # make controls intput
    lags <- 1 + unique(round((N - 1)^0.8*seq(0, 0.7, length.out = 4)^0.9))

    # check k param imput
    if(any(k_params < 1)) {
        cli_abort("K parameter values must be all at least 1.")
    }

    # if just a single number then repeat it
    if(length(k_params) == 1 && floor(k_params) == k_params){
        k_params <- rep(k_params, N-1)
    }else if(length(k_params) != N-1){
        cli_abort("K parameter input must be either 1 value or N-1!")
    }else if(any(floor(k_params) != k_params)){
        # if either the length is not N-1 or its not all integers then throw
        # error
        cli_abort("K parameter values must be all integers")
    }

    # create merge split parameter information

    # check that ms_steps_multiplier is an integer
    assertthat::assert_that(
        floor(ms_steps_multiplier) == ms_steps_multiplier &&
            ms_steps_multiplier > 0,
        msg = "`ms_steps_multiplier` must be a positive integer!"
    )

    # there are N-1 splits so for now just do it
    total_smc_steps <- N-1

    # check if there will be any merge split steps
    any_ms_steps_ran <- ms_freq <= total_smc_steps && run_ms

    merge_split_step_vec <- rep(FALSE, N-1)

    if(any_ms_steps_ran){
        # Now add merge split every `ms_freq` steps

        # insertion trick
        # https://stackoverflow.com/questions/1493969/insert-elements-into-a-vector-at-given-indexes
        ind <- seq(from = ms_freq, to = N-1, by = ms_freq)
        val <- c( merge_split_step_vec, rep(TRUE,length(ind)) )
        id  <- c( seq_along(merge_split_step_vec), ind+0.5 )

        # number of merge split is sum of trues
        merge_split_step_vec <- val[order(id)]
    }

    assertthat::assert_that(sum(!merge_split_step_vec) == total_smc_steps)
    # assert first step is not smc
    assertthat::assert_that(!merge_split_step_vec[1])

    total_ms_steps <- sum(merge_split_step_vec)

    # total number of steps to run
    total_steps <- total_smc_steps + total_ms_steps


    control <- list(
        lags=lags,
        pop_temper = pop_temper,
        k_params = k_params,
        split_district_only = split_district_only,
        merge_split_step_vec = merge_split_step_vec,
        ms_steps_multiplier = ms_steps_multiplier
        )

    # TODO fix this later
    est_k_params <- k_params
    est_k_params <- rep(-1, total_steps)
    est_k_params[!merge_split_step_vec] <- k_params

    for (i in 1:total_steps) {
        if(est_k_params[i] <= 0){
            est_k_params[i] <- est_k_params[i-1]
        }
    }

    assertthat::assert_that(
        all(est_k_params > 0),
        msg = "Something went wrong with est_k_params, fix it!"
    )

    # verbosity stuff
    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0


    # get the map in adjacency form
    map <- validate_redist_map(state_map)
    V <- nrow(state_map)
    adj_list <- get_adj(state_map)




    counties <- rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(counties)) {
        counties <- rep(1, V)
    } else {
        if (any(is.na(counties)))
            cli_abort("County vector must not contain missing values.")

        # handle discontinuous counties
        component <- contiguity(adj_list, vctrs::vec_group_id(counties))
        counties <- dplyr::if_else(component > 1,
                                   paste0(as.character(counties), "-", component),
                                   as.character(counties)) %>%
            as.factor() %>%
            as.integer()
        if (any(component > 1)) {
            cli_warn("Counties were not contiguous; expect additional splits.")
        }
    }

    ndists <- N
    n_steps <- ndists - 1

    # handle particle inits
    if (is.null(init_region_ids_mat)) {
        init_region_ids_mat <- matrix(0L, nrow = V, ncol = nsims)
        n_drawn <- 0L
    } else {
        if (nrow(init_region_ids_mat) != V)
            cli_abort("{.arg init_region_ids_mat} must have as many rows as {.arg map} has precincts.")
        if (ncol(init_region_ids_mat) != nsims)
            cli_abort("{.arg init_region_ids_mat} must have {.arg nsims} columns.")
        n_drawn <- as.integer(max(init_region_ids_mat[, 1]))
    }
    if (is.null(n_steps)) {
        n_steps <- ndists - n_drawn - 1L
    }
    final_dists <- n_drawn + n_steps + 1L
    if (final_dists > ndists) {
        cli_abort("Too many districts already drawn to take {n_steps} steps.")
    }

    # handle region dval inits
    if (is.null(init_dvals_mat)) {
        # initialize it so the first element of every column is ndists
        init_dvals_mat <- matrix(0L, nrow = ndists, ncol = nsims)
        init_dvals_mat[1,] <- ndists
    } else {
        if (nrow(init_dvals_mat) != V)
            cli_abort("{.arg init_dvals_mat} must have as many rows as {.arg map} has precincts.")
        if (ncol(init_dvals_mat) != nsims)
            cli_abort("{.arg init_dvals_mat} must have {.arg nsims} columns.")
    }
    if (!all(colSums(init_dvals_mat) == ndists)) {
        cli_abort("Too many dvals in already drawn plans.")
    }

    # get population stuff
    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]
    if (any(pop >= pop_bounds[3])) {
        too_big <- as.character(which(pop >= pop_bounds[3]))
        cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                population larger than the district target.",
                    "x" = "Redistricting impossible."))
    }

    # compute lags thing
    lags <- 1 + unique(round((N - 1)^0.8*seq(0, 0.7, length.out = 4)^0.9))


    # set up parallel
    ncores_max <- parallel::detectCores()
    ncores_runs <- min(ncores_max, runs)
    ncores_per <- as.integer(ncores)
    if (ncores_per == 0) {
        if (M/100*length(adj_list)/200 < 20) {
            ncores_per <- 1L
        } else {
            ncores_per <- floor(ncores_max/ncores_runs)
        }
    }

    # if sequentially
    if(!multiprocess){
        # either max cores if
        if(ncores == 0){
            ncores_per = ncores_max
        }else{
            ncores_per = ncores
        }
    }


    if (ncores_runs > 1 && multiprocess) {
        `%oper%` <- `%dorng%`
        of <- if (Sys.info()[["sysname"]] == "Windows") {
            tempfile(pattern = paste0("smc_", substr(Sys.time(), 1, 10)), fileext = ".txt")
        } else {
            ""
        }

        if (!silent)
            cl <- makeCluster(ncores_runs, outfile = of, methods = FALSE,
                              useXDR = .Platform$endian != "little")
        else
            cl <- makeCluster(ncores_runs, methods = FALSE,
                              useXDR = .Platform$endian != "little")
        doParallel::registerDoParallel(cl, cores = ncores_runs)
        on.exit(stopCluster(cl))
    } else {
        `%oper%` <- `%do%`
    }


    control[["num_threads"]] <- as.integer(ncores_per)

    t1 <- Sys.time()
    all_out <- foreach(chain = seq_len(runs), .inorder = FALSE, .packages="gredist") %oper% {


        run_verbosity <- if (chain == 1 || !multiprocess) verbosity else 0
        t1_run <- Sys.time()


        algout <- gredist::optimal_gsmc_with_merge_split_plans(
            N=N,
            adj_list=adj_list,
            counties=counties,
            pop=pop,
            target=pop_bounds[2],
            lower=pop_bounds[1],
            upper=pop_bounds[3],
            M=M,
            region_id_mat = init_region_ids_mat,
            region_dvals_mat = init_dvals_mat,
            control = control,
            verbosity=run_verbosity,
            diagnostic_mode=diagnostic_mode)



        if (length(algout) == 0) {
            cli::cli_process_done()
            cli::cli_process_done()
        }


        # add 1 to make it plans 1 indexed
        algout$plans_mat <- algout$plans_mat + 1L
        # add 1 to make parent mat  1-indexed for R indexing
        algout$parent_index <- algout$parent_index + 1L


        # convert the order added results into an actual list of arrays where
        # for each list entry n and column entry i that is a vector of length n
        # mapping the region id to its sorted order. The way to interpret is
        # the index where algout$region_order_added_list[[test_n]][,test_i] == r
        # is the new index region id r was mapped to.
        # In other words which(algout$region_order_added_list[[n]][,i] == r) is
        # the new ordered value r should be set to


        if(!diagnostic_mode){
            # if not diagnostic mode
            # make the region_ids_mat_list input just null since there's nothing else
            algout$region_ids_mat_list <- NULL
            algout$region_dvals_mat_list <- NULL
        }else{
            # make intermediate plans 1 indexed
            for (i in seq_len(length(algout$region_ids_mat_list))) {
                algout$region_ids_mat_list[[i]] <- algout$region_ids_mat_list[[i]] + 1L
            }
        }




        # if no merge split was run them remove those attributes
        if(!any_ms_steps_ran){
            algout$merge_split_success_mat <- NULL
            algout$merge_split_attempt_counts <- NULL
        }

        gc()


        # turn it into a character vector
        algout$step_split_types <- ifelse(
            algout$merge_split_steps, "ms", "smc"
        )

        num_ms_steps <- sum(
            algout$step_split_types == "ms"
        )



        # make parent succesful tries matrix counting the number of
        # times a parent index was successfully sampled
        parent_successful_tries_mat <- apply(
            algout$parent_index, 2, tabulate, nbins = M
        )


        # pull out the log weights
        lr <- algout$log_incremental_weights_mat[,total_smc_steps]

        wgt <- exp(lr - mean(lr))
        n_eff <- length(wgt)*mean(wgt)^2/mean(wgt^2)


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
        }


        t2_run <- Sys.time()

        # get original ancestor matrix from parent index
        algout$original_ancestors_mat <- get_original_ancestors_mat(
            algout$parent_index
        )

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

        if (!is.nan(n_eff) && n_eff/M <= 0.05)
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

        # add diagnostic stuff
        algout$l_diag <- list(
            n_eff = n_eff,
            step_n_eff = algout$step_n_eff,
            adapt_k_thresh = .99, # adapt_k_thresh, NEED TO DEAL WITH
            est_k = est_k_params, # algout$est_k,
            accept_rate = algout$acceptance_rates,
            sd_lp = sd_lp,
            sd_temper = rep(NA, total_steps), # algout$sd_temper,
            unique_survive = nunique_parent_indices,
            ancestors = algout$ancestors,
            seq_alpha = .99,
            pop_temper = pop_temper,
            runtime = as.numeric(t2_run - t1_run, units = "secs"),
            num_threads = ncores_per,
            nunique_original_ancestors = nunique_original_ancestors,
            parent_index_mat = algout$parent_index,
            original_ancestors_mat = algout$original_ancestors_mat,
            region_dvals_mat_list = algout$region_dvals_mat_list,
            log_incremental_weights_mat = algout$log_incremental_weights_mat,
            region_ids_mat_list = algout$region_ids_mat_list,
            draw_tries_mat = algout$draw_tries_mat,
            parent_unsuccessful_tries_mat = algout$parent_unsuccessful_tries_mat,
            parent_successful_tries_mat = parent_successful_tries_mat,
            rs_idx = rs_idx,
            merge_split_steps = algout$merge_split_steps,
            step_split_types = algout$step_split_types,
            merge_split_success_mat = algout$merge_split_success_mat,
            merge_split_attempt_counts = algout$merge_split_attempt_counts,
            num_ms_steps = num_ms_steps,
            ms_steps_multiplier = ms_steps_multiplier
        )


        algout

    }

    if (verbosity >= 2) {
        t2 <- Sys.time()
        cli_text("{format(M*runs, big.mark=',')} plans sampled in
                 {format(t2-t1, digits=2)}")
    }


    plans <- do.call(cbind, lapply(all_out, function(x) x$plans))
    wgt <- do.call(c, lapply(all_out, function(x) x$wgt))
    l_diag <- lapply(all_out, function(x) x$l_diag)
    n_dist_act <- dplyr::n_distinct(plans[, 1]) # actual number (for partial plans)




    out <- new_redist_plans(plans, map, alg_type, wgt, resample,
                            ndists = N,
                            n_eff = all_out[[1]]$n_eff,
                            compactness = 1,
                            constraints = constraints,
                            version = packageVersion("gredist"),
                            diagnostics = l_diag,
                            pop_bounds = pop_bounds)



    if (runs > 1) {
        out <- mutate(out, chain = rep(seq_len(runs), each = n_dist_act*M)) %>%
            dplyr::relocate('chain', .after = "draw")
    }



    exist_name <- attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name) && N == final_dists) {
        ref_name <- if (!is.null(ref_name)) ref_name else exist_name
        out <- add_reference(out, map[[exist_name]], ref_name)
    }

    out
}
