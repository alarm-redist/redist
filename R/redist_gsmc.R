#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2024/08/18
# Purpose: tidy R wrapper to run gSMC with merge split steps
# redistricting code
####################################################


#' gSMC Redistricting Sampler (O'Sullivan, McCartan and Imai ???)
#'
#' `redist_gsmc` uses a Sequential Monte Carlo algorithm (O'Sullivan, McCartan and Imai ???)
#' to generate representative samples of congressional or legislative
#' redistricting plans according to contiguity, population, compactness, and
#' administrative boundary constraints.
#'
#' This function draws samples from a specific target measure controlled by
#' the `map` parameters.
#'
#' @inheritParams redist_smc
#' @param k_params Either a single value to use as the splitting parameter for
#' every round or a vector of length ndists-1 where each value is the one to use for
#' a split.
#' @param multiprocess Whether or not to launch multiple processes (sometimes
#' better to disable to avoid using too much memory. NOTE: Non multiprocessing
#' appears to introduce validation bugs right now)
#'
#' @return `redist_gsmc` returns a [redist_plans] object containing the simulated
#' plans.
#'
#' @export
redist_gsmc <- function(
        map, nsims, counties = NULL,
        runs = 1L,
        n_smc_steps = NULL,
        estimate_cut_k = TRUE,
       manual_k_params = 6, adapt_k_thresh = .9999,
       split_district_only = FALSE, weight_type = "optimal",
       ms_freq = 0,
       ms_steps_multiplier = 1L,
       run_ms = 0 < ms_freq && ms_freq <= ndists,
       merge_prob_type = "uniform",
       resample = TRUE,
       ncores = 0L, multiprocess=FALSE,
       pop_temper = 0,
       init_region_ids_mat = NULL,
       init_region_sizes_mat = NULL,
       min_region_cut_sizes = NULL, max_region_cut_sizes = NULL,
       permitted_split_region_sizes_list = NULL,
       permitted_presplit_region_sizes_list = NULL,
       verbose = FALSE, silent = FALSE, diagnostic_mode = FALSE){


    if(run_ms){
        cli_abort("Merge Split not supported at this moment!")
    }

    # no constraints
    constraints <- list()

    ndists <- attr(map, "ndists")

    # get the map parameters
    map_params <- get_map_parameters(map, counties)

    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds

    # handle particle inits
    if (is.null(init_region_ids_mat)) {
        init_region_ids_mat <- matrix(0L, nrow = V, ncol = nsims)
        init_num_regions <-  1
    } else { # if user input then check its valid
        init_num_regions <- length(unique(init_region_ids_mat[,1]))
        validate_initial_region_id_mat(init_region_ids_mat, V, nsims, init_num_regions)
    }
    if (is.null(n_smc_steps)) {
        n_smc_steps <- ndists - init_num_regions
    }
    final_dists <- init_num_regions + n_smc_steps
    if (final_dists > ndists) {
        cli_abort("Too many districts already drawn to take {n_smc_steps} steps.")
    }

    # handle region dval inits
    if (is.null(init_region_sizes_mat)) {
        # initialize it so the first element of every column is ndists
        init_region_sizes_mat <- matrix(0L, nrow = ndists, ncol = nsims)
        init_region_sizes_mat[1,] <- ndists
    } else {
        validate_initial_region_sizes_mat(init_region_sizes_mat, ndists, nsims, init_num_regions)
    }

    total_smc_steps <- n_smc_steps

    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")
    # check k param input if not estimating
    if(any(manual_k_params < 1)) {
        cli_abort("K parameter values must be all at least 1.")
    }
    # if just a single number then repeat it
    if(length(manual_k_params) == 1 && floor(manual_k_params) == manual_k_params){
        manual_k_params <- rep(manual_k_params, total_smc_steps)
    }else if(length(manual_k_params) != total_smc_steps){
        cli_abort("K parameter input must be either 1 value or number of smc steps!")
    }else if(any(floor(manual_k_params) != manual_k_params)){
        # if either the length is not ndists-1 or its not all integers then throw
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


    # figure out the alg type
    if(split_district_only && any_ms_steps_ran){
        alg_type <- "smc_ms"
    }else if(!split_district_only && any_ms_steps_ran){
        alg_type <- "gsmc_ms"
    }else if(split_district_only && !any_ms_steps_ran){
        alg_type <- "basic_smc"
    }else{
        alg_type <- "gsmc"
    }

    assertthat::assert_that(sum(!merge_split_step_vec) == total_smc_steps)
    # assert first step is not smc
    assertthat::assert_that(!merge_split_step_vec[1])

    total_ms_steps <- sum(merge_split_step_vec)

    # total number of steps to run
    total_steps <- total_smc_steps + total_ms_steps


    # setting the splitting size regime
    if(split_district_only){
        splitting_size_regime = "split_district_only"
    }else if(is.null(permitted_split_region_sizes_list)){
        # this is means allow any valid sizes
        splitting_size_regime = "any_valid_sizes"
    }else{ # else its custom
        if(is.null(min_region_cut_sizes)){
            min_region_cut_sizes <- rep(1, total_smc_steps)
        }
        if(is.null(max_region_cut_sizes)){
            # num regions-1,..., 1
            max_region_cut_sizes <-  seq(
                from = ndists - init_num_regions,
                by = -1,
                length.out = total_smc_steps
                )
        }
        # validate the cut sizes
        OLD_validate_cut_sizes(ndists, total_smc_steps, min_region_cut_sizes, max_region_cut_sizes)
        validate_cut_sizes (ndists, total_smc_steps, init_num_regions, permitted_split_region_sizes_list)

        splitting_size_regime = "custom"
    }


    # compute lags thing
    lags <- 1 + unique(round((ndists - 1)^0.8*seq(0, 0.7, length.out = 4)^0.9))


    control <- list(
        weight_type=weight_type,
        lags=lags,
        pop_temper = pop_temper,
        splitting_method = NAIVE_K_SPLITTING,
        splitting_size_regime = splitting_size_regime,
        k_params = manual_k_params,
        permitted_split_region_sizes_list=permitted_split_region_sizes_list,
        permitted_presplit_region_sizes_list=permitted_presplit_region_sizes_list,
        min_region_cut_sizes=min_region_cut_sizes,
        max_region_cut_sizes=max_region_cut_sizes,
        merge_split_step_vec = merge_split_step_vec,
        ms_steps_multiplier = ms_steps_multiplier,
        adapt_k_thresh = adapt_k_thresh,
        estimate_cut_k=estimate_cut_k,
        merge_prob_type = merge_prob_type
        )

    # verbosity stuff
    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0






    # set up parallel
    ncores_max <- parallel::detectCores()
    ncores_runs <- min(ncores_max, runs)
    ncores_per <- as.integer(ncores)
    if (ncores_per == 0) {
        if (nsims/100*length(adj_list)/200 < 20) {
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


        algout <- gredist::run_redist_gsmc(
            ndists=ndists,
            adj_list=adj_list,
            counties=counties,
            pop=pop,
            step_types=step_types,
            target=pop_bounds[2],
            lower=pop_bounds[1],
            upper=pop_bounds[3],
            region_id_mat = init_region_ids_mat,
            region_sizes_mat = init_region_sizes_mat,
            sampling_space = GRAPH_PLAN_SPACE_SAMPLING,
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
            algout$region_sizes_mat_list <- NULL
        }else{
            # make intermediate plans 1 indexed
            for (i in seq_len(length(algout$region_ids_mat_list))) {
                algout$region_ids_mat_list[[i]] <- algout$region_ids_mat_list[[i]] + 1L
            }
            # add plans as well

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
            algout$parent_index, 2, tabulate, nbins = nsims
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

        # add diagnostic stuff
        algout$l_diag <- list(
            estimate_cut_k=estimate_cut_k,
            n_eff = n_eff,
            step_n_eff = algout$step_n_eff,
            adapt_k_thresh = adapt_k_thresh, # adapt_k_thresh, NEED TO DEAL WITH
            est_k = algout$est_k,
            accept_rate = algout$acceptance_rates,
            sd_lp = sd_lp,
            sd_temper = rep(NA, total_steps),
            unique_survive = nunique_parent_indices,
            ancestors = algout$ancestors,
            seq_alpha = .99,
            pop_temper = pop_temper,
            runtime = as.numeric(t2_run - t1_run, units = "secs"),
            num_threads = ncores_per,
            permitted_split_region_sizes_list=permitted_split_region_sizes_list,
            permitted_presplit_region_sizes_list=permitted_presplit_region_sizes_list,
            nunique_original_ancestors = nunique_original_ancestors,
            parent_index_mat = algout$parent_index,
            original_ancestors_mat = algout$original_ancestors_mat,
            region_sizes_mat_list = algout$region_sizes_mat_list,
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
    t2 <- Sys.time()

    if (verbosity >= 2) {

        cli_text("{format(nsims*runs, big.mark=',')} plans sampled in
                 {format(t2-t1, digits=2)}")
    }


    plans <- do.call(cbind, lapply(all_out, function(x) x$plans))
    wgt <- do.call(c, lapply(all_out, function(x) x$wgt))
    l_diag <- lapply(all_out, function(x) x$l_diag)
    n_dist_act <- dplyr::n_distinct(plans[, 1]) # actual number (for partial plans)




    out <- new_redist_plans(plans, map, alg_type, wgt, resample,
                            ndists = ndists,
                            n_eff = all_out[[1]]$n_eff,
                            compactness = 1,
                            constraints = constraints,
                            version = packageVersion("gredist"),
                            diagnostics = l_diag,
                            pop_bounds = pop_bounds,
                            entire_runtime = t2-t1,
                            min_region_cut_sizes = min_region_cut_sizes,
                            max_region_cut_sizes = max_region_cut_sizes,
                            weight_type = weight_type,
                            merge_prob_type = merge_prob_type)



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

