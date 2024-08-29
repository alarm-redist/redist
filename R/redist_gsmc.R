#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2024/08/18
# Purpose: tidy R wrapper to run gSMC redistricting code
####################################################


#' gSMC Redistricting Sampler
#'
#'
#' @return `redist_smc` returns a [redist_plans] object containing the simulated
#' plans.
#'
#' @export
redist_gsmc <- function(state_map, M, counties = NULL, k_param_val = 6,
                       resample = TRUE, runs = 1L, ncores = 0L,
                       verbose = FALSE, silent = FALSE){
    N <- attr(state_map, "ndists")
    # no constraints
    constraints = list()


    # make controls intput
    lags <- 1 + unique(round((N - 1)^0.8*seq(0, 0.7, length.out = 4)^0.9))
    control = list(lags=lags)

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


    if (ncores_runs > 1) {
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



    t1 <- Sys.time()
    all_out <- foreach(chain = seq_len(runs), .inorder = FALSE, .packages="redist") %oper% {


        run_verbosity <- if (chain == 1) verbosity else 0
        t1_run <- Sys.time()

        algout <- redist::generalized_smc_plans(
            N=N,
            adj_list=adj_list,
            counties=counties,
            pop=pop,
            target=pop_bounds[2],
            lower=pop_bounds[1],
            upper=pop_bounds[3],
            M=M,
            k_param = k_param_val,
            control = control,
            ncores = as.integer(ncores_per),
            verbosity=run_verbosity)


        if (length(algout) == 0) {
            cli::cli_process_done()
            cli::cli_process_done()
        }

        # add a plans matrix
        # IMPORTANT: make it 1 indexed
        algout$plans <- matrix(
            unlist(algout$region_ids_mat[[N-1]]),
            ncol = length(algout$region_ids_mat[[N-1]]),
            byrow = FALSE) + 1

        # make parent mat into matrix
        algout$parent_index <- matrix(
            unlist(algout$parent_index),
            ncol = length(algout$parent_index),
            byrow = FALSE)

        # pull out the log weights
        # IDK why its negative
        lr <- algout$log_incremental_weights_mat[[N-1]]

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
                                 {.url https://github.com/alarm-redist/redist/issues/new}"))
        }



        if (resample) {
            normalized_wgts <- wgt/sum(wgt)
            n_eff <- 1/sum(normalized_wgts^2)

            rs_idx <- resample_lowvar(normalized_wgts)
            n_unique <- dplyr::n_distinct(rs_idx)
            algout$plans <- algout$plans[, rs_idx, drop = FALSE]

            algout$ancestors <- algout$ancestors[rs_idx, , drop = FALSE]
            storage.mode(algout$ancestors) <- "integer"
        }

        storage.mode(algout$plans) <- "integer"
        t2_run <- Sys.time()


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
            est_k = rep(k_param_val, N-1), # algout$est_k,
            accept_rate = algout$acceptance_rates,
            sd_lp = c(
                sapply(algout$log_incremental_weights_mat, sd)  , sd(lr)
                ),
            sd_temper = rep(NA, N-1), # algout$sd_temper,
            unique_survive = c(algout$nunique_parent_indices, n_unique),
            ancestors = algout$ancestors,
            seq_alpha = .99,
            pop_temper = 0,
            runtime = as.numeric(t2_run - t1_run, units = "secs"),
            district_str_labels = algout$final_region_labs,
            nunique_original_ancestors = algout$nunique_original_ancestors,
            parent_index_mat = algout$parent_index,
            region_dvals_mat = algout$region_dvals_mat,
            log_incremental_weights_mat = algout$log_incremental_weights_mat,
            region_ids_mat = algout$region_ids_mat,
            draw_tries_mat = algout$draw_tries_mat
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


    out <- new_redist_plans(plans, map, "gsmc", wgt, resample,
                            ndists = N,
                            n_eff = all_out[[1]]$n_eff,
                            compactness = 1,
                            constraints = constraints,
                            version = packageVersion("redist"),
                            diagnostics = l_diag,
                            pop_bounds = pop_bounds)


    if (runs > 1) {
        out <- mutate(out, chain = rep(seq_len(runs), each = n_dist_act*M)) %>%
            dplyr::relocate('chain', .after = "draw")
    }

    exist_name <- attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name) && ndists == final_dists) {
        ref_name <- if (!is.null(ref_name)) ref_name else exist_name
        out <- add_reference(out, map[[exist_name]], ref_name)
    }

    out
}
