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
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to generate. The chain will run for
#' `warmup+(nsims*thin)` steps.
#' @param chains the number of parallel chains to run. Each chain will have
#' `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#' with its own sampled plan.
#' @param warmup The number of warmup samples to discard. Recommended to be at
#' least the first 20% of samples, and in any case no less than around 100
#' samples, unless initializing from a random plan.
#' @param thin Save every `thin`-th sample after running warump. Defaults to no thinning (1).
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
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
#' the algorithm does not appear to be sampling from the target distribution.
#' Must be between 0 and 1.
#' @param k The number of edges to consider cutting after drawing a spanning
#' tree. Should be selected automatically in nearly all cases.
#' @param init_name a name for the initial plan, or \code{FALSE} to not include
#' the initial plan in the output.  Defaults to the column name of the
#' existing plan, or "\code{<init>}" if the initial plan is sampled.
#' @param verbose Whether to print out intermediate information while sampling.
#' Recommended.
#' @param silent Whether to suppress all diagnostic information.
#' @param ncores The number of clusters to spawn Defaults to the
#' maximum available detected by `parallel::detectCores()`.
#' @param cl_type the cluster type (see [makeCluster()]). Safest is `"PSOCK"`,
#' but `"FORK"` may be appropriate in some settings.
#' @param return_all if `TRUE` return all sampled plans; otherwise, just return
#' the final plan from each chain.
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
#' sampled_basic <- redist_mergesplit(fl_map, 10000, chains = 10)
#'
#' sampled_constr <- redist_mergesplit(fl_map, 10000, chains = 10,
#'     constraints = list(
#'     incumbency = list(strength = 1000, incumbents = c(3, 6, 25))
#'     )
#' )
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_mergesplit <- function(
    map, nsims, counties = NULL,
    warmup = if (is.null(init_plan)) 10 else max(100, nsims %/% 5),
    thin = 1L, chains = 1,
    init_plan = NULL,
    constraints = list(), constraint_fn = function(m) rep(0, ncol(m)),
    sampling_space = c("graph_plan", "spanning_forest", "linking_edge"),
    split_method = c("naive_top_k","uniform_valid_edge", "expo_bigger_abs_dev"),
    split_params = list(adapt_k_thresh = .99),
    merge_prob_type = "uniform", compactness = 1,
    ncores = NULL,
    cl_type = "PSOCK", return_all = TRUE, init_name = NULL,
    verbose = FALSE, silent = FALSE, diagnostic_mode = FALSE
) {
    if (!missing(constraint_fn)) cli_warn("{.arg constraint_fn} is deprecated.")


    # check default inputs
    sampling_space <- rlang::arg_match(sampling_space)
    split_method <- rlang::arg_match(split_method)


    # validate constraints
    constraints <- validate_constraints(map=map, constraints=rlang::enquo(constraints))



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


    thin <- as.integer(thin)

    chains <- as.integer(chains)

    if (compactness < 0)
        cli_abort("{.arg compactness} must be non-negative.")
    if (thin < 1)
        cli_abort("{.arg thin} must be a positive integer.")
    if (nsims < 1)
        cli_abort("{.arg nsims} must be positive.")

    #validate the splitting method and params
    split_params <- validate_sample_space_and_splitting_method(
        sampling_space, split_method, split_params, num_splitting_steps
    )


    exist_name <- attr(map, "existing_col")
    if (is.null(init_plan)) {
        if (!is.null(exist_name)) {
            init_plans <- matrix(rep(vctrs::vec_group_id(get_existing(map)), chains), ncol = chains)
            if (is.null(init_name)) {
                init_names <- rep(exist_name, chains)
            } else {
                init_names <- rep(init_name, chains)
            }
        } else {
            init_plan <- "sample"
        }
    } else if (!is.null(init_plan)) {
        if (is.matrix(init_plan)) {
            stopifnot(ncol(init_plan) == chains)
            init_plans <- init_plan
        } else {
            init_plans <- matrix(rep(as.integer(init_plan), chains), ncol = chains)
        }
        if (is.null(init_name))
            init_names <- paste0("<init> ", seq_len(chains))
        else
            init_names <- rep(init_name, chains)
    }
    if (isTRUE(init_plan == "sample")) {
        if (!silent) cat("Sampling initial plans with SMC\n")
        # heuristic. Do at least 50 plans to not get stuck
        n_smc_nsims <- max(chains, 50)

        init_plans <- get_plans_matrix(
            redist_smc(map, n_smc_nsims, counties, compactness, constraints,
                       resample = TRUE, split_params = split_params,
                       sampling_space = sampling_space, split_method = split_method,
                       ref_name = FALSE, verbose = verbose, silent = silent)
            )[, sample.int(n=n_smc_nsims, size=chains, replace=F), drop=FALSE]

        if (is.null(init_name))
            init_names <- paste0("<init> ", seq_len(chains))
        else
            init_names <- paste(init_name, seq_len(chains))
    }



    # subtract 1 to make it 0 indexed
    init_plans <- init_plans - 1
    # validate initial plans
    validate_initial_region_id_mat(init_plans, V, chains, ndists)

    # check it satsifies population bounds
    # add one because we assume 1 indexed
    init_pop <- pop_tally(init_plans+1, pop, ndists)
    if (any(init_pop < pop_bounds[1]) | any(init_pop > pop_bounds[3]))
        cli_abort("Provided initialization does not meet population bounds.")

    # TODO: add support for multimember district in future
    init_sizes <- matrix(1L, nrow = ndists, ncol = ncol(init_plans))



    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0


    control = list(
        splitting_method=split_method
    )

    # add the splitting parameters
    # need to do it like this because its a list
    control <- c(control, split_params)



    # set up parallel
    if (is.null(ncores)) ncores <- parallel::detectCores()
    ncores <- min(ncores, chains)

    multiprocess <- ncores > 1

    if (multiprocess) {
        `%oper%` <- `%dorng%`

        of <- ifelse(Sys.info()[['sysname']] == 'Windows',
                     tempfile(pattern = paste0('ms_', substr(Sys.time(), 1, 10)), fileext = '.txt'),
                     '')
        # this makes a cluster using socket (NOT FORK) with
        if (!silent)
            cl <- makeCluster(ncores, outfile = of, methods = FALSE,
                              useXDR = .Platform$endian != "little")
        else
            cl <- makeCluster(ncores, methods = FALSE,
                              useXDR = .Platform$endian != "little")

        # this makes it avoid printing the loading required package message each time
        parallel::clusterEvalQ(cl, {
            suppressPackageStartupMessages(library(foreach))
            suppressPackageStartupMessages(library(rngtools))
            suppressPackageStartupMessages(library(redist))
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
        doParallel::registerDoParallel(cl, cores = ncores)
        on.exit(stopCluster(cl))
    } else {
        `%oper%` <- `%do%`
    }


    t1 <- Sys.time()
    out_par <- foreach(chain = seq_len(chains), .inorder = FALSE, .packages="redist") %oper% {
        if(chain == 1){
            is_chain1 <- T
        }
        run_verbosity <- if (is_chain1 || !multiprocess) verbosity else 0
        if (!silent && is_chain1){
            cat("Starting chain ", chain, "\n", sep = "")
            # flush.console()
            }


        t1_run <- Sys.time()
        algout <- ms_plans(
            nsims=nsims, warmup=warmup, thin=thin,
            ndists=ndists, total_seats=total_seats,
            district_seat_sizes=district_seat_sizes,
            adj_list=adj_list, counties=counties, pop=pop,
            target=pop_bounds[2], lower=pop_bounds[1], upper=pop_bounds[3],
            rho=compactness,
            initial_plan=init_plans[, chain, drop=FALSE],
            initial_region_sizes=init_sizes[, chain, drop=FALSE],
            sampling_space_str = sampling_space,
            merge_prob_type = merge_prob_type,
            control=control, constraints=constraints,
            verbosity=run_verbosity, diagnostic_mode=diagnostic_mode
        )
        t2_run <- Sys.time()
        # 1 index the plans
        algout$plans <- algout$plans

        # Internal diagnostics,
        algout$internal_diagnostics <- list(
            log_mh_ratio = algout$log_mh_ratio,
            mh_acceptance = algout$mhdecisions,
            warmup_acceptances = algout$warmup_acceptances,
            post_warump_acceptances = algout$post_warump_acceptances,
            warmup_accept_rate = algout$warmup_acceptances / warmup,
            postwarmup_accept_rate = algout$post_warump_acceptances / (nsims*thin),
            tree_sizes = algout$tree_sizes,
            successful_tree_sizes = algout$successful_tree_sizes,
            proposed_plans = algout$proposed_plans
        )

        algout$l_diag <- list(
            est_k = algout$est_k,
            runtime = as.numeric(t2_run - t1_run, units = "secs"),
            nsims = nsims,
            thin = thin,
            warmup = warmup,
            total_acceptances = (algout$warmup_acceptances + algout$post_warump_acceptances),
            accept_rate = (algout$warmup_acceptances + algout$post_warump_acceptances) / algout$total_steps,
            total_steps = algout$total_steps,
            split_params=split_params
        )

        # Information about the run
        algout$run_information <- list(
            valid_region_sizes_to_split_list=algout$valid_region_sizes_to_split_list,
            valid_split_region_sizes_list=algout$valid_split_region_sizes_list,
            sampling_space=sampling_space,
            split_method = split_method,
            merge_prob_type = merge_prob_type,
            nsims = nsims,
            alg_name = "mergesplit"
        )

        storage.mode(algout$plans) <- "integer"

        algout
    }
    t2 <- Sys.time()


    plans <- lapply(out_par, function(algout) {
        algout$plans
    })
    each_len <- ncol(plans[[1]])
    plans <- do.call(cbind, plans)

    mh <- sapply(out_par, function(algout) {
        mean(as.logical(algout$mhdecisions))
    })
    l_diag <- lapply(out_par, function(algout) algout$l_diag)
    run_information <- lapply(out_par, function(x) x$run_information)
    internal_diagnostics <- lapply(out_par, function(x) x$internal_diagnostics)

    acceptances <- sapply(out_par, function(algout) {
        algout$mhdecisions
    })


    out <- new_redist_plans(plans = plans, map = map, algorithm = "mergesplit",
                            wgt = NULL, resampled = FALSE,
                            compactness = compactness,
                            constraints = constraints,
                            ndists = ndists,
                            mh_acceptance = mh,
                            version = packageVersion("redist"),
                            diagnostics = l_diag,
                            run_information = run_information,
                            internal_diagnostics = internal_diagnostics,
                            pop_bounds = pop_bounds,
                            entire_runtime = t2-t1) %>%
        mutate(chain = rep(seq_len(chains), each = each_len*ndists),
               mcmc_accept = rep(acceptances, each = ndists))

    if (!is.null(init_names) && !isFALSE(init_name)) {
        if (all(init_names[1] == init_names)) {
            out <- add_reference(out, init_plans[, 1], init_names[1])
        } else {
            out <- Reduce(function(cur, idx) {
                add_reference(cur, init_plans[, idx], init_names[idx]) %>%
                    mutate(chain = dplyr::coalesce(chain, idx))
            }, rev(seq_len(chains)), init = out)
        }
    }

    dplyr::relocate(out, chain, .after = "draw")
}

utils::globalVariables("chain")
