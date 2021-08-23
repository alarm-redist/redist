#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: parallel merge-split
####################################################

#' Parallel Merge-Split/Recombination MCMC Redistricting Sampler
#'
#' `redist_mergesplit_parallel()` runs [redist_mergesplit()] on several
#' chains in parallel.
#'
#' @inherit redist_mergesplit details
#'
#' @inheritParams redist_mergesplit
#' @param chains the number of parallel chains to run. Each chain will have
#'   `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#'   with its own sampled plan.
#' @param init_plan The initial state of the map, provided as a single vector
#'   to be shared across all chains, or a matrix with `chains` columns.
#'   If not provided, will default to the reference map of the map object, or if
#'   none exists, will sample a random initial state using redist_smc. You can
#'   also request a random initial state for each chain by setting
#'   init_plan="sample".
#' @param ncores the number of parallel processes to run. Defaults to the
#'   maximum available.
#' @param cl_type the cluster type (see [makeCluster()]). Safest is `"PSOCK"`,
#'  but `"FORK"` may be appropriate in some settings.
#' @param return_all if `TRUE` return all sampled plans; otherwise, just return
#'   the final plan from each chain.
#'
#' @returns A [`redist_plans`] object with all of the simulated plans, and an
#'   additional `chain` column indicating the chain the plan was drawn from.
#'
#' @inherit redist_mergesplit references
#'
#' @examples \dontrun{
#' data(fl25)
#' fl_map = redist_map(fl25, ndists=3, pop_tol=0.1)
#' sampled = redist_mergesplit_parallel(fl_map, nsims=100, chains=100)
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_mergesplit_parallel = function(map, nsims, chains=1, warmup=floor(nsims/2),
                                      init_plan=NULL, counties=NULL, compactness=1,
                                      constraints=list(), constraint_fn=function(m) rep(0, ncol(m)),
                                      adapt_k_thresh=0.975, k=NULL, ncores=NULL,
                                      cl_type="PSOCK", return_all=TRUE, init_name=NULL,
                                      verbose=TRUE, silent=FALSE) {
    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)
    ndists = attr(map, "ndists")

    chains = as.integer(chains)
    stopifnot(chains > 1)

    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")
    if (nsims < 1)
        stop("`nsims` must be positive.")

    exist_name = attr(map, "existing_col")
    counties = rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(init_plan) && !is.null(exist_name)) {
        init_plans = matrix(rep(as.integer(as.factor(get_existing(map))), chains), ncol=chains)
        if (is.null(init_name))
            init_names = rep(exist_name, chains)
        else
            init_names = rep(init_name, chains)
    } else if (!is.null(init_plan)) {
        if (is.matrix(init_plan)) {
            stopifnot(ncol(init_plan) == chains)
            init_plans = init_plan
        } else {
            init_plans = matrix(rep(as.integer(init_plan), chains), ncol=chains)
        }

        if (is.null(init_name))
            init_names = rep(exist_name, chains)
        else
            init_names = rep(init_name, chains)
    }
    if (length(init_plan) == 0L || isTRUE(init_plan == "sample")) {
        if (!silent) cat("Sampling initial plans with SMC")
        init_plans = get_plans_matrix(
            redist_smc(map, chains, counties, compactness, constraints,
                       TRUE, constraint_fn, adapt_k_thresh, ref_name=FALSE,
                       verbose=verbose, silent=silent))
        if (is.null(init_name))
            init_names = paste0("<init> ", seq_len(chains))
        else
            init_names = paste(init_name, seq_len(chains))
    }
    stopifnot(nrow(init_plans) == V)
    stopifnot(max(init_plans) == ndists)

    if (is.null(counties)) {
        counties = rep(1, V)
    } else {
        if (any(is.na(counties)))
            stop("County vector must not contain missing values.")

        # handle discontinuous counties
        component = contiguity(adj, as.integer(as.factor(counties)))
        counties = dplyr::if_else(component > 1,
                                  paste0(as.character(counties), "-", component),
                                  as.character(counties)) %>%
            as.factor() %>%
            as.integer()
    }

    # Other constraints
    constraints = eval_tidy(enquo(constraints), map)
    proc = process_smc_ms_constr(constraints, V)
    constraints = proc$constraints
    n_current = max(constraints$status_quo$current)

    verbosity = 1
    if (verbose) verbosity = 3
    if (silent) verbosity = 0
    if (is.null(k)) k = 0

    pop_bounds = attr(map, "pop_bounds")
    pop = map[[attr(map, "pop_col")]]
    init_pop = pop_tally(init_plans, pop, ndists)
    if (any(init_pop < pop_bounds[1]) | any(init_pop > pop_bounds[3]))
        stop("Provided initialization does not meet population bounds.")
    if (any(pop >= get_target(map)))
        stop("Units ", which(pop >= get_target(map)),
             " have population larger than the district target.\n",
             "Redistricting impossible.")

    # kind of hacky -- extract k=... from outupt
    if (!requireNamespace("utils", quietly=TRUE)) stop()
    out = utils::capture.output({
        x <- ms_plans(1, adj, init_plans[,1], counties, pop, ndists, pop_bounds[2],
                      pop_bounds[1], pop_bounds[3], compactness,
                      0, rep(1, ndists), ndists, 0, 0, 0, 1, rep(0, V),
                      0, 0, 0, rep(1, ndists), 0, beta_fractures = 0, adapt_k_thresh,
                      0L, verbosity=2)
    }, type="output")
    rm(x)
    k = as.integer(stats::na.omit(stringr::str_match(out, "Using k = (\\d+)")[,2]))

    # set up parallel
    if (is.null(ncores)) ncores = parallel::detectCores()
    ncores = min(ncores, chains)
    if (!silent)
        cl = makeCluster(ncores, setup_strategy="sequential", outfile="", methods=FALSE)
    else
        cl = makeCluster(ncores, setup_strategy="sequential", methods=FALSE)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    each_len = if (return_all) nsims - warmup else 1
    plans = foreach(chain=seq_len(chains), .combine=cbind) %dopar% {
        if (!silent) cat("Starting chain ", chain, "\n", sep="")
        algout = ms_plans(N = nsims+1L, l = adj, init = init_plans[, chain],
                          counties = counties, pop = pop, n_distr = ndists,
                          target = pop_bounds[2], lower = pop_bounds[1],
                          upper = pop_bounds[3], rho = compactness,
                          beta_sq = constraints$status_quo$strength,
                          current = constraints$status_quo$current,
                          n_current =  n_current,
                          beta_vra = constraints$vra$strength,
                          tgt_min = constraints$vra$tgt_vra_min,
                          tgt_other = constraints$vra$tgt_vra_other,
                          pow_vra = constraints$vra$pow_vra,
                          min_pop = proc$min_pop,
                          beta_vra_hinge = constraints$hinge$strength,
                          tgts_min = constraints$hinge$tgts_min,
                          beta_inc = constraints$incumbency$strength,
                          incumbents = constraints$incumbency$incumbents,
                          beta_splits = constraints$splits$strength,
                          beta_fractures = constraints$multisplits$strength,
                          thresh = adapt_k_thresh, k = k, verbosity=verbosity)
        if (return_all)
            algout$plans[, -1:-(warmup+1L), drop=FALSE]
        else
            algout$plans[, nsims+1L, drop=FALSE]
    }
    out = new_redist_plans(plans, map, "mergesplit", NULL, FALSE,
                           compactness = compactness,
                           constraints = constraints,
                           adapt_k_thresh = adapt_k_thresh) %>%
        mutate(chain = rep(seq_len(chains), each=each_len*ndists))

    if (!is.null(init_names) && !isFALSE(init_name)) {
        if (all(init_names[1] == init_names)) {
            out = add_reference(out, init_plans[, 1], init_names[1])
        } else {
            out = Reduce(function(cur, idx) {
                add_reference(cur, init_plans[, idx], init_names[idx]) %>%
                    mutate(chain = dplyr::coalesce(chain, idx))
            }, rev(seq_len(chains)), init=out)
        }
    }

    dplyr::relocate(out, chain, .after="draw")
}

utils::globalVariables("chain")
