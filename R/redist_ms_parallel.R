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
#' `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#' with its own sampled plan.
#' @param init_plan The initial state of the map, provided as a single vector
#' to be shared across all chains, or a matrix with `chains` columns.
#' If not provided, will default to the reference map of the map object, or if
#' none exists, will sample a random initial state using redist_smc. You can
#' also request a random initial state for each chain by setting
#' init_plan="sample".
#' @param ncores the number of parallel processes to run. Defaults to the
#' maximum available.
#' @param cl_type the cluster type (see [makeCluster()]). Safest is `"PSOCK"`,
#' but `"FORK"` may be appropriate in some settings.
#' @param return_all if `TRUE` return all sampled plans; otherwise, just return
#' the final plan from each chain.
#'
#' @returns A [`redist_plans`] object with all of the simulated plans, and an
#' additional `chain` column indicating the chain the plan was drawn from.
#'
#' @inherit redist_mergesplit references
#'
#' @examples \dontrun{
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#' sampled <- redist_mergesplit_parallel(fl_map, nsims = 100, chains = 100)
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_mergesplit_parallel <- function(map, nsims, chains = 1,
                                       warmup = if (is.null(init_plan)) 10 else max(100, nsims %/% 5),
                                       thin = 1L, init_plan = NULL, counties = NULL, compactness = 1,
                                       constraints = list(), constraint_fn = function(m) rep(0, ncol(m)),
                                       adapt_k_thresh = 0.99, k = NULL, ncores = NULL,
                                       cl_type = "PSOCK", return_all = TRUE, init_name = NULL,
                                       silly_adj_fix = FALSE,
                                       verbose = FALSE, silent = FALSE) {
    if (!missing(constraint_fn)) cli::cli_warn("{.arg constraint_fn} is deprecated.")

    cli::cli_abort("redist_mergesplit_parallel is deprecated. Please call redist_mergesplit now.")
}

