#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: parallel merge-split
####################################################

#' Parallel Merge-Split/Recombination MCMC Redistricting Sampler
#'
#' `redist_mergesplit_parallel()` is deprecated. Please call `redist_mergesplit`
#' now instead.
#'
#' @inherit redist_mergesplit details
#'
#' @inheritParams redist_mergesplit
#'
#' @returns A [`redist_plans`] object with all of the simulated plans, and an
#' additional `chain` column indicating the chain the plan was drawn from.
#'
#' @inherit redist_mergesplit references
#'
#' @examples \dontrun{
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#' sampled <- redist_mergesplit(fl_map, nsims = 100, chains = 100)
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_mergesplit_parallel <- function(
        map,
        nsims,
        warmup = if (is.null(init_plan)) 10 else max(100, nsims %/% 5),
        thin = 1L,
        chains = 1L,
        init_plan = NULL,
        counties = NULL,
        compactness = 1,
        constraints = list(),
        constraint_fn = function(m) rep(0, ncol(m)),
        sampling_space = c("graph_plan", "spanning_forest", "linking_edge"),
        split_method = NULL,
        split_params = NULL,
        merge_prob_type = "uniform",
        init_seats = NULL,
        ncores = NULL,
        cl_type = "PSOCK",
        return_all = TRUE,
        init_name = NULL,
        verbose = FALSE,
        silent = FALSE,
        diagnostic_mode = FALSE,
        control = list(),
        adapt_k_thresh = .99) {

    cli::cli_abort("redist_mergesplit_parallel is deprecated. Please call redist_mergesplit now.")
}

