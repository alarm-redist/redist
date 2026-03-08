#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: parallel merge-split (deprecated)
####################################################

#' @rdname redist_mergesplit
#' @export
redist_mergesplit_parallel <- function(map, nsims, chains = 1,
                                       warmup = if (is.null(init_plan)) 10 else max(100, nsims %/% 5),
                                       thin = 1L, init_plan = NULL, counties = NULL, compactness = 1,
                                       constraints = list(), constraint_fn = function(m) rep(0, ncol(m)),
                                       adapt_k_thresh = 0.99, k = NULL, ncores = NULL,
                                       cl_type = "PSOCK", return_all = TRUE, init_name = NULL,
                                       silly_adj_fix = FALSE, verbose = FALSE, silent = FALSE) {
    cli::cli_warn(c(
        "{.fn redist_mergesplit_parallel} is deprecated.",
        "i" = "Use {.fn redist_mergesplit} with the {.arg chains} argument instead."
    ), .frequency = "once", .frequency_id = "redist_mergesplit_parallel_deprecated")
    args <- list(map = map, nsims = nsims, warmup = warmup, thin = thin,
                 init_plan = init_plan, chains = chains, counties = counties,
                 compactness = compactness, constraints = constraints,
                 adapt_k_thresh = adapt_k_thresh, k = k, ncores = ncores,
                 cl_type = cl_type, return_all = return_all,
                 init_name = init_name, silly_adj_fix = silly_adj_fix,
                 verbose = verbose, silent = silent)
    if (!missing(constraint_fn)) args$constraint_fn <- constraint_fn
    do.call(redist_mergesplit, args)
}

utils::globalVariables("chain")
