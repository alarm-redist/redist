#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: tidy R wrapper to run Merge-Split/Recom redistricting code
####################################################

#' Merge-Split/Recombination MCMC Redistricting Sampler
#'
#' \code{redist_mergesplit} uses a Markov Chain Monte Carlo algorithm to
#' generate congressional or legislative redistricting plans according to
#' contiguity, population, compactness, and administrative boundary constraints.
#' The MCMC proposal is the same as is used in the SMC sampler; it is similar
#' but not identical to those used in the references.
#'
#' This function draws samples from a specific target measure, controlled by the
#' \code{popcons}, \code{compactness}, \code{constraints}, and
#' \code{constraint_fn} parameters.
#'
#' Higher values of \code{compactness} sample more compact districts;
#' setting this parameter to 1 is computationally efficient and generates nicely
#' compact districts.
#'
#' The \code{constraints} parameter allows the user to apply several common
#' redistricting constraints without implementing them by hand. This parameter
#' is a list, which may contain any of the following named entries:
#' * \code{status_quo}: a list with two entries:
#'   * \code{strength}, a number controlling the tendency of the generated districts
#'   to respect the status quo, with higher values preferring more similar
#'   districts.
#'   * \code{current}, a vector containing district assignments for
#'   the current map.
#' * \code{vra}: a list with five entries:
#'   * \code{strength}, a number controlling the strength of the Voting Rights Act
#'   (VRA) constraint, with higher values prioritizing majority-minority districts
#'   over other considerations.
#'   * \code{tgt_vra_min}, the target percentage of minority voters in minority
#'   opportunity districts. Defaults to 0.55.
#'   * \code{tgt_vra_other} The target percentage of minority voters in other
#'   districts. Defaults to 0.25, but should be set to reflect the total minority
#'   population in the state.
#'   * \code{pow_vra}, which controls the allowed deviation from the target
#'   minority percentage; higher values are more tolerant. Defaults to 1.5
#'   * \code{min_pop}, A vector containing the minority population of each
#'   geographic unit.
#' * \code{incumbency}: a list with two entries:
#'   * \code{strength}, a number controlling the tendency of the generated districts
#'   to avoid pairing up incumbents.
#'   * \code{incumbents}, a vector of precinct indices, one for each incumbent's
#'   home address.
#'
#'
#' @param map A \code{\link{redist_map}} object.
#' @param n_sims The number of samples to draw, including warmup.
#' @param warmup The number of warmup samples to discard.
#' @param init The initial state of the map. If not provided, will default to
#'   the reference map of the \code{map} object.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided, the
#' algorithm will only generate maps which split up to \code{ndists-1} counties.
#' If no county-split constraint is desired, this parameter should be left blank.
#' @param compactness Controls the compactness of the generated districts, with
#' higher values preferring more compact districts. Must be nonnegative. See the
#' 'Details' section for more information, and computational considerations.
#' @param constraints A list containing information on constraints to implement.
#' See the 'Details' section for more information.
#' @param constraint_fn A function which takes in a matrix where each column is
#'  a redistricting plan and outputs a vector of log-weights, which will be
#'  added the the final weights.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
#' the algorithm does not appear to be sampling from the target distribution.
#' Must be between 0 and 1.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @return \code{redist_smc} returns an object of class
#'   \code{\link{redist_plans}} containing the simulated plans.
#'
#' @references
#' Carter, D., Herschlag, G., Hunter, Z., and Mattingly, J. (2019). A
#' merge-split proposal for reversible Monte Carlo Markov chain sampling of
#' redistricting plans. arXiv preprint arXiv:1911.01503.
#'
#' DeFord, D., Duchin, M., and Solomon, J. (2019). Recombination: A family of
#' Markov chains for redistricting. arXiv preprint arXiv:1911.05725.
#'
#' @examples \dontrun{
#' data(fl25)
#'
#' fl_map = redist_map(fl25, n_distr=3, pop_tol=0.1)
#'
#' sampled_basic = redist_mergesplit(fl_map, 10000)
#'
#' sampled_constr = redist_mergesplit(fl_map, 10000, constraints=list(
#'                      incumbency = list(strength=1000, incumbents=c(3, 6, 25))
#'                  ))
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_mergesplit = function(map, n_sims, warmup=floor(n_sims/2),
                             init=NULL, counties=NULL, compactness=1,
                             constraints=list(), constraint_fn=function(m) rep(0, ncol(m)),
                             adapt_k_thresh=0.975, verbose=TRUE, silent=FALSE) {
    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)
    n_distr = attr(map, "n_distr")

    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")
    if (n_sims < 1)
        stop("`n_sims` must be positive.")

    if (is.null(init)) init = as.integer(as.factor(get_existing(map)))
    if (is.null(init)) stop("Must provide an initial map.")
    stopifnot(length(init) == V)
    stopifnot(max(init) == n_distr)

    counties = rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(counties)) {
        counties = rep(1, V)
    } else {
        # handle discontinuous counties
        component = contiguity(adj, as.integer(as.factor(counties)))
        counties = dplyr::if_else(component > 1,
                                  paste0(as.character(counties), "-", component),
                                  as.character(counties)) %>%
            as.factor() %>%
            as.integer()
    }

    # Other constraints
    if (is.null(constraints$status_quo))
        constraints$status_quo = list(strength=0, current=rep(1, V))
    if (is.null(constraints$vra))
        constraints$vra = list(strength=0, tgt_vra_min=0.55, tgt_vra_other=0.25,
                               pow_vra=1.5, min_pop=rep(0, V))
    if (is.null(constraints$incumbency))
        constraints$incumbency = list(strength=0, incumbents=integer())

    if (length(constraints$vra$min_pop) != V)
        stop("Length of minority population vector must match the number of units.")
    if (min(constraints$status_quo$current) == 0)
        constraints$status_quo$current = constraints$status_quo$current + 1
    n_current = max(constraints$status_quo$current)

    verbosity = 1
    if (verbose) verbosity = 3
    if (silent) verbosity = 0

    pop_bounds = attr(map, "pop_bounds")
    pop = map[[attr(map, "pop_col")]]

    plans = ms_plans(n_sims, adj, init, counties, pop, n_distr, pop_bounds[2],
                     pop_bounds[1], pop_bounds[3], compactness,
                     constraints$status_quo$strength, constraints$status_quo$current, n_current,
                     constraints$vra$strength, constraints$vra$tgt_vra_min,
                     constraints$vra$tgt_vra_other, constraints$vra$pow_vra, constraints$vra$min_pop,
                     constraints$incumbency$strength, constraints$incumbency$incumbents,
                     adapt_k_thresh, verbosity)

    new_redist_plans(plans[,-1:-warmup], map, "mergesplit",
                     rep(1, n_sims - warmup), FALSE,
                     compactness = compactness,
                     constraints = constraints,
                     adapt_k_thresh = adapt_k_thresh)
}
