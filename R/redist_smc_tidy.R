#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: tidy R wrapper to run SMC redistricting code
####################################################

#' SMC Redistricting Sampler
#'
#' \code{redist_smc} uses a Sequential Monte Carlo algorithm to
#' generate nearly independent congressional or legislative redistricting
#' plans according to contiguity, population, compactness, and administrative
#' boundary constraints.
#'
#' This function draws nearly-independent samples from a specific target measure,
#' controlled by the \code{popcons}, \code{compactness}, \code{constraints}, and
#' \code{constraint_fn} parameters.
#'
#' Key to ensuring good performance is monitoring the efficiency of the resampling
#' process at each SMC stage.  Unless \code{silent=F}, this function will print
#' out the effective sample size of each resampling step to allow the user to
#' monitor the efficiency.  If \code{verbose=T} the function will also print
#' out information on the \eqn{k_i} values automatically chosen and the
#' acceptance rate (based on the population constraint) at each step.
#'
#' Higher values of \code{compactness} sample more compact districts;
#' setting this parameter to 1 is computationally efficient and generates nicely
#' compact districts.  Values of other than 1 may lead to highly variable
#' importance sampling weights.  By default these weights are truncated using
#' \code{\link{redist_quantile_trunc}} to stabilize the resulting estimates, but
#' if truncation is used, a specific truncation function should probably be
#' chosen by the user.
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
#' @param n_sims The number of samples to draw.
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
#' @param resample Whether to perform a final resampling step so that the
#' generated plans can be used immediately.  Set this to \code{FALSE} to perform
#' direct importance sampling estimates, or to adjust the weights manually.
#' @param constraint_fn A function which takes in a matrix where each column is
#'  a redistricting plan and outputs a vector of log-weights, which will be
#'  added the the final weights.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
#' the algorithm does not appear to be sampling from the target distribution.
#' Must be between 0 and 1.
#' @param seq_alpha The amount to adjust the weights by at each resampling step;
#' higher values prefer exploitation, while lower values prefer exploration.
#' Must be between 0 and 1.
#' @param truncate Whether to truncate the importance sampling weights at the
#' final step by \code{trunc_fn}.  Recommended if \code{compactness} is not 1.
#' Truncation only applied if \code{resample=TRUE}.
#' @param trunc_fn A function which takes in a vector of weights and returns
#' a truncated vector. If \code{\link[loo]{loo}} package is installed (strongly
#' recommended), will default to Pareto-smoothed Importance Sampling (PSIS)
#' rather than naive truncation.
#' @param pop_temper The strength of the automatic population tempering.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @return \code{redist_smc} returns an object of class
#'   \code{\link{redist_plans}} containing the simulated plans.
#'
#' @references
#' McCartan, C., & Imai, K. (2020). Sequential Monte Carlo for Sampling Balanced and Compact Redistricting Plans.
#' Available at \url{https://imai.fas.harvard.edu/research/files/SMCredist.pdf}.
#'
#' @examples \dontrun{
#' data(fl25)
#'
#' fl_map = redist_map(fl25, n_distr=3, pop_tol=0.1)
#'
#' sampled_basic = redist_smc(fl_map, 10000)
#'
#' sampled_constr = redist_smc(fl_map, 10000, constraints=list(
#'                                 incumbency = list(strength=1000, incumbents=c(3, 6, 25))
#'                             ))
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_smc = function(map, n_sims, counties=NULL, compactness=1, constraints=list(),
                      resample=TRUE, constraint_fn=function(m) rep(0, ncol(m)),
                      adapt_k_thresh=0.975, seq_alpha=0.2+0.3*compactness,
                      truncate=(compactness != 1), trunc_fn=redist_quantile_trunc,
                      pop_temper=0, verbose=TRUE, silent=FALSE) {
    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)

    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")
    if (seq_alpha <= 0 | seq_alpha > 1)
        stop("`seq_alpha` parameter must lie in (0, 1].")
    if (n_sims < 1)
        stop("`n_sims` must be positive.")

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
    n_distr = attr(map, "n_distr")

    lp = rep(0, n_sims)
    plans = smc_plans(n_sims, adj, counties, pop, n_distr, pop_bounds[2],
                      pop_bounds[1], pop_bounds[3], compactness,
                      constraints$status_quo$strength, constraints$status_quo$current, n_current,
                      constraints$vra$strength, constraints$vra$tgt_vra_min,
                      constraints$vra$tgt_vra_other, constraints$vra$pow_vra, constraints$vra$min_pop,
                      constraints$incumbency$strength, constraints$incumbency$incumbents,
                      lp, adapt_k_thresh, seq_alpha, pop_temper, verbosity);


    lr = -lp + constraint_fn(plans)
    wgt = exp(lr - mean(lr))
    wgt = wgt / mean(wgt)
    n_eff = length(wgt) * mean(wgt)^2 / mean(wgt^2)

    if (resample) {
        if (!truncate) {
            mod_wgt = wgt
        } else if (suppressMessages(require("loo")) && missing(trunc_fn)) {
            mod_wgt = wgt / sum(wgt)
            mod_wgt = loo::weights.importance_sampling(
                loo::psis(log(mod_wgt), r_eff=NA), log=FALSE)
        } else {
            mod_wgt = trunc_fn(wgt)
        }
        n_eff = length(mod_wgt) * mean(mod_wgt)^2 / mean(mod_wgt^2)
        mod_wgt = mod_wgt / sum(mod_wgt)

        plans = plans[, sample(n_sims, n_sims, replace=T, prob=mod_wgt)]
    }

    if (n_eff/n_sims <= 0.05)
        warning("Less than 5% resampling efficiency. Consider weakening constraints and/or adjusting `seq_alpha`.")

    new_redist_plans(plans, map, "smc", wgt, resample,
                     n_eff = n_eff,
                     compactness = compactness,
                     constraints = constraints,
                     adapt_k_thresh = adapt_k_thresh,
                     seq_alpha = seq_alpha,
                     pop_temper = pop_temper)
}

#' Helper function to truncate importance weights
#'
#' Defined as \code{pmin(x, quantile(x, 1 - length(x)^(-0.5)))}
#'
#' @param x the weights
#'
#' @export
redist_quantile_trunc = function(x) pmin(x, quantile(x, 1 - length(x)^(-0.5)))
