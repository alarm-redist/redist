#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2020/07/08
# Purpose: R wrapper to run SMC redistricting code
####################################################

#' SMC Redistricting Sampler
#'
#' \code{redist.smc} uses a Sequential Monte Carlo algorithm to
#' generate nearly independent congressional or legislative redistricting
#' plans according to contiguity, population, compactness, and administrative
#' boundary constraints.
#'
#' This function draws nearly-independent samples from a specific target measure,
#' controlled by the \code{pop_tol}, \code{compactness}, \code{constraints}, and
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
#' importance sampling weights.  By default these weights are truncated at
#' \code{nsims^0.04 / 100} to stabilize the resulting estimates, but if truncation
#' is used, a specific truncation function should probably be chosen by the user.
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
#' * \code{vra}: a list with three entries:
#'   * \code{strength}, a number controlling the strength of the Voting Rights Act
#'   (VRA) constraint, with higher values prioritizing majority-minority districts
#'   over other considerations.
#'   * \code{tgts_min}, the target percentage(s) of minority voters in minority
#'   opportunity districts. Defaults to \code{c(0.55)}.
#'   * \code{min_pop}, A vector containing the minority population of each
#'   geographic unit.
#' * \code{incumbency}: a list with two entries:
#'   * \code{strength}, a number controlling the tendency of the generated districts
#'   to avoid pairing up incumbents.
#'   * \code{incumbents}, a vector of precinct indices, one for each incumbent's
#'   home address.
#' * \code{vra_old}: a list with five entries, which may be set up using
#'   \code{\link{redist.constraint.helper}}:
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
#'
#' All constraints are fed into a Gibbs measure, with coefficients on each
#' constraint set by the corresponding \code{strength} parameters.
#' The \code{status_quo} constraint adds a term measuring the variation of
#' information distance between the plan and the reference, rescaled to \[0, 1\].
#' The \code{vra} constraint takes a list of target minority percentages. It
#' matches each district to its nearest target percentage, and then applies a
#' penalty of the form \eqn{\sqrt{max(0, tgt - minpct)}}, summing across
#' districts. This penalizes districts which are below their target population.
#' The \code{incumbency} constraint adds a term counting the number of districts
#' containing paired-up incumbents.
#' The \code{vra_old} constraint adds a term of the form
#' \eqn{(|tgtvramin-minpct||tgtvraother-minpct|)^{powvra})}, which
#' encourages districts to have minority percentages near either \code{tgt_vra_min}
#' or \code{tgt_vra_other}. This can be visualized with
#' \code{\link{redist.plot.penalty}}.
#'
#' @param adj An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param total_pop A vector containing the populations of each geographic unit.
#' @param nsims The number of samples to draw.
#' @param ndists The number of districts in each redistricting plan.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which must  be integers ranging from 1
#' to the number of counties.  If provided, the algorithm will only generate
#' maps which split up to \code{ndists-1} counties.  If no county-split
#' constraint is desired, this parameter should be left blank.
#' @param pop_tol The desired population constraint.  All sampled districts
#' will have a deviation from the target district size no more than this value
#' in percentage terms, i.e., \code{pop_tol=0.01} will ensure districts have
#' populations within 1% of the target population.
#' @param compactness Controls the compactness of the generated districts, with
#' higher values preferring more compact districts. Must be nonnegative. See the
#' 'Details' section for more information, and computational considerations.
#' @param constraints A list containing information on constraints to implement.
#' See the 'Details' section for more information.
#' @param resample Whether to perform a final resampling step so that the
#' generated plans can be used immediately.  Set this to \code{FALSE} to perform
#' direct importance sampling estimates, or to adjust the weights manually.
#' @param pop_bounds A numeric vector with three elements \code{c(lower, target, upper)}
#' providing more precise population bounds for the algorithm. Districts
#' will have population between \code{lower} and \code{upper}, with a goal of
#' \code{target}.  If set, overrides \code{popcons}.
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
#' @param trunc_fn A function which takes in a vector of weights and returns
#' a truncated vector. Recommended to specify this manually if truncating weights.
#' @param pop_temper The strength of the automatic population tempering.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended.
#' @param silent Whether to supress all diagnostic information.
#' @param adjobj Deprecated, use adj. An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param popvec Deprecated, use total_pop. A vector containing the populations of each geographic unit.
#' @param popcons The desired population constraint.  All sampled districts
#' will have a deviation from the target district size no more than this value
#' in percentage terms, i.e., \code{popcons=0.01} will ensure districts have
#' populations within 1% of the target population.
#'
#' @return \code{redist.smc} returns an object of class \code{redist}, which
#' is a list containing the following components:
#' \item{aList}{The adjacency list used to sample}
#' \item{cdvec}{The matrix of sampled plans. Each row is a geographical unit,
#' and each column is a sample.}
#' \item{wgt}{The importance sampling weights, normalized to sum to 1.}
#' \item{orig_wgt}{The importance sampling weights before resampling or truncation, normalized to have mean 1.}
#' \item{nsims}{The number of plans sampled.}
#' \item{pct_dist_parity}{The population constraint.}
#' \item{compactness}{The compactness constraint.}
#' \item{counties}{The computed constraint options list (see above).}
#' \item{maxdev}{The maximum population deviation of each sample.}
#' \item{total_pop}{The provided vector of unit populations.}
#' \item{counties}{The provided county vector.}
#' \item{adapt_k_thresh}{The provided control parameter.}
#' \item{seq_alpha}{The provided control vector.}
#' \item{algorithm}{The algorithm used, here \code{"smc"}.}
#'
#' @references
#' McCartan, C., & Imai, K. (2020). Sequential Monte Carlo for Sampling Balanced and Compact Redistricting Plans.
#' Available at \url{https://imai.fas.harvard.edu/research/files/SMCredist.pdf}.
#'
#' @examples \dontrun{
#' data(fl25)
#' data(fl25_adj)
#' data(fl25_enum)
#'
#' sampled_basic = redist.smc(fl25_adj, fl25$pop,
#'                            nsims=10000, ndists=3, pop_tol=0.1)
#'
#' sampled_constr = redist.smc(fl25_adj, fl25$pop,
#'                             nsims=10000, ndists=3, pop_tol=0.1,
#'                             constraints=list(
#'                                 status_quo = list(strength=10, current=fl25_enum$plans[,5118]),
#'                                 incumbency = lsit(strength=1000, incumbents=c(3, 6, 25))
#'                             ))
#' }
#'
#' @concept simulate
#' @md
#' @export
redist.smc = function(adj, total_pop, nsims, ndists, counties=NULL,
                      pop_tol = 0.01, pop_bounds=NULL, compactness=1,
                      constraints=list(),
                      resample=TRUE,
                      constraint_fn=function(m) rep(0, ncol(m)),
                      adapt_k_thresh=0.975, seq_alpha=0.2+0.2*compactness,
                      truncate=(compactness != 1),
                      trunc_fn=function(x) pmin(x, 0.01*nsims^0.4),
                      pop_temper=0, verbose=TRUE, silent=FALSE,
                      adjobj, popvec, popcons) {
    if (!missing(adjobj)) {
        .Deprecated(new = 'adj', old = 'adjobj')
        adj <- adjobj
    }
    if (!missing(popvec)) {
        .Deprecated(new = 'total_pop', old = 'popvec')
        total_pop <- popvec
    }
    if (!missing(popcons)) {
        .Deprecated(new = 'pop_tol', old = 'popcons')
        pop_tol <- popcons
    }


    V = length(total_pop)

    if (missing(adj)) stop("Please supply adjacency matrix or list")
    if (missing(total_pop)) stop("Please supply vector of geographic unit populations")
    if (missing(nsims)) stop("Please supply number of simulations to run algorithm")
    if (pop_tol <= 0) stop("Population constraint must be positive")
    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")
    if (seq_alpha <= 0 | seq_alpha > 1)
        stop("`seq_alpha` parameter must lie in (0, 1].")
    if (nsims < 1)
        stop("`nsims` must be positive.")

    if (is.null(counties)) {
        counties = rep(1, V)
    } else {
        if (length(unique(counties)) != max(counties))
            stop("County numbers must run from 1 to n_county with no interruptions.")

        # handle discontinuous counties
        counties = redist.county.relabel(adj, counties)
        counties <- redist.county.id(counties)
    }

    # Population bounds
    if (!is.null(pop_bounds)) {
        if (length(pop_bounds) != 3)
            stop("`pop_bounds` must be of the form c(lower, target, upper).")
        if (!all(diff(pop_bounds) > 0))
            stop("`pop_bounds` must satisfy lower < target < upper.")
    } else {
        target = sum(total_pop) / ndists
        pop_bounds = target * c(1 - pop_tol, 1, 1 + pop_tol)
    }

    # Other constraints
    proc = process_smc_ms_constr(constraints, V)
    constraints = proc$constraints
    n_current = max(constraints$status_quo$current)

    # sanity-check everything
    preproc = redist.preproc(adj, total_pop, rep(0, V), ndists, pop_tol,
                             temper="none", constraint="none")
    adjlist = preproc$data$adjlist
    class(adjlist) = "list"

    verbosity = 1
    if (verbose) verbosity = 3
    if (silent) verbosity = 0

    lp = rep(0, nsims)
    maps = smc_plans(nsims, adjlist, counties, total_pop, ndists, pop_bounds[2],
                     pop_bounds[1], pop_bounds[3], compactness,
                     constraints$status_quo$strength, constraints$status_quo$current, n_current,
                     constraints$vra_old$strength, constraints$vra_old$tgt_vra_min,
                     constraints$vra_old$tgt_vra_other, constraints$vra_old$pow_vra, proc$min_pop,
                     constraints$vra$strength, constraints$vra$tgts_min,
                     constraints$incumbency$strength, constraints$incumbency$incumbents,
                     lp, adapt_k_thresh, seq_alpha, pop_temper, verbosity);

    dev = max_dev(maps, total_pop, ndists)
    maps = maps

    lr = -lp + constraint_fn(maps)
    wgt = exp(lr - mean(lr))
    wgt = wgt / mean(wgt)
    orig_wgt = wgt
    if (truncate)
        wgt = trunc_fn(wgt)
    wgt = wgt/sum(wgt)
    n_eff = length(wgt) * mean(wgt)^2 / mean(wgt^2)

    if (n_eff/nsims <= 0.05)
        warning("Less than 5% resampling efficiency. Consider weakening constraints and/or adjusting `seq_alpha`.")

    if (resample) {
        maps = maps[, sample(nsims, nsims, replace=T, prob=wgt)]
        wgt = rep(1/nsims, nsims)
    }

    algout = list(
        adj = adjlist,
        plans = maps,
        wgt = wgt,
        orig_wgt = orig_wgt,
        nsims = nsims,
        n_eff = n_eff,
        pct_dist_parity = pop_tol,
        compactness = compactness,
        constraints = constraints,
        maxdev = dev,
        total_pop = total_pop,
        counties = if (max(counties)==1) NULL else counties,
        adapt_k_thresh = adapt_k_thresh,
        seq_alpha = seq_alpha,
        algorithm="smc"
    )
    class(algout) = "redist"

    algout
}

#' Confidence Intervals for Importance Sampling Estimates
#'
#' Builds a confidence interval for a quantity of interest,
#' given importance sampling weights.
#'
#' @param x A numeric vector containing the quantity of interest
#' @param wgt A numeric vector containing the nonnegative importance weights.
#'   Will be normalized automatically.
#' @param conf The confidence level for the interval.
#'
#' @returns A two-element vector of the form [lower, upper] containing
#' the importance sampling confidence interval.
#'
#' @concept post
#' @export
redist.smc_is_ci = function(x, wgt, conf=0.99) {
    wgt = wgt / sum(wgt)
    mu = sum(x*wgt)
    sig = sqrt(sum((x - mu)^2 * wgt^2))
    mu + qnorm(c((1-conf)/2, 1-(1-conf)/2))*sig
}

