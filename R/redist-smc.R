#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2020/07/08
# Purpose: R wrapper to run SMC redistricting code
####################################################

#' @rdname redist_smc
#' @order 2
#'
#' @param adj An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param total_pop A vector containing the populations of each geographic unit.
#' @param nsims The number of samples to draw.
#' @param ndists The number of districts in each redistricting plan.
#' @param pop_tol The desired population constraint.  All sampled districts
#' will have a deviation from the target district size no more than this value
#' in percentage terms, i.e., \code{pop_tol=0.01} will ensure districts have
#' populations within 1% of the target population.
#' @param pop_bounds A numeric vector with three elements \code{c(lower, target, upper)}
#' providing more precise population bounds for the algorithm. Districts
#' will have population between \code{lower} and \code{upper}, with a goal of
#' \code{target}.  If set, overrides \code{pop_tol}.
#'
#' @return \code{redist.smc} (Deprecated) returns an object of class \code{redist}, which
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
                      pop_temper=0, verbose=TRUE, silent=FALSE) {
    .Deprecated("redist_smc")
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
        if (any(is.na(counties)))
            stop("County vector must not contain missing values.")

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
                     constraints$vra$strength, constraints$vra$tgt_vra_min,
                     constraints$vra$tgt_vra_other, constraints$vra$pow_vra, proc$min_pop,
                     constraints$hinge$strength, constraints$hinge$tgts_min,
                     constraints$incumbency$strength, constraints$incumbency$incumbents,
                     constraints$multisplits$strength,
                     lp, adapt_k_thresh, seq_alpha, pop_temper, verbosity);

    dev = max_dev(maps, total_pop, ndists)
    maps = maps

    lr = -lp + constraint_fn(maps)
    wgt = exp(lr - mean(lr, na.rm=TRUE))
    wgt = wgt / mean(wgt, na.rm=TRUE)
    orig_wgt = wgt
    if (truncate)
        wgt = trunc_fn(wgt)
    wgt = wgt/sum(wgt, na.rm=TRUE)
    n_eff = length(wgt) * mean(wgt, na.rm=TRUE)^2 / mean(wgt^2, na.rm=TRUE)
    if (is.nan(n_eff))
        warning("Some invalid plans were generated.")

    if (!is.nan(n_eff) && n_eff/nsims <= 0.05)
        warning("Less than 5% resampling efficiency. Consider weakening constraints and/or adjusting `seq_alpha`.")

    if (resample) {
        maps = maps[, sample(nsims, nsims, replace=TRUE, prob=wgt)]
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

