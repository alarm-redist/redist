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
#' * \code{hinge}: a list with three entries:
#'   * \code{strength}, a number controlling the strength of the constraint, with
#'   higher values prioritizing districts with group populations at least
#'   \code{tgts_min} over other considerations.
#'   * \code{tgts_min}, the target percentage(s) of minority voters in minority
#'   opportunity districts. Defaults to \code{c(0.55)}.
#'   * \code{min_pop}, A vector containing the minority population of each
#'   geographic unit.
#' * \code{incumbency}: a list with two entries:
#'   * \code{strength}, a number controlling the tendency of the generated districts
#'   to avoid pairing up incumbents.
#'   * \code{incumbents}, a vector of precinct indices, one for each incumbent's
#'   home address.
#' * \code{vra}: a list with five entries, which may be set up using
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
#' * \code{multisplits}: a list with one entry:
#'   * \code{strength}, a number controlling the tendency of the generated districts
#'   to avoid splitting counties multiple times.
#'
#' All constraints are fed into a Gibbs measure, with coefficients on each
#' constraint set by the corresponding \code{strength} parameters.
#' The strength can be any real number, with zero corresponding to no constraint.
#' The \code{status_quo} constraint adds a term measuring the variation of
#' information distance between the plan and the reference, rescaled to \[0, 1\].
#' The \code{hinge} constraint takes a list of target minority percentages. It
#' matches each district to its nearest target percentage, and then applies a
#' penalty of the form \eqn{\sqrt{max(0, tgt - minpct)}}, summing across
#' districts. This penalizes districts which are below their target population.
#' The \code{incumbency} constraint adds a term counting the number of districts
#' containing paired-up incumbents.
#' The \code{vra} constraint  (not recommended) adds a term of the form
#' \eqn{(|tgtvramin-minpct||tgtvraother-minpct|)^{powvra})}, which
#' encourages districts to have minority percentages near either \code{tgt_vra_min}
#' or \code{tgt_vra_other}. This can be visualized with
#' \code{\link{redist.plot.penalty}}.
#'
#'
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to draw.
#' @param counties A vector containing county (or other administrative or
#'   geographic unit) labels for each unit, which may be integers ranging from 1
#'   to the number of counties, or a factor or character vector.  If provided,
#'   the algorithm will only generate maps which split up to \code{ndists-1}
#'   counties. If no county-split constraint is desired, this parameter should
#'   be left blank.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative. See
#'   the 'Details' section for more information, and computational
#'   considerations.
#' @param constraints A list containing information on constraints to implement.
#'   See the 'Details' section for more information.
#' @param resample Whether to perform a final resampling step so that the
#'   generated plans can be used immediately.  Set this to \code{FALSE} to
#'   perform direct importance sampling estimates, or to adjust the weights
#'   manually.
#' @param constraint_fn A function which takes in a matrix where each column is
#'  a redistricting plan and outputs a vector of log-weights, which will be
#'  added the the final weights.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#'   value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if the
#'   algorithm does not appear to be sampling from the target distribution. Must
#'   be between 0 and 1.
#' @param seq_alpha The amount to adjust the weights by at each resampling step;
#'   higher values prefer exploitation, while lower values prefer exploration.
#'   Must be between 0 and 1.
#' @param truncate Whether to truncate the importance sampling weights at the
#'   final step by \code{trunc_fn}.  Recommended if \code{compactness} is not 1.
#'   Truncation only applied if \code{resample=TRUE}.
#' @param trunc_fn A function which takes in a vector of weights and returns a
#'   truncated vector. If \code{\link[loo]{loo}} package is installed (strongly
#'   recommended), will default to Pareto-smoothed Importance Sampling (PSIS)
#'   rather than naive truncation.
#' @param pop_temper The strength of the automatic population tempering. Try
#'   values of 0.01-0.05 to start if the algorithm gets stuck on the final few
#'   splits.
#' @param ref_name a name for the existing plan, which will be added as a
#'   reference plan, or \code{FALSE} to not include the initial plan in the
#'   output. Defaults to the column name of the existing plan.
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
#' @examples \donttest{
#' set.seed(1)
#' data(fl25)
#'
#' fl_map = redist_map(fl25, ndists=3, pop_tol=0.1)
#'
#' sampled_basic = redist_smc(fl_map, 10000)
#'
#' sampled_constr = redist_smc(fl_map, 10000, constraints=list(
#'                                 incumbency = list(strength=100, incumbents=c(3, 6, 25))
#'                             ))
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_smc = function(map, nsims, counties=NULL, compactness=1, constraints=list(),
                      resample=TRUE, constraint_fn=function(m) rep(0, ncol(m)),
                      adapt_k_thresh=0.975, seq_alpha=0.2+0.3*compactness,
                      truncate=(compactness != 1), trunc_fn=redist_quantile_trunc,
                      pop_temper=0, ref_name=NULL, verbose=TRUE, silent=FALSE) {
    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)

    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")
    if (seq_alpha <= 0 | seq_alpha > 1)
        stop("`seq_alpha` parameter must lie in (0, 1].")
    if (nsims < 1)
        stop("`nsims` must be positive.")

    counties = rlang::eval_tidy(rlang::enquo(counties), map)
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
        if (any(component > 1)) {
            warning('counties were not contiguous; expect additional splits.')
        }
    }

    # Other constraints
    constraints = eval_tidy(enquo(constraints), map)
    proc = process_smc_ms_constr(constraints, V)
    constraints = proc$constraints
    n_current = max(constraints$status_quo$current)

    verbosity = 1
    if (verbose) verbosity = 3
    if (silent) verbosity = 0

    pop_bounds = attr(map, "pop_bounds")
    pop = map[[attr(map, "pop_col")]]
    ndists = attr(map, "ndists")
    if (any(pop >= get_target(map)))
        stop("Units ", which(pop >= get_target(map)),
             " have population larger than the district target.\n",
             "Redistricting impossible.")

    lp = rep(0, nsims)
    plans = smc_plans(nsims, adj, counties, pop, ndists, pop_bounds[2],
                      pop_bounds[1], pop_bounds[3], compactness,
                      constraints$status_quo$strength, constraints$status_quo$current, n_current,
                      constraints$vra$strength, constraints$vra$tgt_vra_min,
                      constraints$vra$tgt_vra_other, constraints$vra$pow_vra, proc$min_pop,
                      constraints$hinge$strength, constraints$hinge$tgts_min,
                      constraints$incumbency$strength, constraints$incumbency$incumbents,
                      constraints$multisplits$strength,
                      lp, adapt_k_thresh, seq_alpha, pop_temper, verbosity);


    lr = -lp + constraint_fn(plans)
    wgt = exp(lr - mean(lr))
    wgt = wgt / mean(wgt)
    n_eff = length(wgt) * mean(wgt)^2 / mean(wgt^2)

    if (resample) {
        if (!truncate) {
            mod_wgt = wgt
        } else if (requireNamespace("loo", quietly=TRUE) && missing(trunc_fn)) {
            mod_wgt = wgt / sum(wgt)
            mod_wgt = loo::weights.importance_sampling(
                loo::psis(log(mod_wgt), r_eff=NA), log=FALSE)
        } else {
            mod_wgt = trunc_fn(wgt)
        }
        n_eff = length(mod_wgt) * mean(mod_wgt)^2 / mean(mod_wgt^2)
        mod_wgt = mod_wgt / sum(mod_wgt)

        plans = plans[, sample(nsims, nsims, replace=T, prob=mod_wgt), drop=FALSE]
    }

    if (!is.nan(n_eff) && n_eff/nsims <= 0.05)
        warning("Less than 5% resampling efficiency. Consider weakening constraints and/or adjusting `seq_alpha`.")

    out = new_redist_plans(plans, map, "smc", wgt, resample,
                     n_eff = n_eff,
                     compactness = compactness,
                     constraints = constraints,
                     adapt_k_thresh = adapt_k_thresh,
                     seq_alpha = seq_alpha,
                     pop_temper = pop_temper)

    exist_name = attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name)) {
        ref_name = if (!is.null(ref_name)) ref_name else exist_name
        out = add_reference(out, map[[exist_name]], ref_name)
    }

    out
}


# Helper constraint processor.
# Constraint defaults contained HERE.
#
# @param constraints passed into `redist_smc` or `redist_ms`
#
# @return a list with new `constraints` and a minority population vector `min_pop`
process_smc_ms_constr = function(constraints, V) {
    defaults = list(
        status_quo = list(strength=0, current=rep(1, V)),
        hinge = list(strength=0, tgts_min=0.55, min_pop=NULL),
        vra = list(strength=0, tgt_vra_min=0.55, tgt_vra_other=0.25,
                   pow_vra=1.5, min_pop=integer()),
        incumbency = list(strength=0, incumbents=integer()),
        splits = list(strength=0),
        multisplits = list(strength=0)
    )

    for (type in names(constraints)) {
        for (el in names(constraints[[type]])) {
            defaults[[type]][[el]] = constraints[[type]][[el]]
        }
    }

    if (min(defaults$status_quo$current) == 0)
        defaults$status_quo$current = defaults$status_quo$current + 1

    min_pop = rep(0, V)
    if (defaults$hinge$strength > 0) {
        if (defaults$vra$strength > 0)
            stop("Specify one of `vra` or `vra_old` constraints, not both")
        min_pop = defaults$hinge$min_pop
    } else if (defaults$vra$strength > 0) {
        min_pop = defaults$vra$min_pop
    }
    if (length(min_pop) != V)
        stop("Length of minority population vector must match the number of units.")

    list(constraints = defaults,
         min_pop = min_pop)
}

#' Helper function to truncate importance weights
#'
#' Defined as \code{pmin(x, quantile(x, 1 - length(x)^(-0.5)))}
#'
#' @param x the weights
#'
#' @return numeric vector
#'
#' @export
#'
#' @examples
#' redist_quantile_trunc(c(1,2,3,4))
#'
redist_quantile_trunc = function(x) pmin(x, quantile(x, 1 - length(x)^(-0.5)))
