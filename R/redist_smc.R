#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: tidy R wrapper to run SMC redistricting code
####################################################

#' SMC Redistricting Sampler
#'
#' `redist_smc` uses a Sequential Monte Carlo algorithm to
#' generate nearly independent congressional or legislative redistricting
#' plans according to contiguity, population, compactness, and administrative
#' boundary constraints.
#'
#' This function draws nearly-independent samples from a specific target measure,
#' controlled by the `map`, `compactness`, and `constraints` parameters.
#'
#' Key to ensuring good performance is monitoring the efficiency of the resampling
#' process at each SMC stage.  Unless `silent=FALSE`, this function will print
#' out the effective sample size of each resampling step to allow the user to
#' monitor the efficiency.  If `verbose=TRUE` the function will also print
#' out information on the \eqn{k_i} values automatically chosen and the
#' acceptance rate (based on the population constraint) at each step.
#' Users should also check the [plans_diversity()] of the sample.
#'
#' Higher values of `compactness` sample more compact districts;
#' setting this parameter to 1 is computationally efficient and generates nicely
#' compact districts.  Values of other than 1 may lead to highly variable
#' importance sampling weights.  In these cases, these weights are by default
#' truncated using [redist_quantile_trunc()] to stabilize the resulting
#' estimates, but if truncation is used, a specific truncation function should
#' probably be chosen by the user.
#'
#' @param map A [redist_map()] object.
#' @param nsims The number of samples to draw.
#' @param counties A vector containing county (or other administrative or
#'   geographic unit) labels for each unit, which may be integers ranging from 1
#'   to the number of counties, or a factor or character vector.  If provided,
#'   the algorithm will only generate maps which split up to `ndists-1`
#'   counties. If no county-split constraint is desired, this parameter should
#'   be left blank.
#' @param compactness Controls the compactness of the generated districts, with
#'   higher values preferring more compact districts. Must be nonnegative. See
#'   the 'Details' section for more information, and computational
#'   considerations.
#' @param constraints A [redist_constr()] object or a list containing
#'   information on sampling constraints. See [constraints] for more information.
#' @param resample Whether to perform a final resampling step so that the
#'   generated plans can be used immediately.  Set this to `FALSE` to
#'   perform direct importance sampling estimates, or to adjust the weights
#'   manually.
#' @param constraint_fn (Deprecated) A function which takes in a matrix where
#'   each column is a redistricting plan and outputs a vector of log-weights,
#'   which will be added the the final weights.
#' @param init_particles A matrix of partial plans to begin sampling from. For
#'  advanced use only.  The matrix must have `nsims` columns and a row for
#'  every precinct. It is important to ensure that the existing districts meet
#'  contiguity and population constraints, or there may be major issues when
#'  sampling.
#' @param n_steps How many steps to run the SMC algorithm for.
#'   Each step splits off a new district. Defaults to all remaining districts.
#'   If fewer than the number of remaining splits, reference plans are disabled.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#'   value `k_i` for each splitting iteration. Set to 0.9999 or 1 if the
#'   algorithm does not appear to be sampling from the target distribution. Must
#'   be between 0 and 1.
#' @param seq_alpha The amount to adjust the weights by at each resampling step;
#'   higher values prefer exploitation, while lower values prefer exploration.
#'   Must be between 0 and 1.
#' @param truncate Whether to truncate the importance sampling weights at the
#'   final step by `trunc_fn`.  Recommended if `compactness` is not 1.
#'   Truncation only applied if `resample=TRUE`.
#' @param trunc_fn A function which takes in a vector of weights and returns a
#'   truncated vector. If the [loo][loo::loo] package is installed (strongly
#'   recommended), will default to Pareto-smoothed Importance Sampling (PSIS)
#'   rather than naive truncation.
#' @param pop_temper The strength of the automatic population tempering. Try
#'   values of 0.01-0.05 to start if the algorithm gets stuck on the final few
#'   splits.
#' @param final_infl A multiplier for the population constraint on the final
#'   iteration. Used to loosen the constraint when the sampler is getting stuck
#'   on the final split.
#' @param ref_name a name for the existing plan, which will be added as a
#'   reference plan, or `FALSE` to not include the initial plan in the
#'   output. Defaults to the column name of the existing plan.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @return `redist_smc` returns an object of class
#'   [redist_plans()] containing the simulated plans.
#'
#' @references
#' McCartan, C., & Imai, K. (2020). Sequential Monte Carlo for Sampling Balanced and Compact Redistricting Plans.
#' Available at \url{https://imai.fas.harvard.edu/research/files/SMCredist.pdf}.
#'
#' @examples \donttest{
#' data(fl25)
#'
#' fl_map = redist_map(fl25, ndists=3, pop_tol=0.1)
#'
#' sampled_basic = redist_smc(fl_map, 10000)
#'
#' constr = redist_constr(fl_map)
#' constr = add_constr_incumbency(constr, strength=100, incumbents=c(3, 6, 25))
#' sampled_constr = redist_smc(fl_map, 10000, constraints=constr)
#' }
#'
#' @concept simulate
#' @md
#' @order 1
#' @export
redist_smc = function(map, nsims, counties=NULL, compactness=1, constraints=list(),
                      resample=TRUE, constraint_fn=function(m) rep(0, ncol(m)),
                      init_particles=NULL, n_steps=NULL,
                      adapt_k_thresh=0.975, seq_alpha=0.2+0.3*compactness,
                      truncate=(compactness != 1), trunc_fn=redist_quantile_trunc,
                      pop_temper=0, final_infl=1, ref_name=NULL,
                      verbose=TRUE, silent=FALSE) {
    if (!missing(constraint_fn)) cli_warn("{.arg constraint_fn} is deprecated.")

    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)

    if (compactness < 0)
        cli_abort("{.arg compactness} must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")
    if (seq_alpha <= 0 | seq_alpha > 1)
        cli_abort("{.arg seq_alpha} must lie in (0, 1].")
    if (nsims < 1)
        cli_abort("{.arg nsims} must be positive.")

    counties = rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(counties)) {
        counties = rep(1, V)
    } else {
        if (any(is.na(counties)))
            cli_abort("County vector must not contain missing values.")

        # handle discontinuous counties
        component = contiguity(adj, as.integer(as.factor(counties)))
        counties = dplyr::if_else(component > 1,
                                  paste0(as.character(counties), "-", component),
                                  as.character(counties)) %>%
            as.factor() %>%
            as.integer()
        if (any(component > 1)) {
            cli_warn("Counties were not contiguous; expect additional splits.")
        }
    }

    # Other constraints
    if (!inherits(constraints, "redist_constr")) {
        constraints = new_redist_constr(eval_tidy(enquo(constraints), map))
    }
    constraints = as.list(constraints) # drop data attribute

    verbosity = 1
    if (verbose) verbosity = 3
    if (silent) verbosity = 0

    pop_bounds = attr(map, "pop_bounds")
    pop = map[[attr(map, "pop_col")]]
    ndists = attr(map, "ndists")
    if (any(pop >= pop_bounds[3]))
        cli_abort(c("Units {which(pop >= pop_bounds[3])} have population larger than the district target.",
                    "x"="Redistricting impossible."))

    # handle particle inits
    if (is.null(init_particles)) {
        init_particles = matrix(0L, nrow=V, ncol=nsims)
        n_drawn = 0L
    } else {
        if (nrow(init_particles) != V)
            cli_abort("{.arg init_particles} must have as many rows as {.arg map} has precincts.")
        if (ncol(init_particles) != nsims)
            cli_abort("{.arg init_particles} must have {.arg nsims} columns.")
        n_drawn = as.integer(max(init_particles[, 1]))
    }
    if (is.null(n_steps)) {
        n_steps = attr(map, "ndists") - n_drawn - 1L
    }
    final_dists = n_drawn + n_steps + 1L
    if (final_dists > ndists) {
        cli_abort("Too many districts already drawn to take {n_steps} steps.")
    }

    lp = rep(0, nsims)
    plans = smc_plans(nsims, adj, counties, pop, ndists, pop_bounds[2],
                      pop_bounds[1], pop_bounds[3], compactness,
                      init_particles, n_drawn, n_steps, constraints,
                      lp, adapt_k_thresh, seq_alpha, pop_temper, final_infl, verbosity);

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

        plans = plans[, sample(nsims, nsims, replace=TRUE, prob=mod_wgt), drop=FALSE]
    }

    if (!is.nan(n_eff) && n_eff/nsims <= 0.05)
        cli_warn(c("Less than 5% resampling efficiency.",
                   ">"="Consider weakening constraints and/or adjusting {.arg seq_alpha}."))

    out = new_redist_plans(plans, map, "smc", wgt, resample,
                           n_eff = n_eff,
                           compactness = compactness,
                           constraints = constraints,
                           ndists = final_dists,
                           adapt_k_thresh = adapt_k_thresh,
                           seq_alpha = seq_alpha,
                           pop_temper = pop_temper)

    exist_name = attr(map, "existing_col")
    if (!is.null(exist_name) && !isFALSE(ref_name) && ndists == final_dists) {
        ref_name = if (!is.null(ref_name)) ref_name else exist_name
        out = add_reference(out, map[[exist_name]], ref_name)
    }

    out
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
