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
#' but not identical to those used in the references.  1-level hierarchical
#' Merge-split is supported through the \code{counties} parameter; unlike in
#' the SMC algorithm, this does not guarantee a maximum number of county splits.
#'
#' This function draws samples from a specific target measure, controlled by the
#' \code{compactness}, \code{constraints}, and \code{constraint_fn} parameters.
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
#' * \code{hinge}: a list with three entries:
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
#' * \code{splits}: a list with one entry:
#'   * \code{strength}, a number controlling the tendency of the generated districts
#'   to avoid splitting counties.
#' * \code{multisplits}: a list with one entry:
#'   * \code{strength}, a number controlling the tendency of the generated districts
#'   to avoid splitting counties multiple times.
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
#' containing paired-up incumbents. The \code{splits} constraint adds a term
#' counting the number of counties which contain precincts belonging to more
#' than one district.
#' The \code{vra} constraint (not recommended) adds a term of the form
#' \eqn{(|tgtvramin-minpct||tgtvraother-minpct|)^{powvra})}, which
#' encourages districts to have minority percentages near either \code{tgt_vra_min}
#' or \code{tgt_vra_other}. This can be visualized with
#' \code{\link{redist.plot.penalty}}.
#'
#'
#'
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to draw, including warmup.
#' @param warmup The number of warmup samples to discard.
#' @param init_plan The initial state of the map. If not provided, will default to
#'   the reference map of the \code{map} object, or if none exists, will sample
#'   a random initial state using \code{\link{redist_smc}}. You can also request
#'   a random initial state by setting \code{init_plan="sample"}.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided, the
#' algorithm will generate maps tend to follow county lines.  You may combine this
#' with a Gibbs constraint on the number of county splits using the
#' \code{constraints} parameter; see below. If no county-split considerations
#' are desired, this parameter should be left blank.
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
#' @param k The number of edges to consider cutting after drawing a spanning
#'   tree. Should be selected automatically in nearly all cases.
#' @param init_name a name for the initial plan, or \code{FALSE} to not include
#'   the initial plan in the output.  Defaults to the column name of the
#'   existing plan, or "\code{<init>}" if the initial plan is sampled.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended.
#' @param silent Whether to suppress all diagnostic information.
#'
#' @return \code{redist_mergesplit} returns an object of class
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
#' @examples \donttest{
#' data(fl25)
#'
#' fl_map = redist_map(fl25, ndists=3, pop_tol=0.1)
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
#' @order 1
#' @export
redist_mergesplit = function(map, nsims, warmup=floor(nsims/2),
                             init_plan=NULL, counties=NULL, compactness=1,
                             constraints=list(), constraint_fn=function(m) rep(0, ncol(m)),
                             adapt_k_thresh=0.975, k=NULL, init_name=NULL,
                             verbose=TRUE, silent=FALSE) {
    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)
    ndists = attr(map, "ndists")
    warmup = max(warmup, 0L)

    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")
    if (nsims < 1)
        stop("`nsims` must be positive.")

    exist_name = attr(map, "existing_col")
    counties = rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(init_plan) && !is.null(exist_name)) {
        init_plan = as.integer(as.factor(get_existing(map)))
        if (is.null(init_name)) init_name = exist_name
    }
    if (length(init_plan) == 0L || isTRUE(init_plan == "sample")) {
        init_plan = as.integer(get_plans_matrix(
            redist_smc(map, 1, counties, resample=FALSE, ref_name=FALSE, silent=TRUE)))
        if (is.null(init_name)) init_name = "<init>"
    }
    stopifnot(length(init_plan) == V)
    stopifnot(max(init_plan) == ndists)

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
    init_pop = pop_tally(matrix(init_plan, ncol=1), pop, ndists)
    if (any(init_pop < pop_bounds[1]) | any(init_pop > pop_bounds[3]))
        stop("Provided initialization does not meet population bounds.")
    if (any(pop >= get_target(map)))
        stop("Units ", which(pop >= get_target(map)),
             " have population larger than the district target.\n",
             "Redistricting impossible.")

    algout = ms_plans(nsims+1L, adj, init_plan, counties, pop, ndists,
                     pop_bounds[2], pop_bounds[1], pop_bounds[3], compactness,
                     constraints$status_quo$strength, constraints$status_quo$current, n_current,
                     constraints$vra$strength, constraints$vra$tgt_vra_min,
                     constraints$vra$tgt_vra_other, constraints$vra$pow_vra, proc$min_pop,
                     constraints$hinge$strength, constraints$hinge$tgts_min,
                     constraints$incumbency$strength, constraints$incumbency$incumbents,
                     constraints$splits$strength, constraints$multisplits$strength,
                     adapt_k_thresh, k, verbosity)

    plans <- algout$plans
    acceptances = as.logical(algout$mhdecisions)

    out = new_redist_plans(plans[, -1:-(warmup+1), drop=FALSE],
                           map, "mergesplit", NULL, FALSE,
                           compactness = compactness,
                           constraints = constraints,
                           adapt_k_thresh = adapt_k_thresh,
                           mh_acceptance = mean(acceptances))

    if (warmup == 0) {
        out <- out %>% mutate(mcmc_accept = rep(acceptances, each = ndists))
    } else {
        out <- out %>% mutate(mcmc_accept = rep(acceptances[-seq_len(warmup)], each = ndists))
    }


    if (!is.null(init_name) && !isFALSE(init_name)) {
        out = add_reference(out, init_plan, init_name)
    }

    out
}
