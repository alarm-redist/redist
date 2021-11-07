#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/03/20
# Purpose: short-burst redistricting optimization
####################################################

#' Redistricting Optimization through Short Bursts
#'
#' This function uses [redist_mergesplit()] to optimize a redistrict plan
#' according to a user-provided criteria. It does so by running the Markov chain
#' for "short bursts" of usually 10 iterations, and then starting the chain anew
#' from the best plan in the burst, according to the criteria. This implements
#' the ideas in the below-referenced paper, "Voting Rights, Markov Chains, and
#' Optimization by Short Bursts."
#'
#' @param map A \code{\link{redist_map}} object.
#' @param score_fn A function which takes a matrix of plans and returns a score
#' for each plan. Can also be a purrr-style anonymous function. See
#' [`?scorers`][scorers] for some function factories for common scoring rules.
#' @param stop_at A threshold to stop optimization at.
#' @param burst_size The size of each burst. 10 is recommended for mergesplit and 50 for flip.
#' @param max_bursts The maximum number of bursts to run before returning.
#' @param maximize If \code{TRUE}, try to maximize the score; otherwise, try to
#' minimize it.
#' @param init_plan The initial state of the map. If not provided, will default to
#'   the reference map of the \code{map} object, or if none exists, will sample
#'   a random initial state using \code{\link{redist_smc}}. You can also request
#'   a random initial state by setting \code{init_plan="sample"}.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided, the
#' algorithm will only generate maps which split up to \code{ndists-1} counties.
#' If no county-split constraint is desired, this parameter should be left blank.
#' @param compactness Controls the compactness of the generated districts, with
#' higher values preferring more compact districts. Must be non-negative. See
#' \code{\link{redist_mergesplit}} for more information.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
#' the algorithm does not appear to be sampling from the target distribution.
#' Must be between 0 and 1.
#' @param return_all Whether to return all the
#'   Recommended for monitoring purposes.
#' @param backend the MCMC algorithm to use within each burst, either
#'   "mergesplit" or "flip".
#' @param flip_lambda The parameter determining the number of swaps to attempt each iteration of flip mcmc.
#' The number of swaps each iteration is equal to Pois(lambda) + 1. The default is 0.
#' @param flip_eprob  The probability of keeping an edge connected in flip mcmc. The default is 0.05.
#' @param flip_constraints A list of constraints to use for flip mcmc. Can be created with
#' \code{flip_constraints_helper}. Defaults to an edges-removed compactness constraint with weight 0.6.
#' @param verbose Whether to print out intermediate information while sampling.
#'   Recommended for monitoring purposes.
#'
#' @return a \code{redist_plans} object containing the final best plan
#' (or the best plans after each burst, if \code{return_all=TRUE}.
#'
#' @references
#' Cannon, S., Goldbloom-Helzner, A., Gupta, V., Matthews, J. N., & Suwal, B.
#' (2020). Voting Rights, Markov Chains, and Optimization by Short Bursts. arXiv
#' preprint arXiv:2011.02288.
#'
#' @examples \donttest{
#' data(iowa)
#'
#' iowa_map = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.01)
#' redist_shortburst(iowa_map, scorer_frac_kept(iowa_map), max_bursts=50)
#' redist_shortburst(iowa_map, ~ 1 - scorer_frac_kept(iowa_map)(.), max_bursts=50)
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_shortburst = function(map, score_fn=NULL, stop_at=NULL,
                             burst_size = ifelse(backend == 'mergesplit', 10L, 50L),
                             max_bursts=500L, maximize=TRUE, init_plan=NULL,
                             counties=NULL, compactness=1, adapt_k_thresh=0.975,
                             return_all=TRUE, backend="mergesplit",
                             flip_lambda = 0, flip_eprob = 0.05, flip_constraints = list(),
                             verbose=TRUE) {

    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)
    ndists = attr(map, "ndists")

    burst_size = as.integer(burst_size)
    max_bursts = as.integer(max_bursts)
    match.arg(backend, c("flip", "mergesplit"))

    score_fn = rlang::as_closure(score_fn)
    stopifnot(is.function(score_fn))
    if (!is.numeric(stop_at)) {
        stop_at = if (maximize) Inf else -Inf
    }

    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")

    if (burst_size < 1 || max_bursts < 1)
        stop("`burst_size` and `max_bursts` must be positive.")

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

    if (is.null(init_plan)) init_plan = as.integer(as.factor(get_existing(map)))
    if (length(init_plan) == 0L || isTRUE(init_plan == "sample")) {
        init_plan = as.integer(get_plans_matrix(
            redist_smc(map, 1, counties, resample=FALSE, ref_name=FALSE, silent=TRUE)))
    }
    stopifnot(length(init_plan) == V)
    stopifnot(max(init_plan) == ndists)


    if (backend == 'mergesplit') {
        pop_bounds = attr(map, "pop_bounds")
    } else {
        pop_tol <- get_pop_tol(map)
    }

    pop = map[[attr(map, "pop_col")]]
    if (any(pop >= get_target(map)))
        stop("Units ", which(pop >= get_target(map)),
             " have population larger than the district target.\n",
             "Redistricting impossible.")


    if (backend == "mergesplit") {
        # kind of hacky -- extract k=... from outupt
        if (!requireNamespace("utils", quietly=TRUE)) stop()
        out = utils::capture.output({
            x <- ms_plans(1, adj, init_plan, counties, pop, ndists, pop_bounds[2],
                          pop_bounds[1], pop_bounds[3], compactness,
                          0, rep(1, ndists), ndists, 0, 0, 0, 1, rep(0, V),
                          0, 0, 0, rep(1, ndists), 0, 0, adapt_k_thresh, 0L, verbosity=2)
        }, type="output")
        rm(x)
        k = as.integer(stats::na.omit(stringr::str_match(out, "Using k = (\\d+)")[,2]))

        run_burst = function(init) {
            ms_plans(burst_size + 1L, adj, init, counties, pop, ndists,
                     pop_bounds[2], pop_bounds[1], pop_bounds[3], compactness,
                     0, rep(1, ndists), ndists, 0, 0, 0, 1, rep(0, V),
                     0, 0, 0, rep(1, ndists), 0, 0, 1.0, k, verbosity=0)$plans[, -1L]
        }
    } else {
        flip_constraints <- process_flip_constr(constraints = flip_constraints,
                                                nrow(map))

        if (flip_eprob <= 0 || flip_eprob >= 1) {
            stop("flip_eprob must be in the interval (0, 1).")
        }
        if (flip_lambda < 0) {
            stop("flip_lambda must be a nonnegative integer.")
        }

        if (all(flip_constraints$similarity$plan == 1)) {
          if (min(init_plan) == 1) {
            flip_constraints$similarity$plan <- init_plan - 1
          } else {
            flip_constraints$similarity$plan <- init_plan
          }
        }

        run_burst <- function(init) {
            skinny_flips(adj = adj, init_plan = init, total_pop = pop,
                        pop_tol = pop_tol, nsims = burst_size,
                        eprob = flip_eprob, lambda = flip_lambda,
                        constraints = flip_constraints)
        }
    }


    burst = 1
    out_mat = matrix(0L, nrow=V, ncol=max_bursts+1)
    out_mat[, 1] = init_plan

    scores = numeric(max_bursts+1)
    scores[1] = score_fn(out_mat[, 1, drop=FALSE])

    if (verbose) {
        if (backend == 'mergesplit') {
            cat("MERGE-SPLIT SHORT BURSTS\n")
        } else {
            cat('FLIP SHORT BURSTS\n')
        }
        cat("Sampling up to", max_bursts, "bursts of", burst_size,
            "iterations each.\n")
        cat("Burst  Improve?  Score\n")
    }
    report_int = max(round(max_bursts / 10), 1)
    improve_ch = sample(c("\U0001F973", "\U0001F600", "\U0001F60E",
                           "\U0001F642", "\U0001F386", "\U0001F387",
                           "\U0001F942", "\U0001F383", "\U0001FA85",
                           "\U0001F4A5", "\U0001F389", "\U26C4",
                          "\U0001F31F", "\U0001F308"))
    improve_ct = 1L
    for (burst in 1:max_bursts) {
        plans = run_burst(out_mat[, burst])
        plan_scores = score_fn(plans)

        prev_score = scores[burst]
        if (maximize) {
            best_idx = which.max(plan_scores)
            best_score = plan_scores[best_idx]
            condition = best_score > prev_score
        } else {
            best_idx = which.min(plan_scores)
            best_score = plan_scores[best_idx]
            condition = best_score < prev_score
        }

        if (condition) {
            out_mat[, burst+1L] = plans[, best_idx]
            scores[burst+1L] = best_score
            if (verbose) {
                improve_ct = (improve_ct %% length(improve_ch)) + 1L
                cat(sprintf("% 5d     %s     %f\n", burst,
                            improve_ch[improve_ct], best_score))
            }
        } else {
            out_mat[, burst+1L] = out_mat[, burst]
            scores[burst+1L] = prev_score
            if (verbose && burst %% report_int == 0)
                cat(sprintf("% 5d            %f\n", burst, prev_score))
        }


        if (maximize && scores[burst+1L] >= stop_at) break
        if (!maximize && scores[burst+1L] <= stop_at) break
    }

    out_idx = if (return_all) 2:(burst+1L) else burst+1L
    if (maximize)
        converged = scores[burst+1L] >= stop_at
    else
        converged = scores[burst+1L] <= stop_at

    out = new_redist_plans(out_mat[, out_idx, drop=FALSE], map, "shortburst",
                           wgt=NULL, resampled=FALSE,
                           burst_size = burst_size,
                           n_bursts = burst,
                           backend = backend,
                           converged = converged,
                           score_fn = deparse(substitute(score_fn)))
    out$score = rep(scores[out_idx], each=ndists)

    if (return_all) {
        out = add_reference(out, init_plan, "<init>")
        out$score[1:ndists] = scores[1]
    }

    out
}


#' Scoring functions for `redist_shortburst`
#'
#' The output of these functions may be passed into `redist_shortburst()` as
#' `score_fn`.  Scoring functions have type `redist_scorer` and may be combined
#' together using basic arithmetic operations.
#'
#' Function details:
#'
#' - `scorer_group_pct` returns the `k`-th top group percentage across districts.
#' For example, if the group is Democratic voters and `k=3`, then the function
#' returns the 3rd-highest fraction of Democratic voters across all districts.
#' Can be used to target `k` VRA districts or partisan gerrymanders.
#' - `scorer_pop_dev` returns the maximum population deviation within a plan.
#' Smaller values are closer to population parity, so use `maximize=FALSE` with
#' this scorer.
#' - `scorer_splits` returns the fraction of counties that are split within a
#' plan. Higher values have more county splits, so use `maximize=FALSE` with
#' this scorer.
#' - `scorer_frac_kept` returns the fraction of edges kept in each district.
#' Higher values mean more compactness.
#' - `scorer_polsby_popper` returns the `m`-th Polsby Popper score within a plan.
#' Higher scores correspond to more compact districts.  Use `m=ndists/2` to
#' target the median compactness, `m=1` to target the minimum compactness.
#' - `scorer_status_quo` returns 1 - the rescaled variation of information
#' distance between the plan and the `existing_plan`. Larger values indicate the
#' plan is closer to the existing plan.
#'
#' @param map A \code{\link{redist_map}} object.
#'
#' @return A scoring function of class `redist_scorer`. single numeric value, where larger values are better for `frac_kept`,
#' `group_pct`, and `polsby_popper` and smaller values are better for `splits` and `pop_dev`.
#'
#' @examples
#' \donttest{
#' data(iowa)
#' iowa_map = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05, total_pop = pop)
#'
#' scorer_frac_kept(iowa_map)
#' scorer_status_quo(iowa_map)
#' scorer_group_pct(iowa_map, dem_08, tot_08, k=2)
#' 1.5*scorer_frac_kept(iowa_map) + 0.4*scorer_status_quo(iowa_map)
#' 1.5*scorer_frac_kept(iowa_map) + scorer_frac_kept(iowa_map)*scorer_status_quo(iowa_map)
#' }
#'
#' @concept prepare
#' @name scorers
#' @md
NULL

#' @rdname scorers
#' @order 5
#'
#' @export
scorer_frac_kept = function(map) {
    adj = get_adj(map)
    edges = sum(sapply(adj, length)) / 2
    ndists = attr(map, "ndists")

    fn = function(plans) {
        (edges - n_removed(adj, plans, ndists)) / edges
    }
    class(fn) = c("redist_scorer", "function")
    fn
}

#' @rdname scorers
#' @order 1
#'
#' @param group_pop A numeric vector with the population of the group for every precinct.
#' @param total_pop A numeric vector with the population for every precinct.
#' @param k the k-th from the top group fraction to return as the score.
#'
#' @export
scorer_group_pct = function(map, group_pop, total_pop, k=1) {
    group_pop = eval_tidy(enquo(group_pop), map)
    total_pop = eval_tidy(enquo(total_pop), map)
    ndists = attr(map, "ndists")

    if (k == 1) {
        fn = function(plans) {
            colmax(group_pct(plans, group_pop, total_pop, ndists))
        }
    } else if (k == ndists) {
        fn = function(plans) {
            colmin(group_pct(plans, group_pop, total_pop, ndists))
        }
    } else {
        fn = function(plans) {
            group_pct_top_k(plans, group_pop, total_pop, k, ndists)
        }
    }
    class(fn) = c("redist_scorer", "function")
    fn
}

#' @rdname scorers
#' @order 2
#'
#' @export
scorer_pop_dev <- function(map) {
  ndists <- attr(map, 'ndists')
  total_pop = map[[attr(map, "pop_col")]]
  stopifnot(!is.null(total_pop))

  fn = function(plans) {
    max_dev(plans, total_pop, ndists)
  }
  class(fn) <- c("redist_scorer", "function")
  fn
}

#' @rdname scorers
#' @order 3
#'
#' @param counties A numeric vector with an integer from 1:n_counties
#'
#' @export
scorer_splits <- function(map, counties) {
  counties <- eval_tidy(enquo(counties), map)
  counties <- as.integer(as.factor(counties))

  fn = function(plans) {
    nd <- length(unique(plans[, 1]))
    splits(plans, counties, nd, 1)/length(unique(counties))
  }
  class(fn) <- c("redist_scorer", "function")
  fn
}

#' @rdname scorers
#' @order 4
#'
#' @param counties A numeric vector with an integer from 1:n_counties
#'
#' @export
scorer_multisplits <- function(map, counties) {
  counties <- eval_tidy(enquo(counties), map)
  counties <- as.integer(as.factor(counties))

  fn = function(plans) {
    splits(plans, counties, attr(map, 'ndists'), 2)/length(unique(counties))
  }
  class(fn) <- c("redist_scorer", "function")
  fn
}

#' @rdname scorers
#' @order 6
#'
#' @param perim_df perimeter distance dataframe from \code{\link{redist.prep.polsbypopper}}
#' @param areas area of each precinct (ie \code{st_area(map)})
#' @param m the m-th from the bottom Polsby Popper to return as the score. Defaults to 1,
#' the minimum Polsby Popper score
#'
#' @export
scorer_polsby_popper <- function(map, perim_df=NULL, areas=NULL, m = 1) {
  ndists <- attr(map, 'ndists')
  if (is.null(perim_df)) perim_df = redist.prep.polsbypopper(map)
  if (is.null(areas)) areas = sf::st_area(sf::st_geometry(map))

  fn = function(plans) {
    pp <- polsbypopper(
      from = perim_df$origin, to = perim_df$touching, area = areas,
      perimeter = perim_df$edge, dm = plans, nd = ndists
    )

    k_smallest(x = pp, k = m)
  }
  class(fn) <- c("redist_scorer", "function")
  fn
}


#' @rdname scorers
#' @order 7
#'
#' @param existing_plan A vector containing the current plan.
#'
#' @export
scorer_status_quo = function(map, existing_plan=get_existing(map)) {
    exsiting_plan = eval_tidy(enquo(existing_plan), map)
    pop = map[[attr(map, "pop_col")]]
    ndists = attr(map, "ndists")

    stopifnot(!is.null(existing_plan))
    stopifnot(!is.null(pop))
    stopifnot(ndists == length(unique(existing_plan)))

    fn = function(plans) {
        1 - 0.5*var_info_vec(plans, existing_plan, pop) / log(ndists)
    }
    class(fn) = c("redist_scorer", "function")
    fn
}

#' Scoring function arithmetic
#'
#' `redist_scorer` functions may be multiplied by constants and/or added
#' together to form linear combinations.
#'
#' @name scorer-arith
#' @concept prepare
#' @md
#' @returns function of class redist_scorer
NULL

#' @rdname scorer-arith
#'
#' @param x a numeric or a `redist_scorer` function, from [`scorers`]
#' @param fn2 a `redist_scorer` function, from [`scorers`]
#'
#' @export
`*.redist_scorer` = function(x, fn2) {
    stopifnot(is.numeric(x) || inherits(x, "redist_scorer"))
    stopifnot(inherits(fn2, "redist_scorer"))

    if (is.numeric(x)) {
      rlang::fn_body(fn2) = rlang::expr({!!x * !!rlang::fn_body(fn2)})
      return(fn2)
    } else {
      fn = function(plans) { x(plans) * fn2(plans) }
      class(fn) = c("redist_scorer", "function")
    }

    fn
}

#' @rdname scorer-arith
#'
#' @param fn1 a `redist_scorer` function, from [`scorers`]
#' @param fn2 a `redist_scorer` function, from [`scorers`]
#'
#' @export
`+.redist_scorer` = function(fn1, fn2) {
    stopifnot(inherits(fn1, "redist_scorer"))
    stopifnot(inherits(fn2, "redist_scorer"))

    fn = function(plans) { fn1(plans) + fn2(plans) }
    class(fn) = c("redist_scorer", "function")
    fn
}

#' @rdname scorer-arith
#'
#' @param fn1 a `redist_scorer` function, from [`scorers`]
#' @param fn2 a `redist_scorer` function, from [`scorers`]
#'
#' @export
`-.redist_scorer` = function(fn1, fn2) {
    stopifnot(inherits(fn1, "redist_scorer"))
    stopifnot(inherits(fn2, "redist_scorer"))

    fn = function(plans) { fn1(plans) - fn2(plans) }
    class(fn) = c("redist_scorer", "function")
    fn
}
