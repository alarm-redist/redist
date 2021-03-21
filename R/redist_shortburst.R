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
#' @param score_fn A function which takes
#' @param stop_at a threshold to stop optimization at.
#' @param burst_size The size of each burst. 10 is recommended.
#' @param max_bursts The maximum number of bursts to run before returning.
#' @param maximize if \code{TRUE}, try to maximize the score; otherwise, try to
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
#' higher values preferring more compact districts. Must be nonnegative. See
#' \code{\link{redist_mergsplit}} for more information.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
#' the algorithm does not appear to be sampling from the target distribution.
#' Must be between 0 and 1.
#' @param return_all whether to return all the
#'   Recommended for monitoring purposes.
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
#' @examples \dontrun{
#' data(iowa)
#'
#' iowa_map = redist_map(iowa, existing_plan=cd, pop_tol=0.01)
#' redist_shortburst(iowa_map, scorer_frac_kept(iowa_map), max_bursts=50)
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_shortburst = function(map, score_fn=NULL, stop_at=NULL, burst_size=10L,
                             max_bursts=500L, maximize=TRUE, init_plan=NULL,
                             counties=NULL, compactness=1, adapt_k_thresh=0.975,
                             return_all=TRUE, verbose=TRUE) {
    map = validate_redist_map(map)
    V = nrow(map)
    adj = get_adj(map)
    ndists = attr(map, "ndists")

    burst_size = as.integer(burst_size)
    max_bursts = as.integer(max_bursts)

    stopifnot(is.function(score_fn))
    if (!is.numeric(stop_at)) {
        stop_at = if (maximize) Inf else -Inf
    }

    if (compactness < 0) stop("Compactness parameter must be non-negative")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        stop("`adapt_k_thresh` parameter must lie in [0, 1].")
    if (burst_size < 1 || max_bursts < 1)
        stop("`burst_size` and `max_bursts` must be positive.")

    if (is.null(init_plan)) init_plan = as.integer(as.factor(get_existing(map)))
    if (length(init_plan) == 0L || isTRUE(init_plan == "sample")) {
        init_plan = as.integer(get_plan_matrix(
            redist_smc(map, 1, counties, resample=FALSE, silent=TRUE)))
    }
    stopifnot(length(init_plan) == V)
    stopifnot(max(init_plan) == ndists)

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

    pop_bounds = attr(map, "pop_bounds")
    pop = map[[attr(map, "pop_col")]]

    # kind of hacky -- extract k=... from outupt
    out = capture.output({
        x <- ms_plans(1, adj, init_plan, counties, pop, ndists, pop_bounds[2],
                      pop_bounds[1], pop_bounds[3], compactness,
                      0, rep(1, ndists), ndists, 0, 0, 0, 1, rep(0, V),
                      0, 0, 0, rep(1, ndists), adapt_k_thresh, 0L, verbosity=2)
    }, type="output")
    rm(x)
    k = as.integer(na.omit(stringr::str_match(out, "Using k = (\\d+)")[,2]))

    run_burst = function(init) {
        ms_plans(burst_size + 1L, adj, init, counties, pop, ndists,
                 pop_bounds[2], pop_bounds[1], pop_bounds[3], compactness,
                 0, rep(1, ndists), ndists, 0, 0, 0, 1, rep(0, V),
                 0, 0, 0, rep(1, ndists), 1.0, k, verbosity=0)[, -1L]
    }

    burst = 1
    out_mat = matrix(0L, nrow=V, ncol=max_bursts+1)
    out_mat[, 1] = init_plan

    scores = numeric(max_bursts+1)
    scores[1] = score_fn(out_mat[, 1, drop=FALSE])

    if (verbose) {
        cat("MERGE-SPLIT SHORT BURSTS\n")
        cat("Sampling up to", max_bursts, "bursts of", burst_size,
            "iterations each.\n")
        cat("Burst  Improve?  Score\n")
    }
    report_int = round(max_bursts / 10)
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
            if (verbose) cat(sprintf("% 5d     !      %f\n", burst, best_score))
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
                           converged = converged,
                           score_fn = deparse(substitute(score_fn)))
    out$score = rep(scores[out_idx], each=ndists)
    if (return_all) {
        out = add_reference(out, init_plan, "<init>")
        out$score[1:ndists] = scores[1]
    }

    out
}


#' Scoring fuctions for `redist_shortburst`
#'
#' The output of these functions may be passed into `redist_shortburst()` as
#' `score_fn`.
#'
#' Function details:
#'
#' - `scorer_frac_kept` returns the fraction of edges kept in each district.
#' Higher values mean more compactness.
#' - `scorer_group_pct` returns the `k`-th top group percentage across districts.
#' For example, if the group is Democratic voters and `k=3`, then the function
#' returns the 3rd-highest fraction of Democratic voters across all districts.
#' Can be used to target `k` VRA districts or partisan gerrymanders.
#'
#' @param map A \code{\link{redist_map}} object.
#'
#' @return A single numeric value, wherel larger values are better.
#'
#' @name scorers
#' @md
NULL

#' @rdname scorers
#'
#' @export
scorer_frac_kept = function(map) {
    adj = get_adj(map)
    edges = sum(sapply(adj, length)) / 2
    ndists = attr(map, "ndists")

    function(plans) {
        (edges - n_removed(adj, plans, ndists)) / edges
    }
}

#' @rdname scorers
#'
#' @param group_pop A numeric vector with the population of the group for every precinct.
#' @param total_pop A numeric vector with the population for every precinct.
#'
#' @export
scorer_group_pct = function(map, group_pop, total_pop, k=1) {
    group_pop = rlang::eval_tidy(rlang::enquo(group_pop), map)
    total_pop = rlang::eval_tidy(rlang::enquo(total_pop), map)
    ndists = attr(map, "ndists")

    if (k == 1) {
        function(plans) {
            colmax(group_pct(plans, group_pop, total_pop, ndists))
        }
    } else if (k == ndists) {
        function(plans) {
            colmin(group_pct(plans, group_pop, total_pop, ndists))
        }
    } else {
        function(plans) {
            group_pct_top_k(plans, group_pop, total_pop, k, ndists)
        }
    }
}
