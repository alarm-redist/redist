#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/03/20
# Purpose: short-burst redistricting optimization
####################################################

#' Redistricting Optimization through Short Bursts
#'
#' This function uses [redist_mergesplit()] or [redist_flip()] to optimize a
#' redistrict plan according to a user-provided criteria. It does so by running
#' the Markov chain for "short bursts" of usually 10 iterations, and then
#' starting the chain anew from the best plan in the burst, according to the
#' criteria. This implements the ideas in the below-referenced paper, "Voting
#' Rights, Markov Chains, and Optimization by Short Bursts."
#'
#' @param map A [redist_map] object.
#' @param score_fn A function which takes a matrix of plans and returns a score
#'   (or, generally, a row vector) for each plan. Can also be a purrr-style
#'   anonymous function. See [`?scorers`][scorers] for some function factories
#'   for common scoring rules.
#' @param stop_at A threshold to stop optimization at. When `score_fn` returns a
#'   row vector per plan, `maximize` can be an equal-length vector specifying a
#'   threshold for each dimension, which must all be met for the algorithm to
#'   stop.
#' @param burst_size The size of each burst. 10 is recommended for the
#'   `mergesplit` backend and 50 for the `flip` backend. Can also provide
#'   burst schedule function which takes the current iteration (an integer)
#'   and returns the desired burst size. This can be a random function.
#' @param max_bursts The maximum number of bursts to run before returning.
#' @param maximize If `TRUE`, try to maximize the score; otherwise, try to
#'   minimize it. When `score_fn` returns a row vector per plan, `maximize` can
#'   be an equal-length vector specifying whether each dimension should be
#'   maximized or minimized.
#' @param init_plan The initial state of the map. If not provided, will default to
#' the reference map of the `map` object, or if none exists, will sample
#' a random initial state using [redist_smc()]. You can also request
#' a random initial state by setting `init_plan="sample"`.
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided, the
#' algorithm will only generate maps which split up to `ndists-1` counties.
#' If no county-split constraint is desired, this parameter should be left blank.
#' @param constraints A `redist_constr` with Gibbs constraints.
#' @param compactness Controls the compactness of the generated districts, with
#' higher values preferring more compact districts. Must be non-negative. See
#' \code{\link{redist_mergesplit}} for more information.
#' @param adapt_k_thresh The threshold value used in the heuristic to select a
#' value `k_i` for each splitting iteration.
#' @param reversible If `FALSE` and `backend="mergesplit"`, the Markov chain
#' used will not be reversible. This may speed up optimization.
#' @param fixed_k If not `NULL`, will be used to set the `k` parameter for the
#'   `mergesplit` backend. If e.g. `k=1` then the best edge in each spanning
#'   tree will be used.  Lower values may speed up optimization at the
#'   cost of the Markov chain no longer targeting a known distribution.
#'   Recommended only in conjunction with `reversible=FALSE`.
#' @param return_all Whether to return all the burst results or just the best
#' one (generally, the Pareto frontier). Recommended for monitoring purposes.
#' @param thin Save every `thin`-th sample. Defaults to no thinning (1). Ignored
#' if `return_all=TRUE`.
#' @param backend the MCMC algorithm to use within each burst, either
#' "mergesplit" or "flip".
#' @param flip_lambda The parameter determining the number of swaps to attempt each iteration of flip mcmc.
#' The number of swaps each iteration is equal to Pois(lambda) + 1. The default is 0.
#' @param flip_eprob  The probability of keeping an edge connected in flip mcmc. The default is 0.05.
#' @param verbose Whether to print out intermediate information while sampling.
#' Recommended for monitoring purposes.
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
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
#' redist_shortburst(iowa_map, scorer_frac_kept(iowa_map), max_bursts = 50)
#' redist_shortburst(iowa_map, ~ 1 - scorer_frac_kept(iowa_map)(.), max_bursts = 50)
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_shortburst <- function(map, score_fn = NULL, stop_at = NULL,
                              burst_size = ifelse(backend == "mergesplit", 10L, 50L),
                              max_bursts = 500L, maximize = TRUE, init_plan = NULL,
                              counties = NULL,  constraints = redist_constr(map),
                              compactness = 1, adapt_k_thresh = 0.95,
                              reversible=TRUE, fixed_k=NULL,
                              return_all = TRUE, thin = 1L, backend = "mergesplit",
                              flip_lambda = 0, flip_eprob = 0.05,
                              verbose = TRUE) {

    map <- validate_redist_map(map)
    V <- nrow(map)
    adj <- get_adj(map)
    ndists <- attr(map, "ndists")
    thin <- as.integer(thin)

    if (!is.function(burst_size)) {
        per_burst <- as.integer(burst_size)
        burst_size <- function(i) per_burst
    }
    max_bursts <- as.integer(max_bursts)
    match.arg(backend, c("flip", "mergesplit"))

    score_fn <- rlang::as_closure(score_fn)
    stopifnot(is.function(score_fn))

    if (compactness < 0)
        cli_abort("{.arg compactness} must be non-negative.")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")

    if (burst_size(1) < 1 || max_bursts < 1)
        cli_abort("{.arg burst_size} and {.arg max_bursts} must be positive.")
    if (thin < 1 || thin > max_bursts)
        cli_abort("{.arg thin} must be a positive integer, and no larger than {.arg max_bursts}.")

    counties <- rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(counties)) {
        counties <- rep(1, V)
    } else {
        if (any(is.na(counties)))
            cli_abort("County vector must not contain missing values.")

        # handle discontinuous counties
        component <- contiguity(adj, vctrs::vec_group_id(counties))
        counties <- dplyr::if_else(component > 1,
            paste0(as.character(counties), "-", component),
            as.character(counties)) %>%
            as.factor() %>%
            as.integer()

        if (any(component > 1)) {
            cli_warn("Counties were not contiguous; expect additional splits.")
        }
    }

    if (is.null(init_plan)) init_plan <- vctrs::vec_group_id(get_existing(map))
    if (length(init_plan) == 0L || isTRUE(init_plan == "sample")) {
        init_plan <- as.integer(get_plans_matrix(
            redist_smc(map, 10, counties, resample = FALSE, ref_name = FALSE, silent = TRUE))[, 1])
    }

    # check init
    if (length(init_plan) != V)
        cli_abort("{.arg init_plan} must be as long as the number of units as `map`.")
    if (max(init_plan) != ndists)
        cli_abort("{.arg init_plan} must have the same number of districts as `map`.")
    if (any(contiguity(adj, init_plan) != 1))
        cli_warn("{.arg init_plan} should have contiguous districts.")


    if (backend == "mergesplit") {
        pop_bounds <- attr(map, "pop_bounds")
    } else {
        pop_tol <- get_pop_tol(map)
    }

    pop <- map[[attr(map, "pop_col")]]
    if (any(pop >= get_target(map))) {
        too_big <- as.character(which(pop >= pop_bounds[3]))
        cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                    population larger than the district target.",
            "x" = "Redistricting impossible."))
    }

    if (!inherits(constraints, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    constraints <- as.list(constraints)

    if (backend == "mergesplit") {
        control = list(adapt_k_thresh=adapt_k_thresh, do_mh=reversible)
        if (is.null(fixed_k)) {
            x <- ms_plans(1, adj, init_plan, counties, pop, ndists, pop_bounds[2],
                          pop_bounds[1], pop_bounds[3], compactness,
                          list(), control, 0L, 1L, verbosity = 0)
            k <- x$est_k
        } else {
            k = fixed_k
        }

        run_burst <- function(init, steps) {
            ms_plans(steps, adj, init, counties, pop, ndists,
                pop_bounds[2], pop_bounds[1], pop_bounds[3], compactness,
                constraints, control, k, 1L, verbosity = 0)$plans[, -1L]
        }
    } else {

        if (flip_eprob <= 0 || flip_eprob >= 1) {
            cli_abort("{.arg flip_eprob} must be in the interval (0, 1).")
        }
        if (flip_lambda < 0) {
            cli_abort("{.arg flip_lambda} must be a nonnegative integer.")
        }

        run_burst <- function(init, steps) {
            skinny_flips(adj = adj, init_plan = init, total_pop = pop,
                pop_tol = pop_tol, nsims = steps,
                eprob = flip_eprob, lambda = flip_lambda,
                constraints = constraints)
        }
    }


    burst <- 1
    n_out <- max_bursts %/% thin
    out_mat <- matrix(0L, nrow = V, ncol = n_out)
    burst_sizes <- integer(n_out)
    cur_best <- matrix(init_plan, ncol=1)
    rescale <- 1 - maximize * 2

    cur_best_scores <- score_fn(matrix(init_plan, ncol = 1))
    score_init = cur_best_scores
    if (!is.numeric(stop_at)) {
        stop_at <- -Inf
    } else {
        stop_at = rescale * stop_at
    }
    if (!is.matrix(cur_best_scores)) {
        cur_best_scores = matrix(cur_best_scores, ncol=1)
        rownames(cur_best_scores) = "score"
    } else {
        cur_best_scores = t(cur_best_scores)
        if (!is.null(names(rescale))) {
            rescale = rescale[match(rownames(cur_best_scores), names(rescale))]
        }
    }
    cur_best_scores <- cur_best_scores * rescale
    dim_score <- nrow(cur_best_scores)

    scores <- matrix(nrow=n_out, ncol=dim_score)
    colnames(scores) = rownames(cur_best_scores)

    if (verbose) {
        fmt_score <- function(x) {
            paste0(sprintf("%f", x * rescale), collapse=" ")
        }
        if (backend == "mergesplit") {
            cat("MERGE-SPLIT SHORT BURSTS\n")
        } else {
            cat("FLIP SHORT BURSTS\n")
        }
        cat("Sampling up to", max_bursts, "bursts of", burst_size(1),
            "iterations each.\n")
        cat("Burst  Improve? ")
        cur_fmt_len <- nchar(fmt_score(cur_best_scores[, 1]))
        cols <- stringr::str_pad(colnames(scores), round(cur_fmt_len / dim_score))
        cat(cols, "\n")
    }
    report_int <- max(round(max_bursts/10), 1)
    improve_ch <- sample(c("\U0001F973", "\U0001F600", "\U0001F60E",
        "\U0001F642", "\U0001F386", "\U0001F387",
        "\U0001F942", "\U0001F383", "\U0001FA85",
        "\U0001F4A5", "\U0001F389", "\U26C4",
        "\U0001F31F", "\U0001F308"))
    improve_ct <- 1L
    idx <- 1L
    converged <- FALSE
    for (burst in 1:max_bursts) {
        this_burst_size <- burst_size(burst)
        burst_sizes[burst] <- this_burst_size
        if (this_burst_size <= 0) {
            cli_abort(c("Burst size must be at least 1.",
                        "x"="Found {this_burst_size} on iteration {burst}."))
        }
        keep <- seq_len(this_burst_size)
        burst_init = cur_best[, sample.int(ncol(cur_best), 1)]
        plans <- run_burst(burst_init, this_burst_size)[, keep, drop=FALSE]
        plan_scores <- t(matrix(score_fn(plans), ncol=dim_score))
        plan_scores <- plan_scores * rescale

        cur_best <- cbind(cur_best, plans)
        cur_best_scores <- cbind(cur_best_scores, plan_scores)

        dominated <- pareto_dominated(cur_best_scores)
        improved <- any(!tail(dominated, this_burst_size))
        # remove dominated plans
        cur_best <- cur_best[, !dominated, drop=FALSE]
        cur_best_scores <- cur_best_scores[, !dominated, drop=FALSE]

        # add new undominated plans
        out_idx = sample.int(ncol(cur_best), 1) # random plan from frontier
        if (improved) { # improvement
            if (verbose) {
                improve_ct <- (improve_ct %% length(improve_ch)) + 1L

                cat(sprintf("% 5d     %s     %s\n", burst,
                    improve_ch[improve_ct],
                    fmt_score(cur_best_scores[, out_idx])))
            }
        } else if (verbose && burst %% report_int == 0) {
            cat(sprintf("% 5d            %s\n", burst,
                        fmt_score(cur_best_scores[, out_idx])))
        }

        if (burst %% thin == 0) {
            idx <- burst %/% thin
            out_mat[, idx] <- cur_best[, out_idx]
            scores[idx, ] <- cur_best_scores[, out_idx] * rescale

            if (any(colSums(cur_best_scores <= stop_at) == dim_score)) {
                converged = TRUE
                break
            }
        }
    }

    if (return_all) {
        out_idx <- seq_len(idx)
        storage.mode(out_mat) <- "integer"

        pareto_scores = t(cur_best_scores * rescale)
        pareto_scores = pareto_scores[order(pareto_scores[, 1]), , drop=FALSE]

        out <- new_redist_plans(out_mat[, out_idx, drop = FALSE], map, "shortburst",
            wgt = NULL, resampled = FALSE,
            n_bursts = burst,
            backend = backend,
            converged = converged,
            pareto_front = cur_best,
            pareto_scores = pareto_scores,
            version = packageVersion("redist"),
            score_fn = deparse(substitute(score_fn)))
        score_mat = matrix(rep(scores[out_idx, ], each = ndists), ncol = dim_score)
        colnames(score_mat) = colnames(scores)
        out <- dplyr::mutate(out, as.data.frame(score_mat))
        out$burst_size = rep(burst_sizes[out_idx], each = ndists)

        out <- add_reference(out, init_plan, "<init>")
        idx_cols = ncol(out) - dim_score:1
        out[1:ndists, idx_cols] <- matrix(rep(score_init, each = ndists), ncol = dim_score)
    } else {
        out <- new_redist_plans(cur_best, map, "shortburst",
                                wgt = NULL, resampled = FALSE,
                                n_bursts = burst,
                                backend = backend,
                                converged = converged,
                                version = packageVersion("redist"),
                                score_fn = deparse(substitute(score_fn)))
        score_mat = matrix(rep(t(cur_best_scores * rescale), each = ndists),
                           ncol = dim_score)
        colnames(score_mat) = colnames(scores)
        out <- dplyr::mutate(out, as.data.frame(score_mat))
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
#' @return A scoring function of class `redist_scorer` which returns a single numeric value per plan.
#' Larger values are generally better for `frac_kept`, `group_pct`, and `polsby_popper`
#' and smaller values are better for `splits` and `pop_dev`.
#'
#' @examples
#' \donttest{
#' data(iowa)
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)
#'
#' scorer_frac_kept(iowa_map)
#' scorer_status_quo(iowa_map)
#' scorer_group_pct(iowa_map, dem_08, tot_08, k = 2)
#' 1.5*scorer_frac_kept(iowa_map) + 0.4*scorer_status_quo(iowa_map)
#' 1.5*scorer_frac_kept(iowa_map) + scorer_frac_kept(iowa_map)*scorer_status_quo(iowa_map)
#' cbind(
#'     comp = scorer_frac_kept(iowa_map),
#'     sq = scorer_status_quo(iowa_map)
#' )
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
scorer_frac_kept <- function(map) {
    adj <- get_adj(map)
    edges <- sum(sapply(adj, length))/2
    ndists <- attr(map, "ndists")

    fn <- function(plans) {
        (edges - n_removed(adj, plans, ndists))/edges
    }
    class(fn) <- c("redist_scorer", "function")
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
scorer_group_pct <- function(map, group_pop, total_pop, k = 1) {
    group_pop <- eval_tidy(enquo(group_pop), map)
    total_pop <- eval_tidy(enquo(total_pop), map)
    ndists <- attr(map, "ndists")

    if (k == 1) {
        fn <- function(plans) {
            colmax(group_pct(plans, group_pop, total_pop, ndists))
        }
    } else if (k == ndists) {
        fn <- function(plans) {
            colmin(group_pct(plans, group_pop, total_pop, ndists))
        }
    } else {
        fn <- function(plans) {
            group_pct_top_k(plans, group_pop, total_pop, k, ndists)
        }
    }
    class(fn) <- c("redist_scorer", "function")
    fn
}

#' @rdname scorers
#' @order 2
#'
#' @export
scorer_pop_dev <- function(map) {
    ndists <- attr(map, "ndists")
    total_pop <- map[[attr(map, "pop_col")]]
    stopifnot(!is.null(total_pop))

    fn <- function(plans) {
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
    counties <- vctrs::vec_group_id(counties)

    fn <- function(plans) {
        nd <- length(unique(plans[, 1]))
        splits(plans - 1, counties - 1, nd, 1)/length(unique(counties))
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
    counties <- vctrs::vec_group_id(counties)

    fn <- function(plans) {
        splits(plans, counties, attr(map, "ndists"), 2)/length(unique(counties))
    }
    class(fn) <- c("redist_scorer", "function")
    fn
}

#' @rdname scorers
#' @order 6
#'
#' @param perim_df perimeter distance dataframe from [prep_perims()]
#' @param areas area of each precinct (ie `st_area(map)`)
#' @param m the m-th from the bottom Polsby Popper to return as the score. Defaults to 1,
#' the minimum Polsby Popper score
#'
#' @export
scorer_polsby_popper <- function(map, perim_df = NULL, areas = NULL, m = 1) {
    ndists <- attr(map, "ndists")
    if (is.null(perim_df)) perim_df <- redistmetrics::prep_perims(map)
    if (is.null(areas)) areas <- sf::st_area(sf::st_geometry(map))

    fn <- function(plans) {
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
scorer_status_quo <- function(map, existing_plan = get_existing(map)) {
    exsiting_plan <- eval_tidy(enquo(existing_plan), map)
    pop <- map[[attr(map, "pop_col")]]
    ndists <- attr(map, "ndists")

    stopifnot(!is.null(existing_plan))
    stopifnot(!is.null(pop))
    stopifnot(ndists == length(unique(existing_plan)))

    fn <- function(plans) {
        1 - 0.5*var_info_vec(plans, existing_plan, pop)/log(ndists)
    }
    class(fn) <- c("redist_scorer", "function")
    fn
}



#' Combine scoring functions
#'
#' `redist_scorer` functions may be combined together to optimize along multiple
#' dimensions. Rather than linearly combining multiple scorers to form a single
#' objective as with [scorer-arith], these functions allow analysts to approximate
#' the Pareto frontier for a set of scorers.
#'
#' @name scorer-combine
#' @concept prepare
#' @md
#' @returns function of class redist_scorer. Will return a matrix with each
#'   column containing every plan's scores for a particular scoring function.
NULL

#' @rdname scorer-combine
#'
#' @param ... a numeric or a `redist_scorer` function, from [`scorers`]
#' @param deparse.level As in [cbind()].
#'
#' @export
combine_scorers <- function(...) {
    cbind(...)
}

#' @rdname scorer-combine
#' @export
cbind.redist_scorer <- function(..., deparse.level = 1) {
    fns <- list(...)
    stopifnot(all(sapply(fns, function(x) inherits(x, "redist_scorer"))))

    fn <- function(plans) {
        do.call(cbind, c(lapply(fns, function(fn) {
            fn(plans)
        }), list(deparse.level=deparse.level)))
    }
    class(fn) <- c("redist_scorer", "function")
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
`*.redist_scorer` <- function(x, fn2) {
    stopifnot(is.numeric(x) || inherits(x, "redist_scorer"))
    stopifnot(inherits(fn2, "redist_scorer"))

    if (is.numeric(x)) {
        rlang::fn_body(fn2) <- rlang::expr({!!x*!!rlang::fn_body(fn2)})
        return(fn2)
    } else {
        fn <- function(plans) { x(plans)*fn2(plans) }
        class(fn) <- c("redist_scorer", "function")
    }

    fn
}

#' @rdname scorer-arith
#'
#' @param fn1 a `redist_scorer` function, from [`scorers`]
#' @param fn2 a `redist_scorer` function, from [`scorers`]
#'
#' @export
`+.redist_scorer` <- function(fn1, fn2) {
    stopifnot(inherits(fn1, "redist_scorer"))
    stopifnot(inherits(fn2, "redist_scorer"))

    fn <- function(plans) { fn1(plans) + fn2(plans) }
    class(fn) <- c("redist_scorer", "function")
    fn
}

#' @rdname scorer-arith
#'
#' @param fn1 a `redist_scorer` function, from [`scorers`]
#' @param fn2 a `redist_scorer` function, from [`scorers`]
#'
#' @export
`-.redist_scorer` <- function(fn1, fn2) {
    stopifnot(inherits(fn1, "redist_scorer"))
    stopifnot(inherits(fn2, "redist_scorer"))

    fn <- function(plans) { fn1(plans) - fn2(plans) }
    class(fn) <- c("redist_scorer", "function")
    fn
}
