#' Hierarchical split-dimer MCMC redistricting sampler
#'
#' `redist_hier_dimer()` uses a Markov Chain Monte Carlo algorithm based on
#' splitting and rematching district spanning trees.
#' The rematching step uses an exact weighted draw, avoiding enumeration and
#' post hoc filtering of candidate plans.
#'
#' This function draws samples from a target measure controlled by `map`,
#' `hierarchy`, `split_penalty`, `boundary_penalty`, and `constraints`.
#' Sampled plans satisfy the map's population and connectedness requirements.
#' When a hierarchy is supplied, district intersections with hierarchy units
#' must also be connected.
#' Constraints are evaluated on proposed plans and applied with a
#' Metropolis-Hastings correction.
#'
#' Hierarchy levels should be supplied from coarse to fine.
#' See [summary.redist_plans()] for sampler diagnostics.
#'
#' @param map A [redist_map] object.
#' @param nsims The number of samples to draw.
#' @param hierarchy An optional list of map column names or vectors, ordered
#'   from coarse to fine.
#' @param split_penalty Nonnegative target penalties, one for each hierarchy
#'   level.
#' @param cut_bias Nonnegative proposal biases, one for each hierarchy level.
#'   Larger values make the auxiliary cut step more likely to cut tree edges
#'   that cross that level's boundaries, without changing the target
#'   distribution.
#' @param boundary_penalty A nonnegative coefficient on the number of
#'   district-internal adjacency edges.
#'   Larger values favor plans with fewer graph-boundary edges.
#' @param constraints A [redist_constr()] object or a list containing
#'   information on sampling constraints.
#'   See [constraints] for more information.
#'   Constraints are evaluated on each proposed plan and applied with a
#'   Metropolis-Hastings correction.
#' @param l The number of district trees to cut and rematch per transition.
#'   Defaults to the smaller of seven and the number of districts.
#' @param init_plan An initial district assignment.
#'   Defaults to the map's existing plan, or a plan drawn by [redist_smc].
#'   When `hierarchy` is supplied, the initial plan must have connected
#'   intersections with every hierarchy unit.
#' @param init_name A name for the initial plan, or `FALSE` to not include it
#'   as a reference plan in the output.
#' @param thin The number of transitions between saved plans.
#' @param refresh The number of transitions between fresh draws of the
#'   selected district trees.
#'   The default redraws them every transition.
#' @param return_all Whether to return all saved plans or only the last saved
#'   plan.
#' @param verbose Whether to show detailed backend output.
#' @param silent Whether to suppress backend output.
#'
#' @return A [redist_plans] object containing the simulated plans.
#' @export
#'
#' @concept simulate
#' @md
#'
#' @examples
#' adj <- list(1L, c(0L, 2L, 3L), c(1L, 3L), c(1L, 2L))
#' map <- redist_map(pop = rep(1, 4L), ndists = 2L, pop_tol = 0.01, adj = adj)
#'
#' sampled_basic <- redist_hier_dimer(map, nsims = 10L,
#'     init_plan = c(1L, 1L, 2L, 2L), silent = TRUE)
#'
#' constr <- redist_constr(map)
#' constr <- add_constr_incumbency(constr, strength = 1,
#'     incumbents = c(1L, 3L))
#' sampled_constr <- redist_hier_dimer(map, nsims = 10L,
#'     constraints = constr, init_plan = c(1L, 1L, 2L, 2L), silent = TRUE)
#'
redist_hier_dimer <- function(
    map,
    nsims,
    hierarchy = list(),
    split_penalty = NULL,
    cut_bias = NULL,
    boundary_penalty = 0,
    constraints = list(),
    l = NULL,
    init_plan = NULL,
    init_name = NULL,
    thin = 1L,
    refresh = 1L,
    return_all = TRUE,
    verbose = FALSE,
    silent = FALSE
) {
    map <- validate_redist_map(map)
    adj <- get_adj(map)
    n_units <- nrow(map)
    ndists <- attr(map, 'ndists')
    if (nsims < 1) {
        cli::cli_abort('{.arg nsims} must be positive.')
    }
    nsims <- as.integer(nsims)
    if (boundary_penalty < 0) {
        cli::cli_abort('{.arg boundary_penalty} must be non-negative.')
    }
    thin <- as.integer(thin)
    if (thin < 1 || thin > nsims) {
        cli::cli_abort('{.arg thin} must be a positive integer no larger than {.arg nsims}.')
    }
    refresh <- as.integer(refresh)
    if (refresh < 1) {
        cli::cli_abort('{.arg refresh} must be a positive integer.')
    }
    if (is.null(l)) {
        l <- min(7L, ndists)
    }
    if (l < 2 || l > ndists) {
        cli::cli_abort('{.arg l} must be between 2 and the number of districts.')
    }
    l <- as.integer(l)

    n_levels <- length(hierarchy)
    if (is.null(split_penalty)) {
        split_penalty <- rep(0, n_levels)
    }
    if (is.null(cut_bias)) {
        cut_bias <- rep(1, n_levels)
    }
    hierarchy <- prepare_hier_dimer_hierarchy(hierarchy, map, adj)
    if (any(split_penalty < 0)) {
        cli::cli_abort('{.arg split_penalty} must be non-negative.')
    }
    if (length(split_penalty) != ncol(hierarchy)) {
        cli::cli_abort('{.arg split_penalty} must have one entry for each hierarchy level.')
    }
    split_penalty <- as.numeric(split_penalty)
    if (any(cut_bias < 0)) {
        cli::cli_abort('{.arg cut_bias} must be non-negative.')
    }
    if (length(cut_bias) != ncol(hierarchy)) {
        cli::cli_abort('{.arg cut_bias} must have one entry for each hierarchy level.')
    }
    cut_bias <- as.numeric(cut_bias)

    if (!inherits(constraints, 'redist_constr')) {
        constraints <- new_redist_constr(rlang::eval_tidy(rlang::enquo(constraints), map))
    }
    constraints <- as.list(validate_redist_constr(constraints))

    existing_col <- attr(map, 'existing_col')
    if (is.null(init_plan)) {
        if (!is.null(existing_col)) {
            init_plan <- vctrs::vec_group_id(get_existing(map))
            if (is.null(init_name)) init_name <- existing_col
        } else {
            init_plan <- as.integer(get_plans_matrix(
                redist_smc(
                    map,
                    1L,
                    resample = FALSE,
                    ref_name = FALSE,
                    silent = TRUE,
                    ncores = 1L
                )
            )[, 1])
            if (is.null(init_name)) init_name <- '<init>'
        }
    } else if (is.null(init_name)) {
        init_name <- '<init>'
    }
    if (!is.numeric(init_plan) || length(init_plan) != n_units || anyNA(init_plan)) {
        cli::cli_abort('{.arg init_plan} must be a non-missing numeric vector as long as `map`.')
    }
    init_plan <- vctrs::vec_group_id(init_plan)
    if (max(init_plan) != ndists) {
        cli::cli_abort('{.arg init_plan} must have the same number of districts as `map`.')
    }
    if (any(contiguity(adj, init_plan) != 1L)) {
        cli::cli_abort('{.arg init_plan} must have contiguous districts.')
    }

    pop_bounds <- attr(map, 'pop_bounds')
    pop <- map[[attr(map, 'pop_col')]]
    init_pop <- pop_tally(matrix(init_plan, ncol = 1L), pop, ndists)
    if (any(init_pop < pop_bounds[1L]) || any(init_pop > pop_bounds[3L])) {
        cli::cli_abort('Provided initialization does not meet population bounds.')
    }

    verbosity <- 1L
    if (silent) {
        verbosity <- 0L
    } else if (verbose) {
        verbosity <- 3L
    }
    started <- Sys.time()
    algout <- hier_dimer_plans(
        nsims,
        adj,
        init_plan,
        pop,
        hierarchy,
        split_penalty,
        cut_bias,
        boundary_penalty,
        pop_bounds[2L],
        pop_bounds[1L],
        pop_bounds[3L],
        constraints,
        thin,
        refresh,
        l,
        verbosity
    )
    elapsed <- as.numeric(Sys.time() - started, units = 'secs')
    if (is.null(algout$plans) || !is.matrix(algout$plans) || nrow(algout$plans) != n_units) {
        cli::cli_abort('The hierarchical split-dimer backend returned invalid plans.')
    }

    plans <- algout$plans
    storage.mode(plans) <- 'integer'
    if (ncol(plans) > 1L) {
        plans <- plans[, -1L, drop = FALSE]
    }
    if (!return_all) {
        plans <- plans[, ncol(plans), drop = FALSE]
    }
    changed <- as.logical(algout$transition_changed %||% logical())
    if (length(changed) > 0L) {
        changed <- changed[-1L]
        if (!return_all) {
            changed <- tail(changed, 1L)
        }
    }

    transition_change_rate <- NA_real_
    proposal_attempts <- algout$diagnostics$proposal_attempts %||% 0L
    if (proposal_attempts > 0L) {
        transition_change_rate <-
            (algout$diagnostics$proposal_changed %||% 0L) / proposal_attempts
    }

    out <- new_redist_plans(
        plans,
        map,
        'hier_dimer',
        NULL,
        FALSE,
        ndists = ndists,
        compactness = 1,
        boundary_penalty = boundary_penalty,
        l = l,
        hierarchy = hierarchy,
        split_penalty = split_penalty,
        cut_bias = cut_bias,
        constraints = constraints,
        version = packageVersion('redist'),
        diagnostics = c(list(runtime = elapsed), algout$diagnostics %||% list()),
        transition_change_rate = transition_change_rate
    )
    if (length(changed) == ncol(plans)) {
        out <- dplyr::mutate(out, mcmc_changed = rep(changed, each = ndists))
    }
    if (!is.null(init_name) && !isFALSE(init_name)) {
        out <- add_reference(out, init_plan, init_name)
    }
    out
}

prepare_hier_dimer_hierarchy <- function(hierarchy, map, adj) {
    if (!is.list(hierarchy)) {
        cli::cli_abort('{.arg hierarchy} must be a list of vectors or map column names.')
    }
    n_units <- nrow(map)
    if (length(hierarchy) == 0L) {
        return(matrix(integer(), nrow = n_units, ncol = 0L))
    }

    levels <- vector('list', length(hierarchy))
    for (i in seq_along(hierarchy)) {
        level <- hierarchy[[i]]
        if (is.character(level) && length(level) == 1L) {
            if (!(level %in% names(map))) {
                cli::cli_abort('Hierarchy column {.val {level}} is not present in {.arg map}.')
            }
            level <- map[[level]]
        }
        if (length(level) != n_units || anyNA(level)) {
            cli::cli_abort('Hierarchy level {i} must be a non-missing vector as long as `map`.')
        }
        level <- vctrs::vec_group_id(level)
        if (i > 1L) {
            level <- vctrs::vec_group_id(paste(levels[[i - 1L]], level, sep = ':'))
        }
        component <- contiguity(adj, level)
        level <- vctrs::vec_group_id(paste(level, component, sep = ':'))
        levels[[i]] <- as.integer(level)
    }
    out <- do.call(cbind, levels)
    level_names <- names(hierarchy)
    if (is.null(level_names)) {
        level_names <- rep('', length(hierarchy))
    }
    missing_names <- !nzchar(level_names)
    level_names[missing_names] <- paste0('level_', which(missing_names))
    colnames(out) <- make.unique(level_names)
    out
}
