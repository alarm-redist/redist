devtools::load_all()

make_grid_map <- function() {
    grid_sf <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 4, ymax = 4)) |>
        sf::st_as_sfc() |>
        sf::st_make_grid(n = c(4, 4)) |>
        sf::st_sf(geometry = _) |>
        sf::st_set_crs(3857)
    grid_sf$adj <- redist.adjacency(grid_sf)
    grid_sf$pop <- rep(1L, nrow(grid_sf))
    centers <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(grid_sf)))
    east_west <- ifelse(centers[, 1] < 2, 'west', 'east')
    north_south <- ifelse(centers[, 2] < 2, 'south', 'north')
    grid_sf$municipality <- paste(east_west, north_south, sep = '_')
    grid_sf$county <- east_west
    grid_sf$init <- vctrs::vec_group_id(grid_sf$municipality)

    redist_map(
        grid_sf,
        existing_plan = init,
        pop_bounds = c(3.99, 4, 4.01),
        adj = grid_sf$adj
    )
}

coassignment_distance <- function(old, new) {
    old_same <- outer(old, old, '==')
    new_same <- outer(new, new, '==')
    mean(old_same[upper.tri(old_same)] != new_same[upper.tri(new_same)])
}

split_count <- function(plan, group) {
    sum(vapply(split(plan, group), \(x) length(unique(x)) - 1L, integer(1)))
}

intersection_disconnected <- function(adj, plan, group) {
    intersection <- vctrs::vec_group_id(interaction(plan, group, drop = TRUE))
    any(contiguity(adj, intersection) > 1L)
}

plan_checks <- function(plans, map, hierarchy) {
    plan_m <- get_plans_matrix(plans)
    pop_bounds <- attr(map, 'pop_bounds')
    pop <- map[[attr(map, 'pop_col')]]
    ndists <- attr(map, 'ndists')
    adj <- get_adj(map)

    data.frame(
        population_failure = apply(plan_m, 2L, function(plan) {
            totals <- pop_tally(matrix(plan, ncol = 1L), pop, ndists)
            any(totals < pop_bounds[1L] | totals > pop_bounds[3L])
        }),
        contiguity_failure = apply(plan_m, 2L, function(plan) {
            any(contiguity(adj, plan) > 1L)
        }),
        hierarchy_intersection_failure = apply(plan_m, 2L, function(plan) {
            any(vapply(
                hierarchy,
                function(group) {
                    intersection_disconnected(adj, plan, group)
                },
                logical(1)
            ))
        }),
        municipality_splits = apply(
            plan_m,
            2L,
            split_count,
            group = hierarchy[[length(hierarchy)]]
        ),
        county_splits = apply(plan_m, 2L, split_count, group = hierarchy[[1L]])
    )
}

transition_checks <- function(plans) {
    plan_m <- get_plans_matrix(plans)
    if (ncol(plan_m) < 2L) {
        return(data.frame(
            nonidentity = logical(),
            coassignment_distance = numeric()
        ))
    }

    do.call(
        rbind,
        lapply(seq_len(ncol(plan_m) - 1L), function(i) {
            old <- plan_m[, i]
            new <- plan_m[, i + 1L]
            data.frame(
                nonidentity = coassignment_distance(old, new) > 0,
                coassignment_distance = coassignment_distance(old, new)
            )
        })
    )
}

diagnostic_value <- function(plans, name) {
    if (name %in% names(plans) && !is.null(plans[[name]])) {
        return(as.numeric(plans[[name]]))
    }
    if (!is.null(attr(plans, name))) {
        return(as.numeric(attr(plans, name)))
    }

    diagnostics <- attr(plans, 'diagnostics')
    if (is.null(diagnostics) || length(diagnostics) == 0L) {
        return(NA_real_)
    }

    if (!(name %in% names(diagnostics))) {
        return(NA_real_)
    }
    value <- diagnostics[[name]]
    if (is.null(value)) {
        return(NA_real_)
    }
    as.numeric(value)
}

failed_diagnostic <- function(method, seed, error, runtime) {
    data.frame(
        method = method,
        seed = seed,
        completed = FALSE,
        error = error,
        proposal_attempts = NA_real_,
        proposal_changed = NA_real_,
        proposal_identity = NA_real_,
        proposal_nonplanar = NA_real_,
        proposal_numerical_failure = NA_real_,
        proposal_no_matching = NA_real_,
        mean_vertices_changed = NA_real_,
        mean_matching_edges_changed = NA_real_,
        mean_old_matching_prob = NA_real_,
        nonidentity_rate = NA_real_,
        coassignment_distance = NA_real_,
        population_failures = NA_real_,
        contiguity_failures = NA_real_,
        hierarchy_intersection_failures = NA_real_,
        municipality_splits = NA_real_,
        county_splits = NA_real_,
        runtime = runtime
    )
}

run_diagnostic <- function(
    map,
    proposal_hierarchy,
    check_hierarchy,
    split_penalty,
    cut_bias,
    seed,
    method,
    nsims = 200L
) {
    started <- proc.time()[['elapsed']]
    out <- tryCatch(
        {
            set.seed(seed)
            redist_hier_dimer(
                map = map,
                nsims = nsims,
                hierarchy = proposal_hierarchy,
                split_penalty = split_penalty,
                cut_bias = cut_bias,
                init_plan = get_existing(map),
                thin = 1L,
                refresh = 1L,
                return_all = TRUE,
                silent = TRUE
            )
        },
        error = \(e) e
    )
    elapsed <- proc.time()[['elapsed']] - started

    if (inherits(out, 'error')) {
        return(failed_diagnostic(method, seed, conditionMessage(out), elapsed))
    }

    transitions <- transition_checks(out)
    checks <- plan_checks(out, map, check_hierarchy)
    runtime <- diagnostic_value(out, 'runtime')
    if (is.na(runtime)) {
        runtime <- elapsed
    }
    data.frame(
        method = method,
        seed = seed,
        completed = TRUE,
        error = NA_character_,
        proposal_attempts = diagnostic_value(out, 'proposal_attempts'),
        proposal_changed = diagnostic_value(out, 'proposal_changed'),
        proposal_identity = diagnostic_value(out, 'proposal_identity'),
        proposal_nonplanar = diagnostic_value(out, 'proposal_nonplanar'),
        proposal_numerical_failure = diagnostic_value(out, 'proposal_numerical_failure'),
        proposal_no_matching = diagnostic_value(out, 'proposal_no_matching'),
        mean_vertices_changed = diagnostic_value(out, 'mean_vertices_changed'),
        mean_matching_edges_changed = diagnostic_value(out, 'mean_matching_edges_changed'),
        mean_old_matching_prob = diagnostic_value(out, 'mean_old_matching_prob'),
        nonidentity_rate = mean(transitions$nonidentity),
        coassignment_distance = mean(transitions$coassignment_distance),
        population_failures = sum(checks$population_failure),
        contiguity_failures = sum(checks$contiguity_failure),
        hierarchy_intersection_failures = sum(checks$hierarchy_intersection_failure),
        municipality_splits = mean(checks$municipality_splits),
        county_splits = mean(checks$county_splits),
        runtime = runtime
    )
}

run_comparison <- function(
    map,
    hierarchy,
    split_penalty,
    label,
    cut_bias = rep(1, length(hierarchy)),
    seeds = 1001:1010,
    nsims = 200L
) {
    hierarchy_args <- list(
        hierarchy_targeted = list(
            hierarchy = hierarchy,
            split_penalty = split_penalty,
            cut_bias = cut_bias
        ),
        hierarchy_support_only = list(
            hierarchy = hierarchy,
            split_penalty = rep(0, length(hierarchy)),
            cut_bias = cut_bias
        ),
        hierarchy_disabled = list(
            hierarchy = list(),
            split_penalty = numeric(),
            cut_bias = numeric()
        )
    )
    results <- lapply(names(hierarchy_args), function(method) {
        settings <- hierarchy_args[[method]]
        do.call(
            rbind,
            lapply(seeds, function(seed) {
                run_diagnostic(
                    map = map,
                    proposal_hierarchy = settings$hierarchy,
                    check_hierarchy = hierarchy,
                    split_penalty = settings$split_penalty,
                    cut_bias = settings$cut_bias,
                    seed = seed,
                    method = method,
                    nsims = nsims
                )
            })
        )
    })
    results <- do.call(rbind, results)
    results$example <- label
    results
}

data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
iowa_hierarchy <- list(municipality = iowa$region)

grid_map <- make_grid_map()
grid_hierarchy <- list(
    county = grid_map$county,
    municipality = grid_map$municipality
)

n_seeds <- as.integer(Sys.getenv('REDIST_HDIMER_NSEEDS', '10'))
n_sims <- as.integer(Sys.getenv('REDIST_HDIMER_NSIMS', '200'))
diagnostic_seeds <- 1000L + seq_len(n_seeds)

iowa_results <- run_comparison(
    iowa_map,
    hierarchy = iowa_hierarchy,
    split_penalty = 1,
    cut_bias = 3,
    label = 'Iowa: one level',
    seeds = diagnostic_seeds,
    nsims = n_sims
)
grid_results <- run_comparison(
    grid_map,
    hierarchy = grid_hierarchy,
    split_penalty = c(2, 1),
    cut_bias = c(4, 2),
    label = '4 by 4 grid: two nested levels',
    seeds = diagnostic_seeds,
    nsims = n_sims
)

results <- rbind(iowa_results, grid_results)
summary <- aggregate(
    results[setdiff(names(results), c('method', 'seed', 'error', 'example'))],
    by = results[c('example', 'method')],
    FUN = \(x) mean(x, na.rm = TRUE)
)

print(results)
print(summary)

if (
    !all(results$completed) ||
        any(results$population_failures > 0, na.rm = TRUE) ||
        any(results$contiguity_failures > 0, na.rm = TRUE) ||
        any(
            results$hierarchy_intersection_failures[
                results$method != 'hierarchy_disabled'
            ] >
                0,
            na.rm = TRUE
        )
) {
    cli::cli_abort('Hierarchical split-dimer diagnostics found a failed run or invalid plan.')
}
