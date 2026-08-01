if (sys.nframe() == 0L) {
    devtools::load_all()
}

make_grid_case <- function(n) {
    shape <- geomander::checkerboard |>
        dplyr::filter(i < n, j < n) |>
        sf::st_set_crs(3857) |>
        dplyr::mutate(pop = n)

    list(
        map = redist_map(shape, ndists = n, pop_tol = 0.01),
        stripe_plan = as.integer(vctrs::vec_group_id(shape$i)),
        hierarchy = list(
            coarse = (as.integer(vctrs::vec_group_id(shape$i)) + 1L) %/% 2L,
            fine = as.integer(vctrs::vec_group_id(shape$i))
        )
    )
}

coassignment_distance <- function(old, new) {
    choose2 <- \(counts) sum(counts * (counts - 1) / 2)
    disagree <- choose2(tabulate(old)) +
        choose2(tabulate(new)) -
        2 * choose2(tabulate(vctrs::vec_group_id(interaction(old, new))))
    disagree / choose(length(old), 2)
}

transition_metrics <- function(plans) {
    plan_matrix <- get_plans_matrix(plans)
    if (ncol(plan_matrix) < 2L) {
        return(c(nonidentity_rate = NA_real_, mean_coassignment_distance = NA_real_))
    }

    distance <- vapply(
        seq_len(ncol(plan_matrix) - 1L),
        \(i) coassignment_distance(plan_matrix[, i], plan_matrix[, i + 1L]),
        numeric(1)
    )
    c(
        nonidentity_rate = mean(distance > 0),
        mean_coassignment_distance = mean(distance)
    )
}

edge_sequence <- function(plans, map) {
    with_edges <- dplyr::mutate(plans, edges_rem = comp_edges_rem(pl(), map))
    by_plan(with_edges$edges_rem, attr(map, 'ndists'))
}

edge_ess <- function(edges) {
    if (length(edges) < 3L || length(unique(edges)) < 2L) {
        return(0)
    }
    as.numeric(coda::effectiveSize(edges))
}

target_7x7 <- stats::setNames(
    c(
        0.020807687768601,
        0.0710143952315672,
        0.14977414247499,
        0.197857201120775,
        0.2074320577516,
        0.163717974693393,
        0.103810429703219,
        0.0529758100529919,
        0.0222509101447097,
        0.0076427174466794,
        0.00214991687533369,
        0.000477190469004618,
        0.0000799810543944961,
        0.00000903290674648099,
        0.000000552305514866464
    ),
    28:42
)

target_error <- function(edges, n) {
    if (n != 7L) {
        return(c(edge_tv = NA_real_, edge_rmse = NA_real_))
    }
    estimate <- table(factor(edges, levels = names(target_7x7))) / length(edges)
    difference <- as.numeric(estimate) - target_7x7
    c(edge_tv = sum(abs(difference)) / 2, edge_rmse = sqrt(mean(difference^2)))
}

run_sampler <- function(algorithm, case, init_plan, nsims, seed) {
    set.seed(seed)
    gc(FALSE)
    timing <- system.time({
        # compactness = 1 is the unadjusted spanning-tree target for the
        # established MCMC samplers. HDimer fixes that exponent at one.
        plans <- switch(
            algorithm,
            hier_dimer = redist_hier_dimer(
                case$map,
                nsims,
                hierarchy = list(),
                init_plan = init_plan,
                thin = 1L,
                refresh = 1L,
                silent = TRUE
            ),
            hier_dimer_refresh10 = redist_hier_dimer(
                case$map,
                nsims,
                hierarchy = list(),
                init_plan = init_plan,
                thin = 1L,
                refresh = 10L,
                silent = TRUE
            ),
            mergesplit = redist_mergesplit(
                case$map,
                nsims,
                warmup = 0L,
                thin = 1L,
                init_plan = init_plan,
                compactness = 1,
                init_name = FALSE,
                silent = TRUE
            ),
            cyclewalk = redist_cyclewalk(
                case$map,
                nsims,
                chains = 1L,
                warmup = 0L,
                thin = 1L,
                cycle_walk_frac = 0.1,
                init_plan = init_plan,
                compactness = 1,
                ncores = 1L,
                init_name = FALSE,
                silent = TRUE
            ),
            hier_dimer_nested = redist_hier_dimer(
                case$map,
                nsims,
                hierarchy = case$hierarchy,
                split_penalty = c(0, 0),
                cut_bias = c(4, 2),
                init_plan = case$stripe_plan,
                thin = 1L,
                refresh = 1L,
                silent = TRUE
            ),
            hier_dimer_nested_refresh10 = redist_hier_dimer(
                case$map,
                nsims,
                hierarchy = case$hierarchy,
                split_penalty = c(0, 0),
                cut_bias = c(4, 2),
                init_plan = case$stripe_plan,
                thin = 1L,
                refresh = 10L,
                silent = TRUE
            )
        )
    })
    elapsed <- timing[['elapsed']]
    cpu <- timing[['user.self']] + timing[['sys.self']]

    moves <- transition_metrics(plans)
    edges <- edge_sequence(plans, case$map)
    ess <- edge_ess(edges)
    error <- target_error(edges, attr(case$map, 'ndists'))
    if (grepl('nested', algorithm, fixed = TRUE)) {
        error <- c(edge_tv = NA_real_, edge_rmse = NA_real_)
    }
    diagnostics <- attr(plans, 'diagnostics')

    data.frame(
        algorithm = algorithm,
        elapsed = elapsed,
        cpu = cpu,
        saved_plans = length(edges),
        transitions_per_second = nsims / elapsed,
        transitions_per_cpu_second = nsims / cpu,
        changed_transitions_per_second = nsims * moves[['nonidentity_rate']] / elapsed,
        nonidentity_rate = moves[['nonidentity_rate']],
        mean_coassignment_distance = moves[['mean_coassignment_distance']],
        edge_ess = ess,
        edge_ess_per_second = ess / elapsed,
        edge_tv = error[['edge_tv']],
        edge_rmse = error[['edge_rmse']],
        backend_nonplanar = diagnostics$proposal_nonplanar %||% NA_real_,
        backend_numerical_failure = diagnostics$proposal_numerical_failure %||% NA_real_,
        backend_no_matching = diagnostics$proposal_no_matching %||% NA_real_,
        microseconds_tree_sampling = 1e6 *
            (diagnostics$seconds_tree_sampling %||% NA_real_) /
            nsims,
        microseconds_cut = 1e6 * (diagnostics$seconds_cut %||% NA_real_) / nsims,
        microseconds_fragments = 1e6 * (diagnostics$seconds_fragments %||% NA_real_) / nsims,
        microseconds_candidates = 1e6 * (diagnostics$seconds_candidates %||% NA_real_) / nsims,
        microseconds_matching = 1e6 * (diagnostics$seconds_matching %||% NA_real_) / nsims,
        microseconds_relabel = 1e6 * (diagnostics$seconds_relabel %||% NA_real_) / nsims
    )
}

run_case <- function(
    n,
    nsims,
    seeds,
    algorithms = c(
        'hier_dimer', 'hier_dimer_refresh10', 'mergesplit', 'cyclewalk',
        'hier_dimer_nested', 'hier_dimer_nested_refresh10'
    )
) {
    case <- make_grid_case(n)
    set.seed(7000L + n)
    initial <- get_plans_matrix(redist_smc(case$map, length(seeds), silent = TRUE))

    result <- do.call(
        rbind,
        lapply(seq_along(seeds), function(i) {
            do.call(
                rbind,
                lapply(
                    algorithms,
                    run_sampler,
                    case = case,
                    init_plan = initial[, i],
                    nsims = nsims,
                    seed = seeds[[i]]
                )
            )
        })
    )
    result$size <- sprintf('%dx%d', n, n)
    result$seed <- rep(seeds, each = length(algorithms))
    result
}

summarize_result <- function(results) {
    stats::aggregate(
        results[setdiff(names(results), c('algorithm', 'size', 'seed'))],
        by = results[c('size', 'algorithm')],
        FUN = mean,
        na.rm = TRUE
    )
}

if (sys.nframe() == 0L) {
    nsims <- as.integer(Sys.getenv('REDIST_HDIMER_BENCH_NSIMS', '10000'))
    n_replicates <- as.integer(Sys.getenv('REDIST_HDIMER_BENCH_REPS', '5'))
    seeds <- 9000L + seq_len(n_replicates)

    results <- rbind(
        run_case(6L, nsims, seeds),
        run_case(7L, nsims, seeds)
    )
    summary <- summarize_result(results)

    print(results)
    print(summary)

    if (
        any(results$backend_nonplanar > 0, na.rm = TRUE) ||
            any(results$backend_numerical_failure > 0, na.rm = TRUE) ||
            any(results$backend_no_matching > 0, na.rm = TRUE)
    ) {
        cli::cli_abort('The hierarchical split-dimer backend reported a failed matching step.')
    }
}
