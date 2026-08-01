devtools::load_all()
source('explore/benchmark_hier_dimer_grids.R')

grid_edges <- function(map) {
    adjacency <- get_adj(map)
    from <- rep(seq_along(adjacency), lengths(adjacency))
    to <- unlist(adjacency, use.names = FALSE) + 1L
    keep <- from < to
    from <- from[keep]
    to <- to[keep]

    data.frame(
        from = from,
        to = to,
        axis = ifelse(map$i[from] != map$i[to], 'x', 'y')
    )
}

grid_statistics <- function(plans, map, block_size = 10000L) {
    plan_matrix <- get_plans_matrix(plans)
    edges <- grid_edges(map)
    n_plans <- ncol(plan_matrix)
    cut_x <- integer(n_plans)
    cut_y <- integer(n_plans)

    for (first in seq.int(1L, n_plans, by = block_size)) {
        columns <- first:min(first + block_size - 1L, n_plans)
        for (axis in c('x', 'y')) {
            selected <- edges$axis == axis
            count <- colSums(
                plan_matrix[edges$from[selected], columns, drop = FALSE] !=
                    plan_matrix[edges$to[selected], columns, drop = FALSE]
            )
            if (axis == 'x') {
                cut_x[columns] <- count
            } else {
                cut_y[columns] <- count
            }
        }
    }

    data.frame(
        edges = cut_x + cut_y,
        axis = cut_x - cut_y,
        abs_axis = abs(cut_x - cut_y)
    )
}

coassignment_signatures <- function(plans) {
    plan_matrix <- as.matrix(plans)
    if (inherits(plans, 'redist_plans')) {
        plan_matrix <- get_plans_matrix(plans)
    }
    pairs <- which(
        upper.tri(matrix(FALSE, nrow(plan_matrix), nrow(plan_matrix))),
        arr.ind = TRUE
    )
    plan_matrix[pairs[, 1L], , drop = FALSE] == plan_matrix[pairs[, 2L], , drop = FALSE]
}

diverse_initial_plans <- function(case, n_chains, seed) {
    set.seed(seed)
    pool <- get_plans_matrix(redist_smc(case$map, max(200L, 20L * n_chains), silent = TRUE))
    vertical <- case$stripe_plan
    horizontal <- as.integer(vctrs::vec_group_id(case$map$j))
    starts <- cbind(vertical, horizontal)
    pool_signature <- coassignment_signatures(pool)
    selected_signature <- coassignment_signatures(starts)

    minimum_distance <- rep(Inf, ncol(pool))
    for (i in seq_len(ncol(selected_signature))) {
        minimum_distance <- pmin(
            minimum_distance,
            colMeans(pool_signature != selected_signature[, i])
        )
    }
    while (ncol(starts) < n_chains) {
        next_plan <- which.max(minimum_distance)
        starts <- cbind(starts, pool[, next_plan])
        minimum_distance <- pmin(
            minimum_distance,
            colMeans(pool_signature != pool_signature[, next_plan])
        )
        minimum_distance[next_plan] <- -Inf
    }
    starts
}

truth_6x6 <- function(case) {
    path <- '../mmss/data-out/enum_6x6_into_6.rds'
    if (!file.exists(path)) {
        return(NULL)
    }

    enumeration <- readRDS(path)
    log_weight <- comp_log_st(enumeration, case$map)
    weight <- exp(log_weight - max(log_weight))
    weight <- 6 * weight / sum(weight)
    keep <- enumeration[['district']] == 1L
    plan_statistics <- grid_statistics(enumeration, case$map)
    edge <- enumeration[['edges_rem']][keep]
    plan_weight <- weight[keep]

    list(
        edge_distribution = tapply(plan_weight, edge, sum),
        edges = weighted.mean(edge, plan_weight),
        central = weighted.mean(edge == 22, plan_weight),
        upper = weighted.mean(edge >= 27, plan_weight),
        lower = weighted.mean(edge == 18, plan_weight),
        axis = weighted.mean(plan_statistics$axis, plan_weight),
        abs_axis = weighted.mean(plan_statistics$abs_axis, plan_weight)
    )
}

truth_7x7 <- function() {
    edge <- as.integer(names(target_7x7))
    list(
        edge_distribution = target_7x7,
        edges = weighted.mean(edge, target_7x7),
        central = target_7x7[['32']],
        upper = sum(target_7x7[names(target_7x7) >= '35']),
        lower = sum(target_7x7[names(target_7x7) <= '29']),
        axis = 0,
        abs_axis = NA_real_
    )
}

observable <- function(statistics, name, n) {
    central <- 32L
    upper <- 35L
    if (n == 6L) {
        central <- 22L
        upper <- 27L
    }
    switch(
        name,
        edges = statistics$edges,
        central = statistics$edges == central,
        upper = statistics$edges >= upper,
        lower = statistics$edges == 18L & n == 6L | statistics$edges <= 29L & n == 7L,
        axis = statistics$axis,
        abs_axis = statistics$abs_axis
    )
}

mcmc_interval <- function(draws, conf = 0.95) {
    chains <- coda::mcmc.list(lapply(seq_len(ncol(draws)), function(i) {
        coda::mcmc(draws[, i])
    }))
    standard_error <- unname(suppressWarnings(summary(chains)$statistics['Time-series SE']))
    estimate <- unname(mean(draws))
    half_width <- stats::qt(
        1 - (1 - conf) / 2,
        df = length(draws) - 1L
    ) *
        standard_error
    c(estimate = estimate, width = 2 * half_width, lower = estimate - half_width, upper = estimate + half_width)
}

summarize_draws <- function(chains, prefixes, runtime, cpu_runtime, algorithm, n, truth) {
    metrics <- c('edges', 'central', 'upper', 'lower', 'axis', 'abs_axis')
    do.call(
        rbind,
        lapply(seq_along(prefixes), function(budget) {
            prefix <- prefixes[[budget]]
            do.call(
                rbind,
                lapply(metrics, function(metric) {
                    draws <- vapply(
                        chains,
                        \(chain) observable(chain, metric, n)[seq_len(prefix)],
                        numeric(prefix)
                    )
                    group <- rep(seq_len(ncol(draws)), each = nrow(draws))
                    rhat <- suppressWarnings(tryCatch(
                        diag_rhat(as.vector(draws), group, split = TRUE),
                        error = \(e) NA_real_
                    ))
                    interval <- suppressWarnings(tryCatch(
                        mcmc_interval(draws),
                        error = function(e) {
                            c(
                    estimate = NA_real_,
                    width = NA_real_,
                    lower = NA_real_,
                    upper = NA_real_
                )
                        }
                    ))
                    target <- truth[[metric]] %||% NA_real_

                    data.frame(
                        size = sprintf('%dx%d', n, n),
                        algorithm = algorithm,
                        budget = budget,
                        transitions = prefix,
                        seconds_per_chain = mean(runtime) * prefix / max(prefixes),
                        cpu_seconds_per_chain = mean(cpu_runtime) * prefix / max(prefixes),
                        metric = metric,
                        rhat = rhat,
                        estimate = interval[['estimate']],
                        ci_width = interval[['width']],
                        truth = target,
                        absolute_error = abs(interval[['estimate']] - target),
                        covers_truth = interval[['lower']] <= target &
                            interval[['upper']] >= target
                    )
                })
            )
        })
    )
}

histogram_error <- function(chains, prefixes, runtime, cpu_runtime, algorithm, n, truth) {
    target <- truth$edge_distribution
    bins <- as.integer(names(target))
    do.call(
        rbind,
        lapply(seq_along(prefixes), function(budget) {
            prefix <- prefixes[[budget]]
            edges <- unlist(lapply(chains, \(chain) chain$edges[seq_len(prefix)]))
            estimate <- table(factor(edges, levels = bins)) / length(edges)
            difference <- as.numeric(estimate) - target
            data.frame(
                size = sprintf('%dx%d', n, n),
                algorithm = algorithm,
                budget = budget,
                transitions = prefix,
                seconds_per_chain = mean(runtime) * prefix / max(prefixes),
                cpu_seconds_per_chain = mean(cpu_runtime) * prefix / max(prefixes),
                edge_tv = sum(abs(difference)) / 2,
                edge_rmse = sqrt(mean(difference^2))
            )
        })
    )
}

run_chain <- function(algorithm, case, init_plan, nsims, seed) {
    set.seed(seed)
    gc(FALSE)
    timing <- system.time({
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
            )
        )
    })
    list(
        statistics = grid_statistics(plans, case$map),
        runtime = timing[['elapsed']],
        cpu_runtime = timing[['user.self']] + timing[['sys.self']]
    )
}

run_convergence_case <- function(n, schedules, n_chains) {
    case <- make_grid_case(n)
    starts <- diverse_initial_plans(case, n_chains, 8000L + n)
    truth <- truth_7x7()
    if (n == 6L) {
        truth <- truth_6x6(case)
    }
    if (is.null(truth)) {
        cli::cli_abort('The exact 6x6 enumeration is required for this benchmark.')
    }

    result <- list()
    histogram <- list()
    for (algorithm in names(schedules)) {
        prefixes <- schedules[[algorithm]]
        chains <- vector('list', n_chains)
        runtime <- numeric(n_chains)
        cpu_runtime <- numeric(n_chains)
        for (chain in seq_len(n_chains)) {
            cli::cli_inform(
                '{n}x{n} {algorithm} chain {chain}/{n_chains}',
                .frequency = 'once',
                .frequency_id = paste(n, algorithm, chain, sep = '-')
            )
            sampled <- run_chain(
                algorithm,
                case,
                starts[, chain],
                max(prefixes),
                seed = 100000L * n + 1000L * match(algorithm, names(schedules)) + chain
            )
            chains[[chain]] <- sampled$statistics
            runtime[[chain]] <- sampled$runtime
            cpu_runtime[[chain]] <- sampled$cpu_runtime
        }
        result[[algorithm]] <- summarize_draws(
            chains,
            prefixes,
            runtime,
            cpu_runtime,
            algorithm,
            n,
            truth
        )
        histogram[[algorithm]] <- histogram_error(
            chains,
            prefixes,
            runtime,
            cpu_runtime,
            algorithm,
            n,
            truth
        )
        rm(chains)
        gc(FALSE)
    }
    list(
        metrics = do.call(rbind, result),
        histogram = do.call(rbind, histogram)
    )
}

scale_schedule <- function(schedule, scale) {
    lapply(schedule, function(x) {
        value <- pmax(10L, as.integer(round(x * scale)))
        value + value %% 2L
    })
}

n_chains <- as.integer(Sys.getenv('REDIST_HDIMER_CONV_CHAINS', '8'))
scale <- as.numeric(Sys.getenv('REDIST_HDIMER_CONV_SCALE', '1'))
schedules <- scale_schedule(
    list(
    hier_dimer = c(1000L, 2500L, 5000L, 10000L),
    hier_dimer_refresh10 = c(1300L, 3250L, 6500L, 13000L),
    mergesplit = c(10000L, 22500L, 45000L, 90000L),
    cyclewalk = c(20000L, 50000L, 100000L, 200000L)
),
    scale
)

benchmark <- lapply(c(6L, 7L), run_convergence_case, schedules, n_chains)
metrics <- do.call(rbind, lapply(benchmark, `[[`, 'metrics'))
histogram <- do.call(rbind, lapply(benchmark, `[[`, 'histogram'))

curve <- merge(
    stats::aggregate(
        rhat ~ size +
            algorithm +
            budget +
            transitions +
            seconds_per_chain +
            cpu_seconds_per_chain,
        metrics,
        function(x) {
            if (all(is.na(x))) {
                return(NA_real_)
            }
            max(x, na.rm = TRUE)
        }
    ),
    histogram,
    by = c(
        'size', 'algorithm', 'budget', 'transitions', 'seconds_per_chain',
        'cpu_seconds_per_chain'
    )
)
names(curve)[names(curve) == 'rhat'] <- 'maximum_rhat'

print(curve)
print(metrics[metrics$budget == max(metrics$budget), ])

output_path <- Sys.getenv('REDIST_HDIMER_CONV_OUT', '')
if (nzchar(output_path)) {
    saveRDS(
        list(curve = curve, metrics = metrics, histogram = histogram),
        output_path
    )
}
