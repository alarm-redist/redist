test_that('hierarchy preprocessing splits disconnected units and is laminar', {
    adj <- list(c(1L), c(0L, 2L), c(1L, 3L), c(2L, 4L), c(3L))
    map <- redist_map(pop = rep(1, 5L), ndists = 2L, pop_tol = 0.5, adj = adj)
    map$coarse <- c('a', 'a', 'bridge', 'a', 'a')
    map$fine <- c('x', 'x', 'bridge', 'x', 'x')

    hierarchy <- prepare_hier_dimer_hierarchy(
        list('coarse', 'fine'),
        map,
        get_adj(map)
    )

    expect_equal(ncol(hierarchy), 2L)
    expect_identical(hierarchy[1L] == hierarchy[4L], FALSE)
    parent_counts <- vapply(
        split(hierarchy[, 1], hierarchy[, 2]),
        \(x) length(unique(x)),
        integer(1)
    )
    expect_identical(unname(parent_counts), rep(1L, length(parent_counts)))
})

test_that('hierarchical split-dimer draws satisfy every hard constraint', {
    adj <- list(
        c(1L, 3L), c(0L, 2L, 4L), c(1L, 5L),
        c(0L, 4L), c(1L, 3L, 5L), c(2L, 4L)
    )
    map <- redist_map(pop = rep(1, 6L), ndists = 3L, pop_tol = 0.01, adj = adj)
    coarse <- c(1L, 1L, 2L, 1L, 1L, 2L)
    fine <- c(1L, 2L, 3L, 1L, 2L, 3L)
    init <- c(1L, 1L, 2L, 3L, 3L, 2L)

    set.seed(20260729)
    out <- redist_hier_dimer(
        map,
        60L,
        hierarchy = list(coarse, fine),
        split_penalty = c(1, 0.5),
        init_plan = init,
        silent = TRUE
    )
    sampled <- subset_sampled(out)
    plans <- get_plans_matrix(sampled)
    processed <- attr(out, 'hierarchy')
    population_valid <- apply(plans, 2L, \(plan) {
        all(pop_tally(matrix(plan, ncol = 1L), rep(1, 6L), 3L) == 2)
    })
    district_contiguous <- apply(plans, 2L, \(plan) {
        all(contiguity(adj, plan) == 1L)
    })
    hierarchy_contiguous <- apply(plans, 2L, \(plan) {
        all(vapply(
            seq_len(ncol(processed)),
            \(level) {
                intersections <- vctrs::vec_group_id(
                    interaction(plan, processed[, level], drop = TRUE)
                )
                all(contiguity(adj, intersections) == 1L)
            },
            logical(1)
        ))
    })

    expect_equal(dim(plans), c(6L, 60L))
    expect_identical(unname(population_valid), rep(TRUE, 60L))
    expect_identical(unname(district_contiguous), rep(TRUE, 60L))
    expect_identical(unname(hierarchy_contiguous), rep(TRUE, 60L))
    expect_false(anyNA(sampled$mcmc_changed))
    expect_equal(attr(out, 'diagnostics')$proposal_attempts, 60L)
    expect_equal(attr(out, 'diagnostics')$proposal_nonplanar, 0L)
    expect_equal(attr(out, 'diagnostics')$proposal_numerical_failure, 0L)
    expect_equal(attr(out, 'diagnostics')$proposal_subset_dp, 60L)
})

test_that('sampler retains the polynomial matching path above the DP threshold', {
    side <- 4L
    adj <- lapply(seq_len(side^2), \(vertex) {
        row <- (vertex - 1L) %/% side + 1L
        column <- (vertex - 1L) %% side + 1L
        neighbors <- c(
            if (row > 1L) vertex - side,
            if (row < side) vertex + side,
            if (column > 1L) vertex - 1L,
            if (column < side) vertex + 1L
        )
        as.integer(neighbors - 1L)
    })
    map <- redist_map(pop = rep(1, side^2), ndists = 8L, pop_tol = 0.01, adj = adj)
    init <- rep(seq_len(8L), each = 2L)

    set.seed(8172)
    out <- redist_hier_dimer(
        map,
        20L,
        l = 8L,
        init_plan = init,
        silent = TRUE
    )
    plans <- get_plans_matrix(subset_sampled(out))
    contiguous <- apply(plans, 2L, \(plan) all(contiguity(adj, plan) == 1L))

    expect_equal(attr(out, 'diagnostics')$proposal_subset_dp, 0L)
    expect_equal(attr(out, 'diagnostics')$proposal_nonplanar, 0L)
    expect_equal(attr(out, 'diagnostics')$proposal_numerical_failure, 0L)
    expect_identical(unname(contiguous), rep(TRUE, 20L))
})

test_that('unconstrained tiny-grid target agrees with exact enumeration', {
    adj <- list(
        c(1L, 3L), c(0L, 2L, 4L), c(1L, 5L),
        c(0L, 4L), c(1L, 3L, 5L), c(2L, 4L)
    )
    map <- redist_map(pop = rep(1, 6L), ndists = 3L, pop_tol = 0.01, adj = adj)
    init <- c(1L, 1L, 2L, 3L, 3L, 2L)

    set.seed(9182)
    out <- redist_hier_dimer(
        map,
        900L,
        l = 2L,
        refresh = 3L,
        init_plan = init,
        silent = TRUE
    )
    plans <- get_plans_matrix(subset_sampled(out))
    keys <- apply(plans, 2L, plan_partition_key)
    frequencies <- table(keys) / length(keys)
    saved <- cbind(init, plans)
    changed_label_counts <- vapply(
        seq_len(ncol(saved) - 1L),
        \(i) {
            changed_units <- saved[, i] != saved[, i + 1L]
            length(unique(c(saved[changed_units, i], saved[changed_units, i + 1L])))
        },
        integer(1)
    )

    # The 2-by-3 grid has exactly three domino tilings. Every district tree and
    # every plan has the same target weight, so their probabilities are 1/3.
    expect_length(frequencies, 3L)
    expect_lt(max(abs(as.numeric(frequencies) - 1 / 3)), 0.06)
    expect_lte(max(changed_label_counts), 2L)
    expect_equal(attr(out, 'diagnostics')$l, 2L)
})

test_that('redist_constr changes the target with an MH correction', {
    adj <- list(
        c(1L, 3L), c(0L, 2L, 4L), c(1L, 5L),
        c(0L, 4L), c(1L, 3L, 5L), c(2L, 4L)
    )
    map <- redist_map(pop = rep(1, 6L), ndists = 3L, pop_tol = 0.01, adj = adj)
    oracle_plans <- list(
        c(1L, 1L, 2L, 3L, 3L, 2L),
        c(1L, 2L, 3L, 1L, 2L, 3L),
        c(1L, 2L, 2L, 1L, 3L, 3L)
    )
    oracle_keys <- vapply(oracle_plans, plan_partition_key, character(1))
    oracle_penalty <- vapply(
        oracle_plans,
        \(plan) as.numeric(plan[[1L]] == plan[[2L]]),
        numeric(1)
    )
    oracle_probabilities <- exp(-1.5 * oracle_penalty)
    oracle_probabilities <- oracle_probabilities / sum(oracle_probabilities)

    constraints <- add_constr_custom(
        redist_constr(map),
        strength = 1.5,
        fn = function(plan, distr) {
            as.numeric(all(c(1L, 2L) %in% which(plan == distr)))
        }
    )

    set.seed(20260801)
    out <- redist_hier_dimer(
        map,
        4000L,
        l = 2L,
        refresh = 3L,
        init_plan = oracle_plans[[1L]],
        constraints = constraints,
        silent = TRUE
    )
    sampled_keys <- apply(get_plans_matrix(subset_sampled(out)), 2L, plan_partition_key)
    observed <- table(factor(sampled_keys, levels = oracle_keys)) / length(sampled_keys)

    expect_lt(sum(abs(as.numeric(observed) - oracle_probabilities)) / 2, 0.08)
    expect_gt(attr(out, 'diagnostics')$constraint_proposals, 0L)
    expect_gt(attr(out, 'diagnostics')$constraint_rejections, 0L)
    expect_lt(attr(out, 'diagnostics')$constraint_acceptance, 1)
    expect_true('constraints' %in% attr(out, 'redist_attr'))
})

test_that('tree-weighted target agrees with an exact varying-weight oracle', {
    adj <- list(
        c(1L, 3L, 4L), c(0L, 2L, 4L), c(1L, 5L),
        c(0L, 4L), c(0L, 1L, 3L, 5L), c(2L, 4L)
    )
    map <- redist_map(pop = rep(1, 6L), ndists = 2L, pop_tol = 0.01, adj = adj)
    boundary_penalty <- 0.35

    is_connected <- function(vertices) {
        reached <- vertices[[1L]]
        frontier <- reached
        while (length(frontier)) {
            next_vertices <- unique(unlist(adj[frontier + 1L], use.names = FALSE))
            next_vertices <- intersect(next_vertices, vertices)
            next_vertices <- setdiff(next_vertices, reached)
            reached <- c(reached, next_vertices)
            frontier <- next_vertices
        }
        length(reached) == length(vertices)
    }
    tree_count <- function(vertices) {
        if (length(vertices) == 1L) {
            return(1)
        }
        graph <- matrix(0, length(vertices), length(vertices))
        for (i in seq_along(vertices)) {
            neighbors <- intersect(adj[[vertices[[i]] + 1L]], vertices)
            graph[i, match(neighbors, vertices)] <- 1
        }
        laplacian <- diag(rowSums(graph)) - graph
        round(det(laplacian[-1L, -1L, drop = FALSE]))
    }
    subsets <- combn(0:5, 3L, simplify = FALSE)
    subsets <- Filter(
        \(vertices) {
            0L %in% vertices && is_connected(vertices) && is_connected(setdiff(0:5, vertices))
        },
        subsets
    )
    oracle_plans <- lapply(subsets, \(vertices) {
        as.integer(!(0:5 %in% vertices)) + 1L
    })
    oracle_weights <- vapply(
        seq_along(subsets),
        \(i) {
            vertices <- subsets[[i]]
            plan <- oracle_plans[[i]]
            internal_edges <- sum(vapply(
                0:5,
                \(vertex) {
                    sum(
                        adj[[vertex + 1L]] > vertex &
                            plan[adj[[vertex + 1L]] + 1L] == plan[[vertex + 1L]]
                    )
                },
                integer(1)
            ))
            tree_count(vertices) *
                tree_count(setdiff(0:5, vertices)) *
                exp(boundary_penalty * internal_edges)
        },
        numeric(1)
    )
    names(oracle_weights) <- vapply(oracle_plans, plan_partition_key, character(1))
    oracle_probabilities <- oracle_weights / sum(oracle_weights)

    set.seed(2029)
    out <- redist_hier_dimer(
        map,
        4000L,
        boundary_penalty = boundary_penalty,
        init_plan = oracle_plans[[1L]],
        silent = TRUE
    )
    sampled_keys <- apply(get_plans_matrix(subset_sampled(out)), 2L, plan_partition_key)
    observed <- table(factor(sampled_keys, levels = names(oracle_probabilities))) /
        length(sampled_keys)

    expect_gt(min(observed), 0)
    expect_lt(sum(abs(as.numeric(observed) - oracle_probabilities)) / 2, 0.08)
})

test_that('hierarchy target is invariant to cut bias and tree refresh', {
    adj <- list(
        c(1L, 3L, 4L), c(0L, 2L, 4L), c(1L, 5L),
        c(0L, 4L), c(0L, 1L, 3L, 5L), c(2L, 4L)
    )
    map <- redist_map(pop = rep(1, 6L), ndists = 2L, pop_tol = 0.01, adj = adj)
    hierarchy <- matrix(c(1L, 1L, 1L, 2L, 2L, 2L), ncol = 1L)
    split_penalty <- 0.7
    boundary_penalty <- 0.25

    subsets <- combn(0:5, 3L, simplify = FALSE)
    subsets <- Filter(\(vertices) 0L %in% vertices, subsets)
    oracle_plans <- lapply(subsets, \(vertices) {
        as.integer(!(0:5 %in% vertices)) + 1L
    })
    oracle_weights <- vapply(
        seq_along(subsets),
        \(i) {
            vertices <- subsets[[i]]
            complement <- setdiff(0:5, vertices)
            first_trees <- hierarchical_tree_count(adj, vertices, hierarchy)
            second_trees <- hierarchical_tree_count(adj, complement, hierarchy)
            if (first_trees == 0 || second_trees == 0) {
                return(0)
            }
            plan <- oracle_plans[[i]]
            internal_edges <- sum(vapply(
                0:5,
                \(vertex) {
                    sum(
                        adj[[vertex + 1L]] > vertex &
                            plan[adj[[vertex + 1L]] + 1L] == plan[[vertex + 1L]]
                    )
                },
                integer(1)
            ))
            split_count <- sum(vapply(
                unique(hierarchy[, 1L]),
                \(unit) {
                    length(unique(plan[hierarchy[, 1L] == unit])) - 1L
                },
                integer(1)
            ))
            first_trees *
                second_trees *
                exp(boundary_penalty * internal_edges - split_penalty * split_count)
        },
        numeric(1)
    )
    supported <- oracle_weights > 0
    oracle_plans <- oracle_plans[supported]
    oracle_weights <- oracle_weights[supported]
    names(oracle_weights) <- vapply(oracle_plans, plan_partition_key, character(1))
    oracle_probabilities <- oracle_weights / sum(oracle_weights)

    tuning <- list(
        frequent_unbiased = list(cut_bias = 0, refresh = 1L, seed = 1301L),
        delayed_biased = list(cut_bias = 7, refresh = 4L, seed = 1302L)
    )
    sampled_probabilities <- lapply(tuning, \(config) {
        set.seed(config$seed)
        out <- redist_hier_dimer(
            map,
            6000L,
            hierarchy = list(region = hierarchy[, 1L]),
            split_penalty = split_penalty,
            cut_bias = config$cut_bias,
            boundary_penalty = boundary_penalty,
            init_plan = oracle_plans[[1L]],
            refresh = config$refresh,
            silent = TRUE
        )
        sampled_keys <- apply(get_plans_matrix(subset_sampled(out)), 2L, plan_partition_key)
        as.numeric(
            table(factor(sampled_keys, levels = names(oracle_probabilities))) /
                length(sampled_keys)
        )
    })
    total_variation <- vapply(
        sampled_probabilities,
        \(probabilities) {
            sum(abs(probabilities - oracle_probabilities)) / 2
        },
        numeric(1)
    )

    expect_gt(length(oracle_probabilities), 2L)
    expect_lt(max(total_variation), 0.08)
})

test_that('hierarchical split-dimer is reproducible and honors thinning', {
    adj <- list(
        c(1L, 3L), c(0L, 2L, 4L), c(1L, 5L),
        c(0L, 4L), c(1L, 3L, 5L), c(2L, 4L)
    )
    map <- redist_map(pop = rep(1, 6L), ndists = 3L, pop_tol = 0.01, adj = adj)
    init <- c(1L, 1L, 2L, 3L, 3L, 2L)

    set.seed(42)
    first <- redist_hier_dimer(map, 20L, thin = 2L, init_plan = init, silent = TRUE)
    set.seed(42)
    second <- redist_hier_dimer(map, 20L, thin = 2L, init_plan = init, silent = TRUE)
    set.seed(42)
    final <- redist_hier_dimer(
        map,
        20L,
        thin = 2L,
        init_plan = init,
        return_all = FALSE,
        silent = TRUE
    )

    expect_identical(get_plans_matrix(first), get_plans_matrix(second))
    expect_identical(first$mcmc_changed, second$mcmc_changed)
    sampled <- subset_sampled(first)
    saved <- cbind(init, get_plans_matrix(sampled))
    expected_changed <- vapply(
        seq_len(ncol(saved) - 1L),
        \(i) any(saved[, i] != saved[, i + 1L]),
        logical(1)
    )
    saved_changed <- sampled$mcmc_changed[
        seq.int(1L, length(sampled$mcmc_changed), by = 3L)
    ]
    expect_identical(saved_changed, expected_changed)
    expect_equal(ncol(get_plans_matrix(sampled)), 10L)
    expect_equal(ncol(get_plans_matrix(subset_sampled(final))), 1L)
})
