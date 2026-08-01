test_that('exact matching partition functions agree with brute force', {
    cases <- list(
        cycle4 = list(
            n = 4L,
            u = c(1L, 2L, 3L, 4L),
            v = c(2L, 3L, 4L, 1L),
            w = c(2, 3, 5, 7)
        ),
        complete4 = list(
            n = 4L,
            u = c(1L, 1L, 1L, 2L, 2L, 3L),
            v = c(2L, 3L, 4L, 3L, 4L, 4L),
            w = c(2, 3, 5, 7, 11, 13)
        ),
        cycle6 = list(
            n = 6L,
            u = c(1L, 2L, 3L, 4L, 5L, 6L),
            v = c(2L, 3L, 4L, 5L, 6L, 1L),
            w = c(2, 3, 5, 7, 11, 13)
        ),
        disconnected = list(
            n = 6L,
            u = c(1L, 2L, 3L, 4L, 5L),
            v = c(2L, 3L, 4L, 1L, 6L),
            w = c(2, 3, 5, 7, 11)
        ),
        parallel = list(
            n = 2L,
            u = c(1L, 1L),
            v = c(2L, 2L),
            w = c(1, 3)
        )
    )

    oracle <- lapply(cases, function(case) {
        oracle <- matching_oracle(case$n, case$u, case$v, case$w)
        sum(vapply(oracle, `[[`, numeric(1), 'weight'))
    })
    outputs <- lapply(cases, function(case) {
        hdimer_matching_diagnostic(case$n, case$u, case$v, case$w)
    })

    expect_identical(
        unname(vapply(outputs, `[[`, logical(1), 'success')),
        rep(TRUE, length(cases))
    )
    expect_identical(
        unname(vapply(outputs, `[[`, character(1), 'status')),
        rep('success', length(cases))
    )
    expect_identical(
        unname(vapply(outputs, `[[`, logical(1), 'planarity_checked')),
        rep(TRUE, length(cases))
    )
    expect_identical(
        unname(vapply(outputs, `[[`, logical(1), 'planar')),
        rep(TRUE, length(cases))
    )
    expect_identical(
        unname(vapply(outputs, `[[`, logical(1), 'numerical_failure')),
        rep(FALSE, length(cases))
    )
    expect_identical(
        unname(vapply(outputs, `[[`, logical(1), 'used_subset_dp')),
        rep(TRUE, length(cases))
    )
    expect_equal(
        vapply(outputs, `[[`, numeric(1), 'log_partition'),
        log(unlist(oracle)),
        tolerance = 1e-10
    )
    expect_lt(max(vapply(outputs, `[[`, numeric(1), 'max_probability_error')), 1e-8)
})

test_that('larger matching graphs retain the polynomial Kasteleyn path', {
    n_vertices <- 16L
    edge_u <- seq_len(n_vertices)
    edge_v <- c(2:n_vertices, 1L)
    out <- hdimer_matching_diagnostic(n_vertices, edge_u, edge_v, rep(1, n_vertices))

    expect_identical(out$success, TRUE)
    expect_identical(out$planar, TRUE)
    expect_identical(out$used_subset_dp, FALSE)
    expect_equal(out$log_partition, log(2), tolerance = 1e-10)
})

test_that('the Kasteleyn backend samples the exact weighted law', {
    columns <- seq_len(8L)
    edge_u <- c(columns, seq_len(7L), 8L + seq_len(7L))
    edge_v <- c(8L + columns, 2:8, 8L + 2:8)
    weights <- seq(0.8, 2.1, length.out = length(edge_u))
    oracle <- matching_oracle(16L, edge_u, edge_v, weights)
    oracle_weights <- vapply(oracle, `[[`, numeric(1), 'weight')
    expected_edge_probability <- sum(
        oracle_weights[vapply(oracle, \(matching) 1L %in% matching$edges, logical(1))]
    ) /
        sum(oracle_weights)

    partition <- hdimer_matching_diagnostic(16L, edge_u, edge_v, weights)
    set.seed(7301)
    sampled_edge <- replicate(
        2000L,
        1L %in%
            hdimer_matching_diagnostic(
                16L,
                edge_u,
                edge_v,
                weights
            )$sampled_edge_indices
    )

    expect_identical(partition$success, TRUE)
    expect_identical(partition$used_subset_dp, FALSE)
    expect_equal(partition$log_partition, log(sum(oracle_weights)), tolerance = 1e-10)
    expect_lt(abs(mean(sampled_edge) - expected_edge_probability), 0.04)
})

test_that('perfect-matching samples have the exact weighted law', {
    edge_u <- c(1L, 2L, 3L, 4L)
    edge_v <- c(2L, 3L, 4L, 1L)
    weights <- c(2, 3, 5, 7)

    set.seed(117)
    draws <- replicate(1500, {
        sort(hdimer_matching_diagnostic(4L, edge_u, edge_v, weights)$sampled_edge_indices)
    })
    first_matching <- colMeans(draws == c(1L, 3L)) == 1
    observed <- mean(first_matching)
    expected <- weights[[1L]] *
        weights[[3L]] /
        (weights[[1L]] * weights[[3L]] + weights[[2L]] * weights[[4L]])

    expect_lt(abs(observed - expected), 0.04)

    set.seed(918)
    parallel_draws <- replicate(1200, {
        hdimer_matching_diagnostic(
            2L,
            c(1L, 1L),
            c(2L, 2L),
            c(1, 3)
        )$sampled_edge_indices
    })
    expect_lt(abs(mean(parallel_draws == 2L) - 0.75), 0.04)
})

test_that('matching failures retain correct planarity diagnostics', {
    isolated <- hdimer_matching_diagnostic(2L, integer(), integer(), numeric())
    star <- hdimer_matching_diagnostic(
        4L,
        c(1L, 1L, 1L),
        c(2L, 3L, 4L),
        rep(1, 3)
    )
    k33 <- hdimer_matching_diagnostic(
        6L,
        rep(1:3, each = 3L),
        rep(4:6, times = 3L),
        rep(1, 9)
    )

    expect_identical(isolated$success, FALSE)
    expect_identical(isolated$planar, TRUE)
    expect_identical(isolated$numerical_failure, FALSE)
    expect_identical(star$success, FALSE)
    expect_identical(star$planar, TRUE)
    expect_identical(star$numerical_failure, FALSE)
    expect_identical(k33$success, FALSE)
    expect_identical(k33$planar, FALSE)
})

test_that('matching probabilities are invariant to a global weight scale', {
    edge_u <- c(1L, 2L, 3L, 4L)
    edge_v <- c(2L, 3L, 4L, 1L)
    weights <- c(2, 3, 5, 7)
    base <- hdimer_matching_diagnostic(4L, edge_u, edge_v, weights)
    scaled <- hdimer_matching_diagnostic(4L, edge_u, edge_v, weights * 1e200)

    expect_identical(scaled$success, TRUE)
    expect_identical(scaled$numerical_failure, FALSE)
    expect_equal(
        scaled$log_partition - base$log_partition,
        2 * log(1e200),
        tolerance = 1e-10
    )
})

test_that('matching diagnostics distinguish unchecked and numerical failures', {
    invalid <- hdimer_matching_diagnostic(2L, 1L, 3L, 1)
    numerical <- hdimer_matching_diagnostic(
        2L,
        c(1L, 1L),
        c(2L, 2L),
        c(.Machine$double.xmin, .Machine$double.xmax)
    )

    expect_identical(invalid$status, 'invalid_input')
    expect_identical(invalid$planarity_checked, FALSE)
    expect_identical(invalid$planar, NA)
    expect_identical(numerical$status, 'numerical_failure')
    expect_identical(numerical$success, FALSE)
    expect_identical(numerical$numerical_failure, TRUE)
})
