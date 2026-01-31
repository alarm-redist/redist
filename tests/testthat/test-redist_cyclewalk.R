test_that("redist_cyclewalk runs on Iowa data", {
    result <- redist_cyclewalk(ia, nsims = 50)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 51)
    expect_equal(nrow(get_plans_matrix(result)), 99)
})

test_that("redist_cyclewalk respects population bounds", {
    set.seed(1)
    result <- redist_cyclewalk(ia, nsims = 20)

    pop_bounds <- attr(ia, "pop_bounds")
    pops <- result |>
        dplyr::group_by(draw) |>
        dplyr::summarize(
            min_pop = min(total_pop),
            max_pop = max(total_pop)
        )

    expect_true(all(pops$min_pop >= pop_bounds[1]))
    expect_true(all(pops$max_pop <= pop_bounds[3]))
})

test_that("redist_cyclewalk with verbose output works", {
    expect_output(
        redist_cyclewalk(ia, nsims = 10, verbose = TRUE),
        "Forest Walk"
    )
})

test_that("redist_cyclewalk with silent output works", {
    expect_silent(redist_cyclewalk(ia, nsims = 10, silent = TRUE))
})

test_that("cycle walk can produce different plans", {
    set.seed(1)
    result <- redist_cyclewalk(ia, nsims = 100, verbose = FALSE)

    plans_mat <- get_plans_matrix(result)

    # Count unique plans
    n_unique <- ncol(unique(plans_mat, MARGIN = 2))

    expect_true(n_unique >= 10)
    expect_true(ncol(plans_mat) == 101)
})

test_that("all generated plans are contiguous", {
    result <- redist_cyclewalk(ia, nsims = 50, verbose = FALSE)
    plans_mat <- get_plans_matrix(result)
    cont <- apply(plans_mat, 2, \(x) max(contiguity(adj = ia$adj, x)))
    expect_equal(max(cont), 1)
})

test_that("redist_cyclewalk runs with constraints", {
    skip_on_cran()

    constr <- redist_constr(ia) |>
        add_constr_pop_dev(strength = 10)

    result <- redist_cyclewalk(ia, nsims = 50, constraints = constr, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 51)

    constr <- redist_constr(ia) |>
        add_constr_splits(strength = 5, admin = region)

    result <- redist_cyclewalk(ia, nsims = 50, constraints = constr, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 51)
})

test_that("all plans have valid district structure", {
    result <- redist_cyclewalk(ia, nsims = 100, verbose = FALSE)
    plans_mat <- get_plans_matrix(result)

    ndists <- attr(ia, "ndists")

    expect_true(all(apply(plans_mat, 2, function(x) length(unique(x))) == ndists))
    expect_true(all(apply(plans_mat, 2, function(x) x >= 1 & x <= ndists)))
})

test_that("longer chain runs without crashing", {
    skip_on_cran()

    result <- redist_cyclewalk(ia, nsims = 500, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 501)
})

test_that("thinning works correctly", {
    # With thin=5 and nsims=100, should get 100/5 + 1 = 21 plans
    result <- redist_cyclewalk(ia, nsims = 100, thin = 5, verbose = FALSE)

    expect_equal(ncol(get_plans_matrix(result)), 21)
})

test_that("warmup burns initial samples", {
    result <- redist_cyclewalk(ia, nsims = 50, warmup = 10, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
    # 1 reference + 40 sampled
    expect_equal(ncol(get_plans_matrix(result)), 41)
})

test_that("init_plan parameter works", {
    result <- redist_cyclewalk(ia, nsims = 20, init_plan = ia$cd_2010, verbose = FALSE)

    plans_mat <- get_plans_matrix(result)

    expect_true(all(vctrs::vec_group_id(ia$cd_2010) == vctrs::vec_group_id(plans_mat[, 1])))
})

test_that("cycle walk produces valid plans on 4x4 grid", {
    skip_on_cran()

    grid <- create_4x4_grid()
    set.seed(42)

    result <- redist_cyclewalk(
        grid,
        nsims = 500,
        init_plan = grid$init,
        verbose = FALSE
    )

    cut_edge_counts <- comp_edges_rem(result, shp = grid)  |>
        by_plan(ndists = 4)

    # these are the only valid configurations
    expect_true(all(cut_edge_counts %in% c(8, 10, 11, 12)))

    # Should have some variation (not all the same plan)
    expect_true(length(unique(cut_edge_counts)) >= 2)
})
