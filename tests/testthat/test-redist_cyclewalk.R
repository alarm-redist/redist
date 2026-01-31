# Tests for redist_cyclewalk

library(testthat)
library(redist)

# load Iowa test data
data(iowa)
ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)

test_that("redist_cyclewalk runs on Iowa data", {
    skip_on_cran()

    result <- redist_cyclewalk(ia, nsims = 50)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 51)
    expect_equal(nrow(get_plans_matrix(result)), 99)
})

test_that("redist_cyclewalk respects population bounds", {
    skip_on_cran()

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
    skip_on_cran()

    expect_output(
        redist_cyclewalk(ia, nsims = 10, verbose = TRUE),
        "Forest Walk"
    )
})

test_that("redist_cyclewalk with silent output works", {
    skip_on_cran()

    expect_silent(redist_cyclewalk(ia, nsims = 10, silent = TRUE))
})

test_that("partition initializes correctly", {
    skip_on_cran()


    result <- redist_cyclewalk(ia, nsims = 10)

    ndists <- attr(ia, "ndists")
    expect_equal(ndists, 4)

    pops <- result |>
        dplyr::filter(draw == "1") |>
        dplyr::pull(total_pop)

    expect_equal(length(pops), 4)
    expect_true(all(pops > 0))
})

test_that("cycle walk can produce different plans", {
    skip_on_cran()


    # Run with enough iterations for cycle walk to produce new plans
    result <- redist_cyclewalk(ia, nsims = 100, verbose = FALSE)

    plans_mat <- get_plans_matrix(result)

    # Count unique plans
    n_unique <- ncol(unique(plans_mat, MARGIN = 2))

    # With 10% cycle walk rate and 100 iterations, we should get multiple unique plans
    # (approximately 10 cycle walk attempts, most should succeed and produce changes)
    # Requiring at least 2 unique plans ensures the algorithm is actually working
    expect_true(n_unique >= 2,
                info = paste("Expected at least 2 unique plans but got", n_unique,
                            "out of", ncol(plans_mat), "total plans"))
    expect_true(ncol(plans_mat) > 1)
})

test_that("all generated plans are contiguous", {
    skip_on_cran()

    result <- redist_cyclewalk(ia, nsims = 50, verbose = FALSE)
    plans_mat <- get_plans_matrix(result)
    cont <- apply(plans_mat, 2, \(x) max(contiguity(adj = ia$adj, x)))
    expect_equal(max(cont), 1)
})

test_that("redist_cyclewalk works with pop_dev constraint", {
    skip_on_cran()


    # Create constraint object with population deviation penalty
    constr <- redist_constr(ia) |>
        add_constr_pop_dev(strength = 10)

    result <- redist_cyclewalk(ia, nsims = 50, constraints = constr, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 51)
})

test_that("redist_cyclewalk works with splits constraint", {
    skip_on_cran()


    # Create constraint object with county splits penalty
    constr <- redist_constr(ia) |>
        add_constr_splits(strength = 5, admin = region)

    result <- redist_cyclewalk(ia, nsims = 50, constraints = constr, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 51)
})

# ============================================================
# Iteration 6: Validation tests
# ============================================================

test_that("all plans have valid district structure", {
    skip_on_cran()


    result <- redist_cyclewalk(ia, nsims = 100, verbose = FALSE)
    plans_mat <- get_plans_matrix(result)

    # Check each plan has valid structure
    ndists <- attr(ia, "ndists")

    expect_true(all(apply(plans_mat, 2, function(x) length(unique(x))) == ndists))
    expect_true(all(apply(plans_mat, 2, function(x) x >= 1 & x <= ndists)))
})

test_that("longer chain runs without crashing", {
    skip_on_cran()


    # Run a longer chain to test stability
    result <- redist_cyclewalk(ia, nsims = 500, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
    expect_equal(ncol(get_plans_matrix(result)), 501)
})

test_that("thinning works correctly", {
    skip_on_cran()


    # With thin=5 and nsims=100, should get 100/5 + 1 = 21 plans
    result <- redist_cyclewalk(ia, nsims = 100, thin = 5, verbose = FALSE)

    # Reference plan + 100/5 = 21 plans
    expect_equal(ncol(get_plans_matrix(result)), 21)
})

test_that("warmup burns initial samples", {
    skip_on_cran()


    # With warmup, the first plan shouldn't necessarily be the init
    result <- redist_cyclewalk(ia, nsims = 50, warmup = 10, verbose = FALSE)

    expect_s3_class(result, "redist_plans")
     # 1 reference + 40 sampled
    expect_equal(ncol(get_plans_matrix(result)), 41)
})

test_that("init_plan parameter works", {
    skip_on_cran()

    result <- redist_cyclewalk(ia, nsims = 20, init_plan = ia$cd_2010, verbose = FALSE)

    plans_mat <- get_plans_matrix(result)

    expect_true(all(vctrs::vec_group_id(ia$cd_2010) == vctrs::vec_group_id(plans_mat[, 1])))
})

# make 4x4 grid graph with unit population
create_4x4_grid <- function() {
    grid_sf <- st_bbox(c(xmin = 0, ymin = 0, xmax = 4, ymax = 4)) |>
        sf::st_as_sfc() |>
        sf::st_make_grid(n = c(4, 4)) |>
        sf::st_sf(geometry = _)

    new_order <- c(13, 14, 9, 10, 15, 16, 11, 12, 5, 6, 1, 2, 7, 8, 3, 4)
    grid_sf <- grid_sf[new_order, ]

    grid_sf$adj <- redist.adjacency(grid_sf)
    grid_sf$pop <- rep(1L, 16)
    grid_sf$init <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L)
    grid_sf$id <- 1:16

    # make pop = 4 only possible choice
    redist_map(
        grid_sf,
        existing_plan = init,
        pop_bounds = c(3.99, 4, 4.01),
        adj = grid_sf$adj
    ) |>
      suppressWarnings()
}

test_that("cycle walk produces valid plans on 4x4 grid", {
    skip_on_cran()

    grid <- create_4x4_grid()
    set.seed(42)

    # Run a moderate number of iterations
    result <- redist_cyclewalk(
        grid,
        nsims = 500,
        init_plan = grid$init,
        verbose = FALSE
    )

    # Count cut edges for each plan
    cut_edge_counts <- comp_edges_rem(result, shp = grid)  |>
        by_plan(ndists = 4)

    # All plans should have 8, 10, 11, or 12 cut edges
    # (these are the only valid configurations for 4x4 grid with 4 districts)
    expect_true(all(cut_edge_counts %in% c(8, 10, 11, 12)),
                info = paste("Invalid cut edge counts found:",
                            paste(unique(cut_edge_counts[!cut_edge_counts %in% c(8, 10, 11, 12)]),
                                  collapse = ", ")))

    # Should have some variation (not all the same plan)
    expect_true(length(unique(cut_edge_counts)) >= 2,
                info = "Expected multiple different cut edge counts")
})
