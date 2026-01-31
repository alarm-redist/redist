# Long-running distribution tests for redist_cyclewalk
# These tests verify that cycle walk samples from the correct distribution
# Ported from Julia CycleWalk.jl test suite
#
# These tests are skipped by default. To run them:
#   testthat::test_file("tests/testthat/test_redist_cyclewalk_long.R")
#
# Or from command line:
#   Rscript -e "devtools::test(filter='cyclewalk_long')"

library(testthat)
library(redist)

skip_on_cran()
skip_on_ci()

# setup ----

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

# check if observed is close to expected
# For large probabilities (>1%), use 10% relative tolerance
# For small probabilities, use 40% relative tolerance
is_close <- function(obs, exp) {
    if (exp > 0.01) {
        ratio <- obs / exp
        ratio >= 0.9 && ratio <= 1.1
    } else {
        ratio <- obs / exp
        ratio >= 0.6 && ratio <= 1.4
    }
}

# ============================================================
# Distribution tests
# ============================================================

test_that("cycle walk distribution matches expected (unweighted, gamma=0)", {
    skip_on_cran()
    grid <- create_4x4_grid()
    set.seed(123)

    # Julia test setup:
    # - cycle_steps=200_000 outer iterations
    # - instep=10 (10 MCMC steps per outer loop)
    # - Total: 2 million MCMC steps, 200K recorded samples
    #
    # Our approach: use thin=10 to match Julia's thinning
    # nsims=200000, thin=10 means we run 200K*10 = 2M MCMC steps
    # and record 200K/10 = 20K samples (after thinning)
    #
    # Actually thin works differently - it keeps every 10th sample from nsims
    # So nsims=200000, thin=10 -> ~20K final samples from 200K steps
    #
    # To get 2M steps with 200K samples, we need nsims that accounts for thinning
    # Let's use nsims=2000000 (2M steps) with thin=10 for 200K samples
    result <- redist_cyclewalk(
        grid,
        nsims = 2000000,
        thin = 10,
        warmup = 10000,
        init_plan = grid$init,
        verbose = FALSE
    )

    # Count cut edges for each plan
    cut_edge_counts <- comp_edges_rem(result, shp = grid)  |>
        by_plan(ndists = 4)

    # Expected distribution (from Julia tests, gamma=0, unweighted):
    # These counts represent the number of valid partitions with each cut edge count
    # multiplied by their symmetry counts
    # 8 cut edges: 256/654 ≈ 0.391
    # 10 cut edges: (128+64+32)/654 = 224/654 ≈ 0.343
    # 11 cut edges: 96/654 ≈ 0.147
    # 12 cut edges: 78/654 ≈ 0.119
    expected <- c(
        `8` = 256/654,
        `10` = 224/654,
        `11` = 96/654,
        `12` = 78/654
    )

    # Compute observed distribution
    observed_counts <- table(cut_edge_counts)
    total <- sum(observed_counts)
    observed <- sapply(c("8", "10", "11", "12"), function(k) {
        if (k %in% names(observed_counts)) observed_counts[[k]]/total else 0
    })

    cat("\nObserved distribution:\n")
    print(round(observed, 4))
    cat("\nExpected distribution:\n")
    print(round(expected, 4))

    # Check that we only have valid cut edge counts
    expect_true(all(cut_edge_counts %in% c(8, 10, 11, 12)),
    info = paste("Invalid cut edge counts found:",
    paste(unique(cut_edge_counts[!cut_edge_counts %in% c(8, 10, 11, 12)]),
    collapse = ", ")))

    # Check each probability matches expected
    for (ce in names(expected)) {
        expect_true(
            is_close(observed[ce], expected[ce]),
            info = paste("Cut edges =", ce,
            ": observed =", round(observed[ce], 4),
            ", expected =", round(expected[ce], 4),
            ", ratio =", round(observed[ce]/expected[ce], 3))
        )
    }
})

test_that("cycle walk distribution with spanning forest weighting (gamma=1)", {
    skip_on_cran()

    grid <- create_4x4_grid()
    set.seed(456)

    # Create constraint with spanning forest weighting
    # This should change the distribution to weight by 1/spanning_trees
    constr <- redist_constr(grid) |>
        add_constr_log_st(strength = 1.0)

    # Julia uses cycle_steps=100_000 with instep=10 = 1M MCMC steps
    # We use nsims=1000000, thin=10 for similar behavior
    result <- redist_cyclewalk(
        grid,
        nsims = 1000000,
        thin = 10,
        warmup = 10000,
        init_plan = grid$init,
        constraints = constr,
        verbose = FALSE
    )

    cut_edge_counts <- comp_edges_rem(result, shp = grid)  |>
        by_plan(ndists = 4)

    # Expected distribution with gamma=1 (from Julia tests):
    # Now weighted by spanning_trees^(1-gamma) = spanning_trees^0 = 1
    # So it's just counting unique partitions by symmetry
    # 8 cut edges: 1/117 (just 1 partition type with 256 spanning trees)
    # 10 cut edges: 14/117 (3 types * varying symmetries)
    # 11 cut edges: 24/117
    # 12 cut edges: 78/117
    expected <- c(
        `8` = 1/117,
        `10` = 14/117,
        `11` = 24/117,
        `12` = 78/117
    )

    # Compute observed distribution
    observed_counts <- table(cut_edge_counts)
    total <- sum(observed_counts)
    observed <- sapply(c("8", "10", "11", "12"), function(k) {
        if (k %in% names(observed_counts)) observed_counts[[k]]/total else 0
    })

    cat("\nObserved distribution (gamma=1):\n")
    print(round(observed, 4))
    cat("\nExpected distribution (gamma=1):\n")
    print(round(expected, 4))

    # For now, just check we get valid cut edge counts
    # The distribution may differ from Julia due to algorithm differences
    expect_true(all(cut_edge_counts %in% c(8, 10, 11, 12)))

    # Check each probability matches expected
    for (ce in names(expected)) {
        expect_true(
            is_close(observed[ce], expected[ce]),
            info = paste("Cut edges =", ce,
                         ": observed =", round(observed[ce], 4),
                         ", expected =", round(expected[ce], 4),
                         ", ratio =", round(observed[ce]/expected[ce], 3))
        )
    }
})

test_that("longer chain produces stable distribution", {
    skip_on_cran()

    grid <- create_4x4_grid()
    set.seed(789)

    # Run 2M MCMC steps like Julia test (2M steps with thin=10 = 200K samples)
    result <- redist_cyclewalk(
        grid,
        nsims = 2000000,
        thin = 10,
        warmup = 10000,
        init_plan = grid$init,
        verbose = FALSE
    )

    cut_edge_counts <- comp_edges_rem(result, shp = grid)  |>
        by_plan(ndists = 4)

    expected <- c(
        `8` = 256/654,
        `10` = 224/654,
        `11` = 96/654,
        `12` = 78/654
    )

    # Compute observed distribution
    observed_counts <- table(cut_edge_counts)
    total <- sum(observed_counts)
    observed <- sapply(c("8", "10", "11", "12"), function(k) {
        if (k %in% names(observed_counts)) observed_counts[[k]]/total else 0
    })

    cat("\nObserved distribution (long chain):\n")
    print(round(observed, 4))
    cat("\nExpected distribution:\n")
    print(round(expected, 4))

    # With longer chain, should have tighter tolerance
    is_close_tight <- function(obs, exp) {
        ratio <- obs / exp
        ratio >= 0.95 && ratio <= 1.05
    }

    # Check each probability matches expected (tighter tolerance)
    for (ce in names(expected)) {
        expect_true(
            is_close_tight(observed[ce], expected[ce]),
            info = paste("Cut edges =", ce,
            ": observed =", round(observed[ce], 4),
            ", expected =", round(expected[ce], 4),
            ", ratio =", round(observed[ce]/expected[ce], 3))
        )
    }
})
