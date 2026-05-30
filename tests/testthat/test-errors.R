# Argument-validation and edge-case tests for the sampling algorithms
# and validate_* helpers. These don't run full simulations — they just
# trigger the cli::cli_abort paths.

iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) |> suppressMessages()
n_iowa <- nrow(iowa_map)
ndists_iowa <- attr(iowa_map, "ndists")


# ---- redist_mergesplit ----

test_that("redist_mergesplit rejects invalid arguments", {
    expect_error(
        redist_mergesplit(iowa_map, nsims = 0, warmup = 0, silent = TRUE),
        "positive"
    )
    expect_error(
        redist_mergesplit(iowa_map, nsims = 10, compactness = -1, silent = TRUE),
        "non-negative"
    )
    expect_error(
        redist_mergesplit(
            iowa_map,
            nsims = 10,
            split_params = list(adapt_k_thresh = 1.5),
            silent = TRUE
        ),
        "0, 1"
    )
    expect_error(
        redist_mergesplit(iowa_map, nsims = 10, warmup = 0, thin = 0, silent = TRUE),
        "positive integer"
    )
})

test_that("redist_mergesplit rejects malformed init_plan", {
    expect_error(
        redist_mergesplit(iowa_map, nsims = 5, warmup = 0, init_plan = 1:3, silent = TRUE),
        "as many rows"
    )
    wrong_ndists <- rep(1L, n_iowa)
    expect_error(
        redist_mergesplit(
            iowa_map,
            nsims = 5,
            warmup = 0,
            init_plan = wrong_ndists,
            silent = TRUE
        ),
        "init_num_regions|districts"
    )
})


# ---- redist_flip ----

test_that("redist_flip rejects invalid arguments", {
    expect_error(
        redist_flip(iowa_map, nsims = 5, constraints = list(), silent = TRUE),
        "redist_constr"
    )
    expect_error(
        redist_flip(iowa_map, nsims = 5, thin = -1, silent = TRUE),
        "nonnegative"
    )
    expect_error(
        redist_flip(iowa_map, nsims = 5, init_plan = rep(1L, n_iowa), silent = TRUE),
        "same number of districts"
    )
})


# ---- redist_shortburst ----

test_that("redist_shortburst rejects invalid arguments", {
    expect_error(
        redist_shortburst(
            iowa_map,
            scorer_frac_kept(iowa_map),
            compactness = -1,
            max_bursts = 1,
            verbose = FALSE
        ),
        "non-negative"
    )
    expect_error(
        redist_shortburst(
            iowa_map,
            scorer_frac_kept(iowa_map),
            max_bursts = 0,
            verbose = FALSE
        ),
        "positive"
    )
    expect_error(
        redist_shortburst(
            iowa_map,
            scorer_frac_kept(iowa_map),
            constraints = list(),
            max_bursts = 1,
            verbose = FALSE
        ),
        "redist_constr"
    )
})


# ---- redist_cyclewalk ----

test_that("redist_cyclewalk rejects invalid arguments", {
    expect_error(redist_cyclewalk(iowa_map, nsims = 10, compactness = -1), "non-negative")
    expect_error(
        redist_cyclewalk(iowa_map, nsims = 5, warmup = 10),
        "thin"
    )
    expect_error(
        redist_cyclewalk(iowa_map, nsims = 50, thin = 1, cycle_walk_frac = 1.5),
        "between 0 and 1"
    )
    expect_error(
        redist_cyclewalk(iowa_map, nsims = 50, thin = 1, chains = 0),
        "positive"
    )
})


# ---- validate_redist_map ----

test_that("validate_redist_map flags non-redist_map and missing adjacency", {
    expect_error(validate_redist_map(list()), "redist_map")
    expect_error(validate_redist_map(iris), "redist_map")

    # a redist_map-classed object that's missing the adjacency column attribute
    no_adj <- iowa_map
    attr(no_adj, "adj_col") <- NULL
    expect_error(validate_redist_map(no_adj), "adjacency")
})


# ---- validate_edge_weights (cyclewalk) ----

test_that("validate_edge_weights rejects non-list edge_weights", {
    # exercised via redist_cyclewalk; passing a non-list edge_weights
    # bypasses the convenience numeric path and hits the validator.
    expect_error(
        redist_cyclewalk(iowa_map, nsims = 50, thin = 1, edge_weights = "bad"),
        "must be a list"
    )
})
