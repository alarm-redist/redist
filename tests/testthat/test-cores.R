# Tests for cores/coarsening helpers: redist.identify.cores,
# redist.coarsen.adjacency, redist.uncoarsen, redist.reduce.adjacency.

test_that("redist.identify.cores returns one core id per precinct", {
    iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) |>         suppressMessages()
    core <- suppressWarnings(
        redist.identify.cores(adj = get_adj(iowa_map), plan = iowa_map$cd_2010)
    )
    expect_length(core, nrow(iowa_map))
    expect_true(all(is.finite(core)))
    expect_true(min(core) >= 1)

    # `simplify = FALSE` returns the raw core list with extra components
    core_full <- suppressWarnings(redist.identify.cores(
        adj = get_adj(iowa_map),
        plan = iowa_map$cd_2010,
        simplify = FALSE
    ))
    expect_type(core_full, "list")
    expect_true(all(c("dm", "k", "conncomp", "group_id") %in% names(core_full)))
    expect_equal(core_full$dm, iowa_map$cd_2010)
    expect_equal(core_full$group_id, core)
})

test_that("redist.identify.cores accepts a plan matrix and uses the first column", {
    iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) |>         suppressMessages()
    mat <- cbind(iowa_map$cd_2010, iowa_map$cd_2010)
    core1 <- suppressWarnings(redist.identify.cores(
        adj = get_adj(iowa_map),
        plan = mat
    ))
    core2 <- suppressWarnings(redist.identify.cores(
        adj = get_adj(iowa_map),
        plan = iowa_map$cd_2010
    ))
    expect_equal(core1, core2)
})

test_that("redist.identify.cores validates adj input", {
    expect_error(redist.identify.cores(adj = "not-a-list", plan = 1:3), "list")
})

test_that("redist.coarsen.adjacency groups precincts and shrinks adjacency", {
    adj0 <- list(c(1L, 2L), c(0L, 2L), c(0L, 1L, 3L), 2L) # 0-indexed
    # group the first two precincts together
    groups <- c(0L, 0L, 1L, 2L)
    coarse <- redist.coarsen.adjacency(adj0, groups)

    expect_type(coarse, "list")
    expect_length(coarse, length(unique(groups)))
    # every neighbor index must be in [0, length-1] and not self-referential
    for (i in seq_along(coarse)) {
        expect_true(all(coarse[[i]] %in% (seq_along(coarse) - 1L)))
        expect_false((i - 1L) %in% coarse[[i]])
    }
})

test_that("redist.coarsen.adjacency validates inputs", {
    adj0 <- list(c(1L, 2L), c(0L, 2L), c(0L, 1L, 3L), 2L)

    # length mismatch
    expect_error(redist.coarsen.adjacency(adj0, groups = 1:3), "different sizes")

    # 1-indexed adj should fail the zero-index check
    adj1 <- lapply(adj0, function(x) x + 1L)
    expect_error(redist.coarsen.adjacency(adj1, c(0L, 0L, 1L, 2L)), "0-indexed")
})

test_that("redist.uncoarsen reverses a coarsening on a plan matrix", {
    # 4 precincts collapsed into 3 groups (precincts 1+2 share a group)
    groups <- c(1L, 1L, 2L, 3L)
    coarse_plans <- matrix(
        c(1L, 2L, 1L,
                              1L, 2L, 2L),
        ncol = 2
    )
    expanded <- redist.uncoarsen(coarse_plans, group_index = groups)
    expect_equal(dim(expanded), c(length(groups), ncol(coarse_plans)))
    # precincts 1 and 2 (same group) get the same district in every column
    expect_equal(expanded[1, ], expanded[2, ])
    # column 1: groups 1->1, 2->2, 3->1
    expect_equal(expanded[, 1], c(1L, 1L, 2L, 1L))
    expect_equal(expanded[, 2], c(1L, 1L, 2L, 2L))
})

test_that("redist.reduce.adjacency subsets to kept rows and reindexes", {
    adj <- list(c(1L, 2L), c(0L, 2L), c(0L, 1L, 3L), 2L)
    keep <- c(1L, 3L) # 1-indexed in the original
    reduced <- redist.reduce.adjacency(adj, keep_rows = keep)

    expect_length(reduced, length(keep))
    # neighbor indices must be 0..length(keep)-1 and exclude self
    for (i in seq_along(reduced)) {
        expect_true(all(reduced[[i]] %in% (seq_along(reduced) - 1L)))
        expect_false((i - 1L) %in% reduced[[i]])
    }
})
