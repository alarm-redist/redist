# Tests for various utility helpers that don't fit elsewhere.

test_that("redist.reorder renumbers districts to start at 1 in order of first appearance", {
    plans <- matrix(
        c(4L, 5L, 2L, 1L, 3L,
                       5L, 4L, 3L, 2L, 1L),
        ncol = 2
    )
    out <- redist.reorder(plans)
    expect_true(is.matrix(out))
    expect_equal(dim(out), dim(plans))
    # the first precinct is always district 1 after renumbering
    expect_equal(out[1, ], c(1L, 1L))
    # each column is a relabeling of the input column that preserves the partition
    for (j in seq_len(ncol(plans))) {
        rel <- match(plans[, j], unique(plans[, j]))
        expect_equal(out[, j], rel)
    }
})

test_that("redist.reorder rejects non-numeric input", {
    expect_error(redist.reorder("oops"), "matrix or integer vector")
})


test_that("redist.random.subgraph returns a connected subset of the requested size", {
    set.seed(1)
    sub <- redist.random.subgraph(iowa, n = 10)
    expect_named(sub, c("shp", "keep_rows"))
    expect_equal(nrow(sub$shp), 10)
    expect_length(sub$keep_rows, 10)
    expect_true(all(sub$keep_rows >= 1 & sub$keep_rows <= nrow(iowa)))
    expect_equal(anyDuplicated(sub$keep_rows), 0L)
    # the induced subgraph should be connected
    sub_adj <- redist.reduce.adjacency(redist.adjacency(iowa), sub$keep_rows)
    expect_true(all(contiguity(sub_adj, rep(1L, length(sub_adj))) == 1))
})

test_that("redist.random.subgraph rejects bad n", {
    expect_error(redist.random.subgraph(iowa, n = 0), "positive")
    expect_error(redist.random.subgraph(iowa, n = nrow(iowa) + 1), "smaller")
})


test_that("redist.subset returns shp + adj for kept rows and computes sub_pop_tol", {
    set.seed(2)
    keep <- 1:25
    sub <- redist.subset(
        shp = iowa,
        keep_rows = keep,
        total_pop = iowa$pop,
        ndists = 4,
        pop_tol = 0.05,
        sub_ndists = 1
    )
    expect_named(sub, c("shp", "adj", "keep_rows", "sub_ndists", "sub_pop_tol"))
    expect_length(keep, nrow(sub$shp))
    expect_length(sub$adj, length(keep))
    expect_equal(sub$sub_ndists, 1)
    expect_true(is.finite(sub$sub_pop_tol))

    # without all of total_pop/ndists/pop_tol/sub_ndists, the tol fields are NA
    sub2 <- redist.subset(shp = iowa, keep_rows = keep)
    expect_true(is.na(sub2$sub_pop_tol))
    expect_true(is.na(sub2$sub_ndists))
})

test_that("redist.subset enforces required arguments", {
    expect_error(redist.subset(), "shp")
    expect_error(redist.subset(shp = iris), "sf")
})


test_that("freeze can hold multiple districts in place", {
    data(fl25_adj)
    plan <- plans_10[, 1]
    # freeze districts 2 and 3 together
    out <- redist.freeze(
        adj = fl25_adj,
        plan = plan,
        freeze_row = plan == 2 | plan == 3
    )
    expect_length(out, length(plan))
    # within each frozen district, all units collapse to a single id
    expect_length(unique(out[plan == 2]), 1)
    expect_length(unique(out[plan == 3]), 1)
})
