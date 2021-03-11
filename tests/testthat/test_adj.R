context("Adjacency Lists")

test_that("Adjacency lists are reduced correctly", {
    g = list(c(1L, 2L), c(0L, 2L), c(0L, 1L, 3L), 2L)
    g_exp = list(c(1L, 2L), 0L, 0L)

    expect_equal(redist.reduce.adjacency(g, c(3, 1, 4)), g_exp)
})
