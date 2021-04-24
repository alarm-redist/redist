test_that("Adjacency lists are reduced correctly", {
    g = list(c(1L, 2L), c(0L, 2L), c(0L, 1L, 3L), 2L)
    g_exp = list(c(1L, 2L), 0L, 0L)

    expect_equal(redist.reduce.adjacency(g, c(3, 1, 4)), g_exp)
})

test_that("redist.adjacency works", {
  adj <- redist.adjacency(fl25)
  expected <- list(1:14, c(0L, 3L, 15L, 16L, 17L, 18L, 19L, 20L, 21L), 
                c(0L, 3L, 5L, 20L, 22L), c(0L, 1L, 2L, 20L), c(0L, 6L, 21L, 23L, 24L), 
                c(0L, 2L, 8L, 22L), c(0L, 4L, 7L), c(0L, 6L, 9L), c(0L, 5L, 13L), 
                c(0L, 7L, 10L), c(0L, 9L, 11L), c(0L, 10L, 12L), 
                c(0L, 11L, 14L), c(0L, 8L, 14L), c(0L, 12L, 13L), c(1L, 16L, 17L), 
                c(1L, 15L, 18L), c(1L, 15L, 20L), c(1L, 16L, 19L), 
                c(1L, 18L, 21L), c(1L, 2L, 3L, 17L, 22L), c(1L, 4L, 19L, 23L, 24L
              ), c(2L, 5L, 20L), c(4L, 21L, 24L), c(4L, 21L, 23L))
  expect_equal(adj, expected)
})
