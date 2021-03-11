data(fl25)
adj <- redist.adjacency(fl25)
test_that("rsg works", {
  set.seed(1)
  out <- redist.rsg(adj = adj, total_pop = fl25$pop, pop_tol = 0.1, ndists = 3,
                    verbose = FALSE)
  expected <- list(plan = c(0L, 0L, 0L, 0L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 
                2L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L), 
       district_list = list(
                  c(17L, 15L, 16L, 0L, 3L, 20L, 2L, 19L, 18L, 1L, 22L), 
                  c(6L, 4L, 24L, 23L, 21L, 7L, 9L), c(14L, 12L, 11L, 10L, 13L, 8L, 5L)), 
       district_pop = c(59739, 57131, 58173))
  expect_equal(out, expected)
})
