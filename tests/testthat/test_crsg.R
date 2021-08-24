test_that("crsg works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  out <- redist.crsg(adj = adj, total_pop = fl25$pop, pop_tol = 0.1, ndists = 3,
                     shp = fl25, verbose = FALSE)
  out$district_list = lapply(out$district_list, sort)
  expected <- list(plan = c(1, 2, 2, 2, 3, 2, 3, 3, 1, 3, 3, 1, 1, 1, 1, 2, 
                            2, 2, 2, 2, 2, 3, 2, 3, 3),
                   district_list = list(sort(c(12L, 14L, 13L, 11L, 8L, 0L)),
                                        sort(c(18L, 19L, 1L, 16L, 15L, 17L, 20L, 3L, 22L, 2L, 5L)),
                                        sort(c(4L, 23L, 6L, 7L, 21L, 24L, 9L, 10L))),
                   district_pop = c(58683, 52653, 63707))
  expect_equal(out, expected)
})
