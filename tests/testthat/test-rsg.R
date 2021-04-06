test_that("rsg works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  out <- redist.rsg(adj = adj, total_pop = pop, pop_tol = 0.2, ndists = 3,
                    verbose = FALSE)

  par <- redist.parity(out$plan, total_pop = pop)

  expect_equal(range(out$plan), c(1,3))
  expect_true(par <= 0.2)
  expect_true(all(names(out) %in% c('plan', 'district_list', 'district_pop')))
})
