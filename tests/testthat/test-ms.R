test_that("redist.mergesplit works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

  nsims <- 10
  out <- redist.mergesplit(adj = adj, total_pop = pop, init_plan = plans_10[, 1],
                     nsims = nsims, ndists = 3, verbose = FALSE, pop_tol = 0.1)
  par <- redist.parity(out$plans, total_pop = pop)

  expect_equal(out$nsims, nsims)
  expect_equal(range(out$plans), c(1, 3))
  expect_true(all(par <= 0.1))
})
