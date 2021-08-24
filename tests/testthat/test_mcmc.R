test_that("mcmc works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

  nsims <- 10
  capture.output(
  out <- redist.flip(adj = adj, total_pop = pop, init_plan = plans_10[, 1],
                     nsims = nsims, ndists = 3, verbose = FALSE, pop_tol = 0.1)
  )
  par <- redist.parity(out$plans, total_pop = pop)

  expect_equal(out$nsims, nsims)
  expect_equal(range(out$plans), c(1, 3))
  expect_true(all(par <= 0.1))

})

test_that("mcmc countysplit works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")


  cty <- rep(1, 25)
  cty[1:4] <- 2

  nsims <- 10
  capture.output(
  out <- redist.flip(adj = adj, total_pop = pop, init_plan = plans_10[, 1],
                     nsims = nsims, ndists = 3, verbose = FALSE, pop_tol = 0.2,
                     counties = cty, constraint = 'countysplit', constraintweights = 5)
  )
  par <- redist.parity(out$plans, total_pop = pop)

  expect_equal(out$nsims, nsims)
  expect_equal(range(out$plans), c(1, 3))
  expect_true(all(par <= 0.2))
})


test_that("mcmc hinge works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

  nsims <- 10
  capture.output(
  out <- redist.flip(adj = adj, total_pop = pop, init_plan = plans_10[, 1],
                     nsims = nsims, ndists = 3, verbose = FALSE, pop_tol = 0.1,
                     constraint = 'hinge', constraintweights = 5,
                     group_pop = fl25$HispPop, minorityprop = c(0.4, 0.3, 0.2))
  )
  par <- redist.parity(out$plans, total_pop = pop)

  expect_equal(out$nsims, nsims)
  expect_equal(range(out$plans), c(1, 3))
  expect_true(all(par <= 0.1))
})


test_that('mcmc flip wrapper works',{
  data(fl25)
  nsims <- 10

  capture.output(
  sims <- redist_flip(map = fl_map, nsims = 10)
  )

  expect_s3_class(sims, "redist_plans")
})


test_that('log-st works',{
  iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)

  cons <- flip_constraints_helper(map = iowa_map, constraintweight = 1,
                                  compactness_metric = 'log-st', counties = name)

  capture.output(
  test <- redist_flip(iowa_map, nsims = 10,  constraints = cons)
  )

  expect_s3_class(test, 'data.frame')

})
