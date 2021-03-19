test_that("mcmc works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  nsims <- 10
  out <- redist.mcmc(adj = adj, total_pop = pop, init_plan = algdat.p10$cdmat[,1],
                     nsims = nsims, ndists = 3, verbose = FALSE, pop_tol = 0.1)
  par <- redist.parity(out$plans, total_pop = pop)
  
  expect_equal(out$nsims, nsims)
  expect_equal(range(out$plans), c(0, 2))
  expect_true(all(par <= 0.1))

})

test_that("mcmc countysplit works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  
  cty <- rep(1, 25)
  cty[1:4] <- 2
  
  nsims <- 10
  out <- redist.mcmc(adj = adj, total_pop = pop, init_plan = algdat.p10$cdmat[,1],
                     nsims = nsims, ndists = 3, verbose = FALSE, pop_tol = 0.1,
                     counties = cty, constraint = 'countysplit', constraintweights = 5)
  par <- redist.parity(out$plans, total_pop = pop)
  
  expect_equal(out$nsims, nsims)
  expect_equal(range(out$plans), c(0, 2))
  expect_true(all(par <= 0.1))
  
})


test_that("mcmc hinge works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  

  nsims <- 10
  out <- redist.mcmc(adj = adj, total_pop = pop, init_plan = algdat.p10$cdmat[,1],
                     nsims = nsims, ndists = 3, verbose = FALSE, pop_tol = 0.1,
                     constraint = 'hinge', constraintweights = 5,
                     group_pop = fl25$HispPop, minorityprop = c(0.4, 0.3, 0.2))
  par <- redist.parity(out$plans, total_pop = pop)
  
  expect_equal(out$nsims, nsims)
  expect_equal(range(out$plans), c(0, 2))
  expect_true(all(par <= 0.1))
  
})