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

test_that("redist_mergesplit_parallel works", {
    skip_on_os("windows")
    data(fl25)
    fl_map = redist_map(fl25, ndists=3, pop_tol=0.1)

    N = 20
    chains = 3

    pl1 = redist_mergesplit_parallel(fl_map, nsims=N, chains=chains, ncores=2,
                                     silent=TRUE)
    pl2 = redist_mergesplit_parallel(fl_map, nsims=N, chains=chains, ncores=2,
                                     warmup=0, init_name=F, silent=TRUE)
    pl3 = redist_mergesplit_parallel(fl_map, nsims=N, chains=chains, ncores=2,
                                     return_all=F, init_name=F, silent=TRUE)

    expect_equal(get_n_ref(pl1), chains)
    expect_equal(get_n_ref(pl2), 0)
    expect_equal(get_n_ref(pl3), 0)

    expect_equal(nrow(pl1), 3*chains*(N/2+1))
    expect_equal(nrow(pl2), 3*chains*N)
    expect_equal(nrow(pl3), 3*chains)
})
