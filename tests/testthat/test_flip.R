test_that("flip works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

  capture.output(
  out <- redist_flip(fl_map, init_plan = plans_10[, 1],
                     nsims = 10, verbose = FALSE)
  )
  par <- redist.parity(get_plans_matrix(out), total_pop = pop)

  expect_equal(range(get_plans_matrix(out)), c(1, 3))
  expect_true(all(par <= 0.1))

})

test_that("flip countysplit works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")


  cty <- rep(1, 25)
  cty[1:4] <- 2

  cons <- redist_constr(fl_map) %>%
      add_constr_splits(
          strength = 10,
          admin = cty
      )

  capture.output(
  out <- redist_flip(fl_map, init_plan = plans_10[, 1],
                     nsims = 10, verbose = FALSE)
  )
  par <- redist.parity(get_plans_matrix(out), total_pop = pop)

  expect_equal(range(get_plans_matrix(out)), c(1, 3))
  expect_true(all(par <= 0.2))
})


test_that("flip hinge works", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

  cons <- redist_constr(fl_map) %>%
      add_constr_grp_hinge(
          strength = 5,
          group_pop = fl25$HispPop,
          tgts_group = c(0.4, 0.3, 0.2)
      )

  capture.output(
  out <- redist_flip(fl_map, init_plan = plans_10[, 1],
                     nsims = 10, verbose = FALSE,
                     constraints = cons)
  )
  par <- redist.parity(get_plans_matrix(out), total_pop = pop)

  expect_equal(range(get_plans_matrix(out)), c(1, 3))
  expect_true(all(par <= 0.1))
})


test_that('flip flip wrapper works',{
  capture.output(
  sims <- redist_flip(map = fl_map, nsims = 10)
  )

  expect_s3_class(sims, "redist_plans")
})


test_that('log-st works',{
  iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)

  cons <- redist_constr(iowa_map) %>%
      add_constr_log_st(strength = 1,
                        admin = region)

  capture.output(
  test <- redist_flip(iowa_map, nsims = 10,  constraints = cons)
  )

  expect_s3_class(test, 'data.frame')

})
