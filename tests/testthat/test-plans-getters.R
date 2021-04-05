test_that("get_mh_acceptance_rate works", {
  data('iowa')
  iowa_map <- redist_map(iowa, pop_tol = 0.01, existing_plan = cd_2010)
  
  # should exist and be number
  out <- redist_flip(iowa_map, nsims = 5, verbose = FALSE)
  mh <- get_mh_acceptance_rate(out)
  expect_true(is.numeric(mh))
  expect_true(mh >= 0)
  
  out <- redist_mergesplit(iowa_map, nsims = 5, silent = TRUE)
  mh <- get_mh_acceptance_rate(out)
  expect_true(is.numeric(mh))
  expect_true(mh >= 0)
  
  # should return NA_real_
  out <- redist_smc(iowa_map, nsims = 5, silent = TRUE)
  mh <- get_mh_acceptance_rate(out)
  expect_true(is.na(mh))
  expect_true(is.numeric(mh))
  
})
