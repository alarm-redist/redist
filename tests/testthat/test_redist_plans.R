test_that("constructor works", {
    fl = redist_map(fl25, ndists=3, pop_tol=0.1)
    x = redist_plans(plans_10, fl, "enumpart")
    expect_s3_class(x, "redist_plans")
    expect_s3_class(x, "tbl")
    expect_true(is.matrix(as.matrix(x)))
    expect_equal(attr(x, "prec_pop"), fl25$pop)
    expect_output(print(x), "drawn using Enumpart")
    expect_output(print(x), "927 sampled plans with 3 districts")
})


test_that("get_mh_acceptance_rate works", {
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
