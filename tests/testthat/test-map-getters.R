test_that("get_pop_tol works", {
  data(iowa)
  iowa_map <- redist_map(iowa, existing_plan = cd_2010)
  tol <- get_pop_tol(iowa_map)

  expect_equal(tol, 0.000053506567)
})
