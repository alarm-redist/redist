test_that("set_pop_tol works", {
  data(iowa)
  iowa_map <- redist_map(iowa, existing_plan = cd_2010)
  iowa_map <- set_pop_tol(iowa_map, 0.05)
  
  tol <- get_pop_tol(iowa_map)
  
  expect_equal(tol, 0.05)
})
