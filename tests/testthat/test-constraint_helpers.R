test_that("flip constraint helper works", {
  data(iowa)
  iowa_map <- redist_map(iowa, existing_plan = cd_2010)
  
  cons <- flip_constraints_helper(iowa_map)
  
  expect_named(cons)
  expect_true(class(cons) == 'list')
})
