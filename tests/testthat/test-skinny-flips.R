test_that("skinny flips works", {
  
  cons <- process_flip_constr(list(), rep(0, 25), rep(1,25))
  
  plans <- skinny_flips(adj = adj, init_plan = plans_10[, 1],
               total_pop = fl25$pop, pop_tol = 0.15,
               nsims = 5, eprob = 0.05, lambda = 1, constraints = cons)
  
  
  expect_true(class(plans)[1] ==  "matrix" & class(plans)[2] == "array" )
  expect_true(min(plans) == 1)
  
})
