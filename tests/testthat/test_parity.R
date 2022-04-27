test_that("parity works", {
  out <- redist.parity(plans = plans_10[, 1:3],
                       total_pop = pop)

  expected <- c(0.0661665990642299, 0.0363453551413081, 0.0363453551413081)
  expect_equal(out, expected = expected)
})

