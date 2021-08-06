test_that("parity works", {
  out <- redist.parity(plans = plans_10[, 1:3],
                       total_pop = pop)

  expected <- c(0.0661665990642299, 0.0363453551413081, 0.0363453551413081)
  expect_equal(out, expected = expected)
})

test_that("parallel parity works", {
  skip_on_os('windows')
    out1 <- redist.parity(plans=plans_10[, 1:600], total_pop=pop, ncores=1)
    out2 <- redist.parity(plans=plans_10[, 1:600], total_pop=pop, ncores=2)

    expect_equal(out1, out2)
})
