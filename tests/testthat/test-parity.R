test_that("parity works", {
  data(fl25)
  data("algdat.p10")
  
  out <- redist.parity(plans = algdat.p10$cdmat[,1:3],
                       total_pop = fl25$pop)
  
  expected <- c(0.0661665990642299, 0.0363453551413081, 0.0363453551413081)
  expect_equal(out, expected = expected)
})
