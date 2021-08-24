test_that("redist.competitiveness works", {
  out <- redist.competitiveness(plans_10[, 1:10], fl25$obama, fl25$mccain)
  expected <- c(0.104659637804475, 0.112209300916562, 0.1139660196912, 0.114984476024382, 
                0.113789960239524, 0.115931134105923, 0.111970617564315, 0.113200399522692, 
                0.117699805030613, 0.111751822732955)
  
  expect_equal(out, expected)
})
