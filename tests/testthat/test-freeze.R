test_that("redist.freeze works", {
  data(fl25)
  adj <- redist.adjacency(fl25)
  out <- redist.freeze(adj = adj, plan = algdat.p10$cdmat[, 1], freeze_row =  (2 == algdat.p10$cdmat[, 1]))
  expected <- c(1L, 2L, 2L, 2L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 
                2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)
  
  expect_equal(out, expected)
  
  
  out <- redist.freeze(adj = adj, plan = algdat.p10$cdmat[, 1], 
                       freeze_row =  (2 == algdat.p10$cdmat[, 1])| (3 == algdat.p10$cdmat[, 1]))
  expected <- c(1L, 2L, 2L, 2L, 2L, 3L, 4L, 5L, 6L, 5L, 5L, 5L, 5L, 5L, 5L, 
                2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)
  
  expect_equal(out, expected)
  
  out <- redist.freeze(adj = adj, 
                       freeze_row =  (2 == algdat.p10$cdmat[, 1])| (3 == algdat.p10$cdmat[, 1]))
  
  
  })

test_that("freeze works", {
  data(fl25)
  data("algdat.p10")
  fl25$plan <- algdat.p10$cdmat[,1]
  st_crs(fl25) <- 4269
  flmap <- redist_map(fl25, existing_plan = plan)
  out <- flmap %>% freeze(plan == 2, .data = .)
  
  expected <- c(1L, 2L, 2L, 2L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 
                2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)
  
  expect_equal(out, expected)
})
