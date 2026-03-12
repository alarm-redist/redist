# Skipped on CRAN and CI due to runtime

skip_on_cran()
skip_on_ci()

test_that('MMSS distribution matches expected (log st distribution)', {
  set.seed(123)

  result <- redist_mmss(grid,
    nsims = 100000, warmup = 1000, l = 3, thin = 10,
    init_plan = grid$init, compactness = 1, verbose = FALSE,
    exact_mh = TRUE
  )

  cut_edge_counts <- comp_edges_rem(result, shp = grid) |> by_plan(ndists = 4)
  observed_counts <- table(cut_edge_counts)
  total <- sum(observed_counts)
  observed <- sapply(c('8', '10', '11', '12'), function(k) {
    if (k %in% names(observed_counts)) observed_counts[[k]] / total else 0
  })

  # Expected: 8→256/654, 10→224/654, 11→96/654, 12→78/654
  expect_true(all(cut_edge_counts %in% c(8, 10, 11, 12)))
  expect_true(is_close(observed['8'], 256 / 654))
  expect_true(is_close(observed['10'], 224 / 654))
  expect_true(is_close(observed['11'], 96 / 654))
  expect_true(is_close(observed['12'], 78 / 654))
})

# test_that('MMSS distribution without spanning forest weighting (uniform over plans)', {
#   set.seed(456)
#
#   result <- redist_mmss(grid,
#     nsims = 100000, warmup = 1000, l = 3, thin = 1,
#     init_plan = grid$init, compactness = 0, verbose = FALSE
#   )
#
#   cut_edge_counts <- comp_edges_rem(result, shp = grid) |> by_plan(ndists = 4)
#   observed_counts <- table(cut_edge_counts)
#   total <- sum(observed_counts)
#   observed <- sapply(c('8', '10', '11', '12'), function(k) {
#     if (k %in% names(observed_counts)) observed_counts[[k]] / total else 0
#   })
#
#   # Expected: 8→1/117, 10→14/117, 11→24/117, 12→78/117
#   expect_true(all(cut_edge_counts %in% c(8, 10, 11, 12)))
#   expect_true(is_close(observed['8'], 1 / 117))
#   expect_true(is_close(observed['10'], 14 / 117))
#   expect_true(is_close(observed['11'], 24 / 117))
#   expect_true(is_close(observed['12'], 78 / 117))
# })
