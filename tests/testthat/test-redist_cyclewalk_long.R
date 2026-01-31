# Long-running distribution tests for redist_cyclewalk
# Verifies cycle walk samples from the correct distribution (ported from Julia CycleWalk.jl)
# Skipped on CRAN and CI due to runtime

skip_on_cran()
skip_on_ci()

test_that('cycle walk distribution matches expected (gamma=0)', {
  set.seed(123)

  # gamma=0 in Julia = compactness=0 in R (uniform over spanning forests)
  result <- redist_cyclewalk(grid,
    nsims = 200000, instep = 10, warmup = 1000,
    init_plan = grid$init, compactness = 0, verbose = FALSE
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

test_that('cycle walk distribution with spanning forest weighting (gamma=1)', {
  set.seed(456)

  # gamma=1 in Julia = compactness=1 in R (uniform over partitions)
  result <- redist_cyclewalk(grid,
    nsims = 100000, instep = 10, warmup = 1000,
    init_plan = grid$init, compactness = 1, verbose = FALSE
  )

  cut_edge_counts <- comp_edges_rem(result, shp = grid) |> by_plan(ndists = 4)
  observed_counts <- table(cut_edge_counts)
  total <- sum(observed_counts)
  observed <- sapply(c('8', '10', '11', '12'), function(k) {
    if (k %in% names(observed_counts)) observed_counts[[k]] / total else 0
  })

  # Expected: 8→1/117, 10→14/117, 11→24/117, 12→78/117
  expect_true(all(cut_edge_counts %in% c(8, 10, 11, 12)))
  expect_true(is_close(observed['8'], 1 / 117))
  expect_true(is_close(observed['10'], 14 / 117))
  expect_true(is_close(observed['11'], 24 / 117))
  expect_true(is_close(observed['12'], 78 / 117))
})

test_that('longer chain produces stable distribution', {
  set.seed(789)

  result <- redist_cyclewalk(grid,
    nsims = 200000, instep = 10, warmup = 1000,
    init_plan = grid$init, compactness = 0, verbose = FALSE
  )

  cut_edge_counts <- comp_edges_rem(result, shp = grid) |> by_plan(ndists = 4)
  observed_counts <- table(cut_edge_counts)
  total <- sum(observed_counts)
  observed <- sapply(c('8', '10', '11', '12'), function(k) {
    if (k %in% names(observed_counts)) observed_counts[[k]] / total else 0
  })

  # Tighter tolerance for longer chain
  expect_true(is_close_tight(observed['8'], 256 / 654))
  expect_true(is_close_tight(observed['10'], 224 / 654))
  expect_true(is_close_tight(observed['11'], 96 / 654))
  expect_true(is_close_tight(observed['12'], 78 / 654))
})
