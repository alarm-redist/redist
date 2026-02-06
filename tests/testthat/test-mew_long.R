# remove this chunk once merged in w/ cyclewalk ----

# 4x4 grid with unit population for distribution tests
grid_sf <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 4, ymax = 4)) |>
    sf::st_as_sfc() |>
    sf::st_make_grid(n = c(4, 4)) |>
    sf::st_sf(geometry = _)
grid_sf <- grid_sf[c(13, 14, 9, 10, 15, 16, 11, 12, 5, 6, 1, 2, 7, 8, 3, 4), ]
grid_sf$adj <- redist.adjacency(grid_sf)
grid_sf$pop <- rep(1L, 16)
grid_sf$init <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L)
grid <- redist_map(grid_sf, existing_plan = init, pop_bounds = c(3.99, 4, 4.01),
                   adj = grid_sf$adj) |> suppressWarnings()

# end setup ----

test_that('cycle walk distribution matches expected (gamma=0)', {
  set.seed(123)

  result <- redist_cyclewalk(grid,
    nsims = 2000000, warmup = 10000,
    init_plan = grid$init, compactness = 1, verbose = FALSE
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
