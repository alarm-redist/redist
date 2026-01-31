# Tests for redist_cyclewalk edge weights
# fl_map already created in setup.R

test_that('redist_cyclewalk works with NULL edge_weights', {
  set.seed(02139)
  result <- redist_cyclewalk(fl_map, nsims = 20)

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 21)
  expect_equal(nrow(get_plans_matrix(result)), 25)
})

test_that('redist_cyclewalk works with single edge weight', {
  ew <- list(list(edge = c(1, 2), weight = 2.0))

  set.seed(02139)
  result <- redist_cyclewalk(fl_map,
    nsims = 20, edge_weights = ew
  )

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 21)
  expect_equal(nrow(get_plans_matrix(result)), 25)
})

test_that('redist_cyclewalk works with multiple edge weights', {
  ew <- list(
    list(edge = c(1, 2), weight = 2.0),
    list(edge = c(1, 3), weight = 3.0),
    list(edge = c(1, 4), weight = 1.5)
  )

  set.seed(02139)
  result <- redist_cyclewalk(fl_map,
    nsims = 20, edge_weights = ew
  )

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 21)
  expect_equal(nrow(get_plans_matrix(result)), 25)
})

test_that('edge weights can be specified in either direction', {
  skip_on_cran()

  ew1 <- list(list(edge = c(1, 2), weight = 2.0))
  ew2 <- list(list(edge = c(2, 1), weight = 2.0))

  set.seed(02139)
  result1 <- redist_cyclewalk(fl_map,
    nsims = 10, edge_weights = ew1
  )

  set.seed(02139)
  result2 <- redist_cyclewalk(fl_map,
    nsims = 10, edge_weights = ew2
  )

  expect_s3_class(result1, 'redist_plans')
  expect_s3_class(result2, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result1)), 11)
  expect_equal(ncol(get_plans_matrix(result2)), 11)
})

test_that('edge weights prints custom weight count when verbose', {
  skip_on_cran()

  ew <- list(
    list(edge = c(1, 2), weight = 2.0),
    list(edge = c(1, 3), weight = 3.0)
  )

  expect_output(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew,
      verbose = TRUE, silent = FALSE
    ),
    'Using 2 custom edge weights'
  )
})

test_that('edge_weights must be a list', {
  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = 'not a list'
    ),
    'must be a list'
  )
})

test_that('edge_weights entries must be lists', {
  ew <- list('not a list')

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'Entry 1.*must be a list'
  )
})

test_that("edge_weights entries must have 'edge' field", {
  ew <- list(list(weight = 2.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'missing.*edge.*field'
  )
})

test_that("edge_weights entries must have 'weight' field", {
  ew <- list(list(edge = c(1, 2)))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'missing.*weight.*field'
  )
})

test_that('edge field must be numeric vector of length 2', {
  ew <- list(list(edge = c(1, 2, 3), weight = 2.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'must be a numeric vector of length 2'
  )

  ew <- list(list(edge = 'not numeric', weight = 2.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'must be a numeric vector of length 2'
  )
})

test_that('vertices must be in range', {
  ew <- list(list(edge = c(1, 999), weight = 2.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'out of range'
  )

  ew <- list(list(edge = c(0, 2), weight = 2.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'out of range'
  )

  ew <- list(list(edge = c(-1, 2), weight = 2.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'out of range'
  )
})

test_that('edge must exist in adjacency graph', {
  ew <- list(list(edge = c(2, 5), weight = 2.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'not in adjacency graph'
  )
})

test_that('adjacency check works correctly with 0-indexed adj list', {
  ew1 <- list(list(edge = c(1, 2), weight = 2.0))
  ew2 <- list(list(edge = c(2, 1), weight = 2.0))

  expect_silent(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew1, silent = TRUE
    )
  )

  expect_silent(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew2, silent = TRUE
    )
  )
})

test_that('weight must be positive', {
  ew <- list(list(edge = c(1, 2), weight = 0.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'weight must be positive'
  )

  ew <- list(list(edge = c(1, 2), weight = -1.0))

  expect_error(
    redist_cyclewalk(fl_map,
      nsims = 5, edge_weights = ew
    ),
    'weight must be positive'
  )
})

test_that('edge weights work with different compactness values', {
  skip_on_cran()

  ew <- list(list(edge = c(1, 2), weight = 5.0))

  set.seed(02139)
  result1 <- redist_cyclewalk(fl_map,
    nsims = 20, edge_weights = ew,
    compactness = 0.5
  )

  set.seed(02139)
  result2 <- redist_cyclewalk(fl_map,
    nsims = 20, edge_weights = ew,
    compactness = 2.0
  )

  expect_s3_class(result1, 'redist_plans')
  expect_s3_class(result2, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result1)), 21)
  expect_equal(ncol(get_plans_matrix(result2)), 21)
})

test_that('edge weights work with constraints', {
  skip_on_cran()

  ew <- list(list(edge = c(1, 2), weight = 5.0))

  constr <- redist_constr(fl_map) |>
    add_constr_pop_dev(strength = 10)

  result <- redist_cyclewalk(fl_map,
    nsims = 20, edge_weights = ew,
    constraints = constr
  )

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 21)
  expect_equal(nrow(get_plans_matrix(result)), 25)
})

test_that('all plans remain contiguous with edge weights', {
  skip_on_cran()

  ew <- list(
    list(edge = c(1, 2), weight = 10.0),
    list(edge = c(1, 3), weight = 10.0)
  )

  result <- redist_cyclewalk(fl_map,
    nsims = 50, edge_weights = ew
  )

  plans_mat <- get_plans_matrix(result)
  expect_equal(ncol(plans_mat), 51)

  cont <- apply(plans_mat, 2, \(x) max(contiguity(adj = fl_map$adj, x)))
  expect_equal(max(cont), 1)
  expect_equal(min(cont), 1)
})

test_that('all plans respect population bounds with edge weights', {
  skip_on_cran()

  ew <- list(list(edge = c(1, 2), weight = 5.0))

  set.seed(1)
  result <- redist_cyclewalk(fl_map,
    nsims = 50, edge_weights = ew
  )

  expect_equal(ncol(get_plans_matrix(result)), 51)

  pop_bounds <- attr(fl_map, 'pop_bounds')
  pops <- result |>
    dplyr::group_by(draw) |>
    dplyr::summarize(
      min_pop = min(total_pop),
      max_pop = max(total_pop)
    )

  expect_true(all(pops$min_pop >= pop_bounds[1]))
  expect_true(all(pops$max_pop <= pop_bounds[3]))
  expect_equal(nrow(pops), 51)
})

test_that('empty edge_weights list works', {
  skip_on_cran()

  ew <- list()

  result <- redist_cyclewalk(fl_map,
    nsims = 10, edge_weights = ew
  )

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 11)
})

test_that('fractional weights work', {
  skip_on_cran()

  ew <- list(
    list(edge = c(1, 2), weight = 0.5),
    list(edge = c(1, 3), weight = 2.5)
  )

  result <- redist_cyclewalk(fl_map,
    nsims = 20, edge_weights = ew
  )

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 21)
})

test_that('very large weights work', {
  skip_on_cran()

  ew <- list(list(edge = c(1, 2), weight = 1000000.0))

  result <- redist_cyclewalk(fl_map,
    nsims = 10, edge_weights = ew
  )

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 11)
})

test_that('very small weights work', {
  skip_on_cran()

  ew <- list(list(edge = c(1, 2), weight = 0.001))

  result <- redist_cyclewalk(fl_map,
    nsims = 10, edge_weights = ew
  )

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 11)
})

# County-based edge weights tests (uses ia from setup_cyclewalk.R)

test_that('counties parameter auto-generates edge weights', {
  skip_on_cran()

  set.seed(123)
  result <- redist_cyclewalk(ia, 20, counties = iowa$region, silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 21)
})

test_that('edge_weights as number sets county weight multiplier', {
  skip_on_cran()

  set.seed(123)
  result <- redist_cyclewalk(ia, 20, counties = iowa$region, edge_weights = 5, silent = TRUE)
  expect_s3_class(result, 'redist_plans')

  set.seed(123)
  result2 <- redist_cyclewalk(ia, 20, counties = iowa$region, edge_weights = 20, silent = TRUE)
  expect_s3_class(result2, 'redist_plans')
})

test_that('edge_weights as number requires counties', {
  expect_error(
    redist_cyclewalk(fl_map, nsims = 5, edge_weights = 5),
    'requires.*counties'
  )
})

test_that('explicit edge_weights list overrides county auto-generation', {
  skip_on_cran()

  custom_weights <- list(list(edge = c(1, ia$adj[[1]][1] + 1), weight = 3.0))

  set.seed(123)
  result <- redist_cyclewalk(ia, 20,
    counties = iowa$region,
    edge_weights = custom_weights, silent = TRUE
  )
  expect_s3_class(result, 'redist_plans')
})
