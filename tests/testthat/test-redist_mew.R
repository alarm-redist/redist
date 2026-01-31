test_that('MEW basic functionality: structure, contiguity, population, diagnostics', {
  set.seed(42)
  
  # Run algorithm ONCE
  plans <- redist_mew(fl_map, 50, warmup = 10, silent = TRUE)
  plans_mat <- get_plans_matrix(plans)
  
  # Test output structure
  expect_s3_class(plans, 'redist_plans')
  expect_equal(nrow(plans_mat), 25) # 25 precincts
  expect_true(ncol(plans_mat) >= 35) # At least 35 draws after warmup
  
  # Test all districts exist
  expect_equal(range(as.matrix(plans)), c(1, 3))
  
  # Test contiguity
  contig <- apply(plans_mat, 2, function(x) all(contiguity(fl_map$adj, x) == 1))
  expect_true(all(contig), info = 'All MEW plans should be contiguous')
  
  # Test correct number of districts
  n_dists <- apply(plans_mat, 2, function(x) length(unique(x)))
  expect_true(all(n_dists == 3), info = 'All plans should have exactly 3 districts')
  
  # Test population tolerance
  pop_tol <- attr(fl_map, 'pop_tol')
  deviations <- redist.parity(plans_mat, fl25$pop)
  expect_true(all(deviations <= pop_tol),
    info = 'All plans should meet population constraint'
  )
  
  # Test population bounds
  pop_bounds <- attr(fl_map, 'pop_bounds')
  expect_true(all(plans$total_pop >= pop_bounds[1]))
  expect_true(all(plans$total_pop <= pop_bounds[3]))
  
  # Test diagnostics
  diag <- attr(plans, 'diagnostics')
  expect_true(!is.null(diag))
  expect_true('accept_rate' %in% names(diag))
  expect_true('cycle_intersect_rate' %in% names(diag))
  expect_true('avg_proposal_tries' %in% names(diag))
  expect_true(all(diag$accept_rate >= 0 & diag$accept_rate <= 1))
})

test_that('MEW warmup and thinning', {
  skip_on_cran()
  set.seed(456)
  
  # Test warmup
  plans_warmup <- redist_mew(fl_map, 60, warmup = 30, init_plan = plans_10[, 1], silent = TRUE)
  n_draws1 <- length(unique(plans_warmup$draw))
  expect_true(n_draws1 < 60)
  expect_true(n_draws1 >= 25)
  
  # Test thinning on same run
  plans_thin <- redist_mew(fl_map, 60, warmup = 10, thin = 4, init_plan = plans_10[, 1], silent = TRUE)
  n_draws2 <- length(unique(plans_thin$draw))
  expect_true(n_draws2 < n_draws1)
  expect_true(n_draws2 <= 15)
})

test_that('MEW works with different numbers of districts', {
  skip_on_cran()
  
  # Test with k=2
  set.seed(202)
  map2 <- redist_map(fl25, ndists = 2, pop_tol = 0.15) %>% suppressMessages()
  out2 <- redist_mew(map2, 30, warmup = 10, silent = TRUE)
  expect_equal(max(as.matrix(out2)), 2)

  # Test with k=4
  set.seed(100)
  map4 <- redist_map(fl25, ndists = 4, pop_tol = 0.15) %>% suppressMessages()
  out4 <- redist_mew(map4, 30, warmup = 10, silent = TRUE)
  expect_equal(max(as.matrix(out4)), 4)
})

test_that('MEW compactness parameter', {
  skip_on_cran()
  set.seed(303)

  # Run with low and high compactness
  out_low <- redist_mew(fl_map, 50, warmup = 10, compactness = 0.1, silent = TRUE)
  out_high <- redist_mew(fl_map, 50, warmup = 10, compactness = 5, silent = TRUE)

  # Both should produce valid plans
  expect_s3_class(out_low, 'redist_plans')
  expect_s3_class(out_high, 'redist_plans')

  # Higher compactness should lead to fewer cut edges
  cuts_low <- redistmetrics::comp_edges_rem(out_low, fl_map)
  cuts_high <- redistmetrics::comp_edges_rem(out_high, fl_map)
  expect_true(median(cuts_high) < median(cuts_low))
})

test_that('MEW input validation', {
  expect_error(redist_mew(fl_map, 10, compactness = -1, warmup = 0), 'non-negative')
  expect_error(redist_mew(fl_map, 10, nsims = 5, warmup = 10), 'greater than')
  expect_error(redist_mew(fl_map, 100, thin = 0, warmup = 0), 'positive integer')
  expect_error(redist_mew(fl_map, 100, thin = 150, warmup = 0), 'no larger than')
  
  # Invalid init_plan length
  expect_error(
    redist_mew(fl_map, 10, warmup = 0, init_plan = c(1, 2, 3)),
    'as long as the number of units'
  )

  # Wrong number of districts
  bad_plan <- rep(1:2, length.out = 25)
  expect_error(
    redist_mew(fl_map, 10, warmup = 0, init_plan = bad_plan),
    'same number of districts'
  )

  # Non-contiguous plan
  bad_plan <- c(1, 1, 1, 2, 2, 3, 2, 2, 1, 1, 3, 2, 2, 1, 1, 3, 2, 2, 1, 1, 3, 2, 2, 1, 1)
  expect_error(
    redist_mew(fl_map, 10, warmup = 0, init_plan = bad_plan),
    'must be contiguous'
  )
})

test_that('MEW with NULL population tolerance', {
  set.seed(404)
  map_no_tol <- fl_map
  attr(map_no_tol, 'pop_tol') <- NULL

  out <- redist_mew(map_no_tol, 30, warmup = 10, silent = TRUE)
  expect_s3_class(out, 'redist_plans')
})

test_that('MEW reproducibility with set.seed', {
  set.seed(505)
  pl1 <- redist_mew(fl_map, 50, warmup = 10, silent = TRUE)
  set.seed(505)
  pl2 <- redist_mew(fl_map, 50, warmup = 10, silent = TRUE)

  expect_identical(as.matrix(pl1), as.matrix(pl2))
})

test_that('MEW on small grids', {
  skip_on_cran()
  set.seed(606)

  bb <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0)))))
  grid <- sf::st_make_grid(bb, n = 4)
  tb <- sf::st_as_sf(grid) %>%
    rename(geometry = x) %>%
    mutate(pop = c(10, 15, 12, 18, 20, 14, 11, 16, 13, 17, 19, 12, 15, 18, 14, 16))

  grid_map <- redist_map(tb, ndists = 4, pop_tol = 0.25, total_pop = pop) %>%
    suppressWarnings() %>%
    suppressMessages()

  init_plan <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)

  out <- redist_mew(grid_map, 30, warmup = 10, init_plan = init_plan, silent = TRUE)
  expect_s3_class(out, 'redist_plans')
  expect_equal(max(as.matrix(out)), 4)
})

test_that('MEW on Iowa data', {
  skip_on_cran()
  set.seed(909)

  iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.01) %>%
    suppressMessages()

  plans <- redist_mew(iowa_map, 100, warmup = 20, silent = TRUE)

  pop_bounds <- attr(iowa_map, 'pop_bounds')
  expect_true(all(plans$total_pop >= pop_bounds[1]))
  expect_true(all(plans$total_pop <= pop_bounds[3]))
})

test_that('MEW acceptance rate', {
  skip_on_cran()
  set.seed(808)
  out <- redist_mew(fl_map, 100, warmup = 20, silent = TRUE)

  diag <- attr(out, 'diagnostics')
  avg_accept <- mean(diag$accept_rate)

  expect_true(avg_accept >= 0.3)
  expect_true(avg_accept <= 1.0)
})

test_that('MEW integration test with larger sample', {
  skip_on_cran()
  set.seed(999)

  plans <- redist_mew(fl_map, 500, warmup = 50, silent = TRUE)

  # Check output structure
  expect_s3_class(plans, 'redist_plans')
  n_draws <- length(unique(plans$draw))
  expect_true(n_draws == 451)

  # Check diagnostics
  diag <- attr(plans, 'diagnostics')
  expect_true(!is.null(diag))
  avg_accept <- mean(diag$accept_rate, na.rm = TRUE)
  expect_true(is.numeric(avg_accept) && length(avg_accept) == 1)
  expect_true(avg_accept >= 0 && avg_accept <= 1.0)

  # Verify all plans valid
  plans_mat <- get_plans_matrix(plans)
  expect_true(all(apply(plans_mat, 2, function(x) length(unique(x)) == 3)))

  pop_bounds <- attr(fl_map, 'pop_bounds')
  expect_true(all(plans$total_pop >= pop_bounds[1]))
  expect_true(all(plans$total_pop <= pop_bounds[3]))
})

test_that('MEW parallel chains', {
  skip_on_cran()
  set.seed(1010)

  plans <- redist_mew(fl_map, 50, warmup = 10, chains = 2, ncores = 2, silent = TRUE)

  # Check output structure
  expect_s3_class(plans, 'redist_plans')
  expect_true('chain' %in% names(plans))

  # Check 2 chains
  expect_equal(length(unique(plans$chain)), 2)

  # Check equal samples
  chain_counts <- table(plans$chain)
  expect_equal(as.numeric(chain_counts[1]), as.numeric(chain_counts[2]))

  # Check diagnostics
  diag <- attr(plans, 'diagnostics')
  expect_true(is.list(diag))
  expect_equal(length(diag), 2)

  expect_true(all(sapply(diag, function(x) 'accept_rate' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'cycle_intersect_rate' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'runtime' %in% names(x))))

  # All plans valid
  plans_mat <- get_plans_matrix(plans)
  expect_true(all(apply(plans_mat, 2, function(x) length(unique(x)) == 3)))

  # Population constraints met
  pop_bounds <- attr(fl_map, 'pop_bounds')
  expect_true(all(plans$total_pop >= pop_bounds[1]))
  expect_true(all(plans$total_pop <= pop_bounds[3]))
})
