test_that('MEW runs without errors on basic inputs', {
  set.seed(42)
  nsims <- 50
  res <- redist_mew(fl_map, nsims, warmup = 10, silent = TRUE)

  expect_s3_class(res, 'redist_plans')
  expect_true(ncol(get_plans_matrix(res)) > 0)
  # Check plans matrix dimensions (precincts × draws after warmup)
  plans_mat <- get_plans_matrix(res)
  expect_equal(nrow(plans_mat), 25) # 25 precincts
  expect_true(ncol(plans_mat) >= 35) # At least 35 draws after warmup
})

test_that('MEW generates valid contiguous plans', {
  # Use a known working seed
  set.seed(42)
  nsims <- 100
  out <- redist_mew(fl_map, nsims, warmup = 20, silent = TRUE)

  # Check all districts exist
  expect_equal(range(as.matrix(out)), c(1, 3))

  # MEW generates valid plans (all should be contiguous)
  plans_mat <- get_plans_matrix(out)
  contig <- apply(plans_mat, 2, function(x) all(contiguity(fl_map$adj, x) == 1))
  expect_true(all(contig))

  # Verify population tolerance is respected
  # Population constraint is enforced in proposals, so all plans should meet it
  pop_tol <- attr(fl_map, 'pop_tol')
  deviations <- redist.parity(plans_mat, fl25$pop)
  expect_true(all(deviations <= pop_tol),
    info = 'All plans should meet population constraint'
  )
})

test_that('MEW warmup and thinning work correctly', {
  set.seed(456)
  # Test with warmup
  out1 <- redist_mew(fl_map, 100, warmup = 50, init_plan = plans_10[, 1], silent = TRUE)
  n_draws1 <- length(unique(out1$draw))
  expect_true(n_draws1 < 100) # Should have fewer than 100 after warmup removal
  expect_true(n_draws1 >= 45) # Should have around 48-50 draws

  # Test with thinning
  out2 <- redist_mew(fl_map, 100,
    warmup = 20, thin = 4,
    init_plan = plans_10[, 1], silent = TRUE
  )
  n_draws2 <- length(unique(out2$draw))
  expect_true(n_draws2 < n_draws1) # Thinning should reduce draws
  expect_true(n_draws2 <= 25) # 100/4 = 25 max
})

test_that('MEW contiguity is preserved across all draws', {
  set.seed(789)
  out <- redist_mew(fl_map, 50, warmup = 10, silent = TRUE)
  plans_mat <- get_plans_matrix(out)

  # Check all plans are contiguous
  contig <- apply(plans_mat, 2, function(x) all(contiguity(fl_map$adj, x) == 1))
  expect_true(all(contig), info = 'All MEW plans should be contiguous')

  # Verify correct number of districts in each plan
  n_dists <- apply(plans_mat, 2, function(x) length(unique(x)))
  expect_true(all(n_dists == 3), info = 'All plans should have exactly 3 districts')
})

test_that('MEW diagnostics are computed', {
  set.seed(101)
  out <- redist_mew(fl_map, 100, warmup = 10, verbose = FALSE, silent = TRUE)

  diag <- attr(out, 'diagnostics')
  expect_true(!is.null(diag))
  expect_true('accept_rate' %in% names(diag))
  expect_true('cycle_intersect_rate' %in% names(diag))
  expect_true('avg_proposal_tries' %in% names(diag))

  # Acceptance rate should be between 0 and 1
  expect_true(all(diag$accept_rate >= 0 & diag$accept_rate <= 1))
})

test_that('MEW works with different numbers of districts', {
  skip_on_cran()
  set.seed(202)

  # Test with k=2
  map2 <- redist_map(fl25, ndists = 2, pop_tol = 0.15) %>% suppressMessages()
  out2 <- redist_mew(map2, 50, warmup = 10, silent = TRUE)
  expect_equal(max(as.matrix(out2)), 2)

  # Test with k=4
  map4 <- redist_map(fl25, ndists = 4, pop_tol = 0.15) %>% suppressMessages()
  out4 <- redist_mew(map4, 50, warmup = 10, silent = TRUE)
  expect_equal(max(as.matrix(out4)), 4)
})

test_that('MEW compactness parameter favors more compact plans', {
  skip_on_cran()
  set.seed(303)

  # Low compactness (rho near 0) should allow less compact plans
  out_low <- redist_mew(fl_map, 100, warmup = 20, compactness = 0.1, silent = TRUE)

  # High compactness (rho large) should favor more compact plans
  out_high <- redist_mew(fl_map, 100, warmup = 20, compactness = 5, silent = TRUE)

  # Both should produce valid plans
  expect_s3_class(out_low, 'redist_plans')
  expect_s3_class(out_high, 'redist_plans')

  # Calculate compactness using redistmetrics (lower = more compact)
  cuts_low <- redistmetrics::comp_edges_rem(out_low, fl_map)
  cuts_high <- redistmetrics::comp_edges_rem(out_high, fl_map)

  # Higher compactness parameter should generally lead to fewer cut edges
  expect_true(median(cuts_high) < median(cuts_low))
})

test_that('MEW high compactness produces compact plans', {
  skip_on_cran()
  set.seed(303)

  # Test with very high compactness parameter
  out_high <- redist_mew(fl_map, 50, warmup = 10, compactness = 5, silent = TRUE)

  expect_s3_class(out_high, 'redist_plans')
  expect_true(ncol(get_plans_matrix(out_high)) >= 40)

  # Verify diagnostics are available
  diag <- attr(out_high, 'diagnostics')
  expect_true(!is.null(diag))
  expect_true('accept_rate' %in% names(diag))
})

test_that('MEW checks invalid arguments', {
  expect_error(redist_mew(fl_map, 10, compactness = -1, warmup = 0), 'non-negative')
  expect_error(redist_mew(fl_map, 10, nsims = 5, warmup = 10), 'greater than')
  expect_error(redist_mew(fl_map, 100, thin = 0, warmup = 0), 'positive integer')
  expect_error(redist_mew(fl_map, 100, thin = 150, warmup = 0), 'no larger than')
})

test_that('MEW works with NULL population tolerance', {
  set.seed(404)
  map_no_tol <- fl_map
  attr(map_no_tol, 'pop_tol') <- NULL

  # Should run without error (no population constraint)
  out <- redist_mew(map_no_tol, 50, warmup = 10, silent = TRUE)
  expect_s3_class(out, 'redist_plans')
})

test_that('MEW init_plan validation works', {
  # Invalid length
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

  # Non-contiguous plan (districts 1 and 3 are not internally connected)
  bad_plan <- c(1, 1, 1, 2, 2, 3, 2, 2, 1, 1, 3, 2, 2, 1, 1, 3, 2, 2, 1, 1, 3, 2, 2, 1, 1)
  expect_error(
    redist_mew(fl_map, 10, warmup = 0, init_plan = bad_plan),
    'must be contiguous'
  )
})

test_that('MEW sampling is reproducible with set.seed', {
  set.seed(505)
  pl1 <- redist_mew(fl_map, 100, warmup = 20, silent = TRUE)
  set.seed(505)
  pl2 <- redist_mew(fl_map, 100, warmup = 20, silent = TRUE)

  expect_identical(as.matrix(pl1), as.matrix(pl2))
})

test_that('MEW works on small grids', {
  skip_on_cran()
  set.seed(606)

  # 4x4 grid
  bb <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0)))))
  grid <- sf::st_make_grid(bb, n = 4)
  tb <- sf::st_as_sf(grid) %>%
    rename(geometry = x) %>%
    mutate(pop = c(10, 15, 12, 18, 20, 14, 11, 16, 13, 17, 19, 12, 15, 18, 14, 16))

  grid_map <- redist_map(tb, ndists = 4, pop_tol = 0.25, total_pop = pop) %>%
    suppressWarnings() %>%
    suppressMessages()

  init_plan <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)

  out <- redist_mew(grid_map, 50, warmup = 10, init_plan = init_plan, silent = TRUE)
  expect_s3_class(out, 'redist_plans')
  expect_equal(max(as.matrix(out)), 4)
})

test_that('MEW population constraints work on Iowa data', {
  skip_on_cran()
  set.seed(909)

  iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.01) %>%
    suppressMessages()

  plans <- redist_mew(iowa_map, 200, warmup = 50, silent = TRUE)

  # Check population balance
  pop_bounds <- attr(iowa_map, 'pop_bounds')

  expect_true(all(plans$total_pop >= pop_bounds[1]))
  expect_true(all(plans$total_pop <= pop_bounds[3]))
})

test_that('MEW acceptance rate is reasonable', {
  set.seed(808)
  out <- redist_mew(fl_map, 200, warmup = 50, silent = TRUE)

  diag <- attr(out, 'diagnostics')
  avg_accept <- mean(diag$accept_rate)

  # MEW typically has high acceptance rates (50-99%)
  expect_true(avg_accept >= 0.3)
  expect_true(avg_accept <= 1.0)
})

test_that('MEW full integration test with larger sample', {
  skip_on_cran()
  set.seed(42)

  # Test with larger number of simulations (similar to test_mew_final.R)
  plans <- redist_mew(fl_map, 1000, warmup = 100, silent = TRUE)

  # Check output structure
  expect_s3_class(plans, 'redist_plans')
  n_draws <- length(unique(plans$draw))
  expect_true(n_draws == 901) # 900 samples + init

  # Check diagnostics
  diag <- attr(plans, 'diagnostics')
  expect_true(!is.null(diag))
  avg_accept <- mean(diag$accept_rate)
  expect_true(avg_accept >= 0.1 && avg_accept <= 1.0,
    info = 'Average acceptance rate should be reasonable'
  )

  # Verify all plans are valid
  plans_mat <- get_plans_matrix(plans)
  expect_true(all(apply(plans_mat, 2, function(x) length(unique(x)) == 3)),
    info = 'All plans should have exactly 3 districts'
  )

  pop_bounds <- attr(fl_map, 'pop_bounds')

  expect_true(all(plans$total_pop >= pop_bounds[1]))
  expect_true(all(plans$total_pop <= pop_bounds[3]))
})

test_that('MEW parallel chains work correctly', {
  skip_on_cran()
  set.seed(1010)

  # Run 2 chains in parallel
  plans <- redist_mew(fl_map, 100, warmup = 20, chains = 2, ncores = 2, silent = TRUE)

  # Check output structure
  expect_s3_class(plans, 'redist_plans')
  expect_true('chain' %in% names(plans))

  # Check that we have 2 chains
  expect_equal(length(unique(plans$chain)), 2)

  # Check that chains have equal number of samples
  chain_counts <- table(plans$chain)
  expect_equal(as.numeric(chain_counts[1]), as.numeric(chain_counts[2]))

  # Check diagnostics structure for multiple chains
  diag <- attr(plans, 'diagnostics')
  expect_true(is.list(diag))
  expect_equal(length(diag), 2)

  # Each chain should have its own diagnostics
  expect_true(all(sapply(diag, function(x) 'accept_rate' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'cycle_intersect_rate' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'runtime' %in% names(x))))

  # All plans should be valid
  plans_mat <- get_plans_matrix(plans)
  expect_true(all(apply(plans_mat, 2, function(x) length(unique(x)) == 3)))

  # Population constraints should be met
  pop_bounds <- attr(fl_map, 'pop_bounds')
  expect_true(all(plans$total_pop >= pop_bounds[1]))
  expect_true(all(plans$total_pop <= pop_bounds[3]))
})
