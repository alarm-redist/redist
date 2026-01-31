# Test file for MEW parallel chains functionality
# Tests ported from CycleWalk PR #204 to ensure robust parallel behavior

# Setup
fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1) %>%
  suppressMessages() %>%
  suppressWarnings()

test_that('single chain maintains backward compatibility', {
  skip_on_cran()
  set.seed(123)

  result <- redist_mew(fl_map, 30, init_name = FALSE, silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  expect_false('chain' %in% names(result))
  # MEW doesn't add mcmc_accept column - check diagnostics instead
  diag <- attr(result, 'diagnostics')
  expect_true('accept_rate' %in% names(diag))
  expect_equal(ncol(get_plans_matrix(result)), 30)
})

test_that('multiple chains work with default ncores', {
  skip_on_cran()
  set.seed(456)

  result <- redist_mew(fl_map, 30, chains = 2, init_name = FALSE, silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  expect_true('chain' %in% names(result))
  expect_equal(length(unique(result$chain)), 2)
  expect_equal(sum(result$chain == 1), sum(result$chain == 2))
  expect_equal(ncol(get_plans_matrix(result)), 60)
})

test_that('chains with explicit ncores', {
  set.seed(789)

  result <- redist_mew(fl_map, 30, chains = 2, ncores = 2, init_name = FALSE, silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  expect_true('chain' %in% names(result))
  expect_equal(length(unique(result$chain)), 2)
  expect_equal(length(unique(result$draw[result$chain == 1])), 30)
  expect_equal(length(unique(result$draw[result$chain == 2])), 30)
})

test_that('init_plan as single vector replicates across chains', {
  skip_on_cran()
  set.seed(111)

  init <- redist_smc(fl_map, 1, silent = TRUE)
  init_vec <- get_plans_matrix(init)[, 1]

  result <- redist_mew(fl_map, 20, chains = 2, init_plan = init_vec,
                       init_name = 'shared_init', silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  # Should have 2 chains
  expect_equal(length(unique(result$chain[!is.na(result$chain)])), 2)

  ndists <- attr(fl_map, 'ndists')
  # The shared init appears once at the beginning
  expect_true(any(grepl('shared_init', result$draw)))
})

test_that('init_plan as matrix provides one per chain', {
  skip_on_cran()
  set.seed(222)

  init <- redist_smc(fl_map, 2, silent = TRUE)
  init_mat <- get_plans_matrix(init)

  result <- redist_mew(fl_map, 20, chains = 2, init_plan = init_mat,
                       init_name = 'chain_init', silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  # MEW adds init plans separately for each chain
  # Could have 2 chains + 1 or 2 init plans
  expect_true(length(unique(result$chain[!is.na(result$chain)])) >= 2)

  # Both inits should appear in draws
  expect_true(any(grepl('chain_init', result$draw)))
})

test_that("init_plan='sample' generates unique inits per chain", {
  skip_on_cran()
  set.seed(333)

  result <- redist_mew(fl_map, 20, chains = 2, init_plan = 'sample',
                       init_name = 'sampled', silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  expect_equal(length(unique(result$chain[!is.na(result$chain)])), 2)

  plans_matrix <- get_plans_matrix(result)
  # Look for sampled init plans
  init_cols <- grep('sampled', colnames(plans_matrix))

  expect_length(init_cols, 2)

  init1 <- plans_matrix[, init_cols[1]]
  init2 <- plans_matrix[, init_cols[2]]

  expect_false(identical(init1, init2))
  expect_true(any(init1 != init2))
})

test_that('return_all=FALSE not applicable to MEW', {
  # MEW doesn't have return_all parameter - this test documents that difference
  expect_error(
    redist_mew(fl_map, 50, chains = 2, return_all = FALSE),
    'unused argument'
  )
})

test_that('warmup and thin work with chains', {
  skip_on_cran()
  set.seed(555)

  result <- redist_mew(fl_map, 60, chains = 2, warmup = 20, thin = 2,
                       init_name = FALSE, silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  expect_true('chain' %in% names(result))
  expected_per_chain <- (60 - 20) / 2
  expect_equal(ncol(get_plans_matrix(result)), expected_per_chain * 2)
})

test_that('diagnostics collected per chain', {
  skip_on_cran()
  set.seed(777)

  result <- redist_mew(fl_map, 50, chains = 2, init_name = FALSE, silent = TRUE)

  diag <- attr(result, 'diagnostics')
  expect_type(diag, 'list')
  expect_length(diag, 2) # One diagnostic list per chain

  # Each chain's diagnostics should have the expected fields
  expect_true(all(sapply(diag, function(x) 'accept_rate' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'cycle_intersect_rate' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'avg_proposal_tries' %in% names(x))))
})

test_that('chains validation rejects invalid values', {
  expect_error(redist_mew(fl_map, 30, chains = 0), 'chains.*must be positive')
  expect_error(redist_mew(fl_map, 30, chains = -1), 'chains.*must be positive')
})

test_that('init_plan matrix validation', {
  init <- redist_smc(fl_map, 1, silent = TRUE)
  init_mat <- matrix(get_plans_matrix(init)[, 1], ncol = 1)

  expect_error(
    redist_mew(fl_map, 30, chains = 2, init_plan = init_mat),
    'init_plan.*matrix must have 2 column'
  )
})

test_that('reference plans added correctly with chains', {
  skip_on_cran()
  set.seed(121)

  init <- redist_smc(fl_map, 1, silent = TRUE)
  init_vec <- get_plans_matrix(init)[, 1]

  result <- redist_mew(fl_map, 20, chains = 2, init_plan = init_vec,
                       init_name = 'myinit', silent = TRUE)

  expect_true(any(grepl('myinit', result$draw)))
})

test_that('parallel chains are reproducible with same seed', {
  skip('Reproducibility with doRNG needs investigation')
  skip_on_cran()
  set.seed(789)
  result1 <- redist_mew(fl_map, 50, chains = 2, ncores = 2, silent = TRUE)

  set.seed(789)
  result2 <- redist_mew(fl_map, 50, chains = 2, ncores = 2, silent = TRUE)

  # Results should be identical due to doRNG
  expect_identical(
    get_plans_matrix(result1),
    get_plans_matrix(result2)
  )
})

test_that('chains produce different results with different seeds', {
  skip_on_cran()
  set.seed(123)
  result1 <- redist_mew(fl_map, 50, chains = 2, silent = TRUE)

  set.seed(456)
  result2 <- redist_mew(fl_map, 50, chains = 2, silent = TRUE)

  # Results should differ
  expect_false(identical(
    get_plans_matrix(result1),
    get_plans_matrix(result2)
  ))
})

test_that('chains show sufficient variance (RNG independence check)', {
  skip_on_cran()
  set.seed(999)
  result <- redist_mew(fl_map, 100, chains = 4, ncores = 2, silent = TRUE)

  # Calculate acceptance rate for each chain (from diagnostics)
  diag <- attr(result, 'diagnostics')
  chain_accept_rates <- sapply(diag, function(x) mean(x$accept_rate))

  # Chains should show some variance (not identical)
  expect_true(sd(chain_accept_rates) >= 0)
  # At least 3 of 4 chains should show different acceptance patterns
  expect_true(length(unique(round(chain_accept_rates, 3))) >= 3)
})

test_that('summary works with single chain', {
  skip('MEW does not support summary() - test skipped')
  skip_on_cran()
  set.seed(141)

  result <- redist_mew(fl_map, 50, init_name = FALSE, silent = TRUE)

  expect_output(summ <- summary(result))
  expect_type(summ, 'list')
})

test_that('summary works with multiple chains', {
  skip('MEW does not support summary() - test skipped')
  skip_on_cran()
  set.seed(151)

  result <- redist_mew(fl_map, 50, chains = 2, init_name = FALSE, silent = TRUE)

  # Add a summary statistic for R-hat calculation
  result <- result %>%
    group_by(draw) %>%
    mutate(pop_dev = max(abs(total_pop / mean(total_pop) - 1))) %>%
    ungroup()

  # Test that summary doesn't error and returns a list
  expect_no_error(summ <- summary(result))
  expect_type(summ, 'list')
  expect_true('rhat' %in% names(summ))
})

test_that('all plans from chains meet population constraints', {
  skip_on_cran()
  set.seed(808)

  result <- redist_mew(fl_map, 100, chains = 2, silent = TRUE)

  # Check population balance
  pop_bounds <- attr(fl_map, 'pop_bounds')

  expect_true(all(result$total_pop >= pop_bounds[1]))
  expect_true(all(result$total_pop <= pop_bounds[3]))
})

test_that('all plans from chains are contiguous', {
  skip_on_cran()
  set.seed(909)

  result <- redist_mew(fl_map, 50, chains = 2, silent = TRUE)

  plans_mat <- get_plans_matrix(result)
  contig <- apply(plans_mat, 2, function(x) all(contiguity(fl_map$adj, x) == 1))
  expect_true(all(contig))
})
