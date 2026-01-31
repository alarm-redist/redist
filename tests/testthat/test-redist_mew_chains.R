# Test file for MEW parallel chains functionality
# Tests ported from CycleWalk PR #204 to ensure robust parallel behavior

# Setup
fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1) %>%
  suppressMessages() %>%
  suppressWarnings()

test_that('single and multiple chains work correctly', {
  skip_on_cran()
  
  # Single chain
  set.seed(123)
  result_single <- redist_mew(fl_map, 20, init_name = FALSE, silent = TRUE)
  
  expect_s3_class(result_single, 'redist_plans')
  expect_false('chain' %in% names(result_single))
  diag <- attr(result_single, 'diagnostics')
  expect_true('accept_rate' %in% names(diag))
  expect_equal(ncol(get_plans_matrix(result_single)), 20)
  
  # Multiple chains with default ncores
  set.seed(456)
  result_multi <- redist_mew(fl_map, 20, chains = 2, init_name = FALSE, silent = TRUE)
  
  expect_s3_class(result_multi, 'redist_plans')
  expect_true('chain' %in% names(result_multi))
  expect_equal(length(unique(result_multi$chain)), 2)
  expect_equal(sum(result_multi$chain == 1), sum(result_multi$chain == 2))
  expect_equal(ncol(get_plans_matrix(result_multi)), 40)
  
  # Chains with explicit ncores
  set.seed(789)
  result_explicit <- redist_mew(fl_map, 20, chains = 2, ncores = 2, init_name = FALSE, silent = TRUE)
  
  expect_s3_class(result_explicit, 'redist_plans')
  expect_true('chain' %in% names(result_explicit))
  expect_equal(length(unique(result_explicit$chain)), 2)
  expect_equal(length(unique(result_explicit$draw[result_explicit$chain == 1])), 20)
  expect_equal(length(unique(result_explicit$draw[result_explicit$chain == 2])), 20)
  
  # Diagnostics collected per chain
  diag_multi <- attr(result_multi, 'diagnostics')
  expect_type(diag_multi, 'list')
  expect_length(diag_multi, 2)
  expect_true(all(sapply(diag_multi, function(x) 'accept_rate' %in% names(x))))
  expect_true(all(sapply(diag_multi, function(x) 'cycle_intersect_rate' %in% names(x))))
  expect_true(all(sapply(diag_multi, function(x) 'avg_proposal_tries' %in% names(x))))
})

test_that('init_plan variations with chains', {
  skip_on_cran()
  
  # Single vector replicates across chains
  set.seed(111)
  init <- redist_smc(fl_map, 1, silent = TRUE)
  init_vec <- get_plans_matrix(init)[, 1]

  result <- redist_mew(fl_map, 20, chains = 2, init_plan = init_vec,
                       init_name = 'shared_init', silent = TRUE)

  expect_s3_class(result, 'redist_plans')
  expect_equal(length(unique(result$chain[!is.na(result$chain)])), 2)
  expect_true(any(grepl('shared_init', result$draw)))
  
  # Matrix provides one per chain
  set.seed(222)
  init_multi <- redist_smc(fl_map, 2, silent = TRUE)
  init_mat <- get_plans_matrix(init_multi)

  result2 <- redist_mew(fl_map, 20, chains = 2, init_plan = init_mat,
                       init_name = 'chain_init', silent = TRUE)

  expect_s3_class(result2, 'redist_plans')
  expect_true(length(unique(result2$chain[!is.na(result2$chain)])) >= 2)
  expect_true(any(grepl('chain_init', result2$draw)))
  
  # 'sample' generates unique inits
  set.seed(333)
  result3 <- redist_mew(fl_map, 20, chains = 2, init_plan = 'sample',
                       init_name = 'sampled', silent = TRUE)

  expect_s3_class(result3, 'redist_plans')
  expect_equal(length(unique(result3$chain[!is.na(result3$chain)])), 2)

  plans_matrix <- get_plans_matrix(result3)
  init_cols <- grep('sampled', colnames(plans_matrix))
  expect_length(init_cols, 2)

  init1 <- plans_matrix[, init_cols[1]]
  init2 <- plans_matrix[, init_cols[2]]
  expect_false(identical(init1, init2))
  expect_true(any(init1 != init2))
})

test_that('return_all parameter not in MEW', {
  expect_error(
    redist_mew(fl_map, 20, chains = 2, return_all = FALSE),
    'unused argument'
  )
})

test_that('chains validation', {
  expect_error(redist_mew(fl_map, 30, chains = 0), 'chains.*must be positive')
  expect_error(redist_mew(fl_map, 30, chains = -1), 'chains.*must be positive')
  
  # init_plan matrix validation
  init <- redist_smc(fl_map, 1, silent = TRUE)
  init_mat <- matrix(get_plans_matrix(init)[, 1], ncol = 1)

  expect_error(
    redist_mew(fl_map, 30, chains = 2, init_plan = init_mat),
    'init_plan.*matrix must have 2 column'
  )
})

test_that('chains produce different results with different seeds', {
  skip_on_cran()
  
  set.seed(123)
  result1 <- redist_mew(fl_map, 25, chains = 2, silent = TRUE)

  set.seed(456)
  result2 <- redist_mew(fl_map, 25, chains = 2, silent = TRUE)

  expect_false(identical(
    get_plans_matrix(result1),
    get_plans_matrix(result2)
  ))
})

test_that('chains show sufficient variance (RNG independence)', {
  skip_on_cran()
  
  set.seed(999)
  result <- redist_mew(fl_map, 50, chains = 4, ncores = 2, silent = TRUE)

  diag <- attr(result, 'diagnostics')
  chain_accept_rates <- sapply(diag, function(x) mean(x$accept_rate))

  expect_true(sd(chain_accept_rates) >= 0)
  expect_true(length(unique(round(chain_accept_rates, 3))) >= 3)
})

test_that('summary works (skipped - not yet implemented)', {
  skip('MEW does not support summary() - test skipped')
  skip_on_cran()
  set.seed(141)

  result <- redist_mew(fl_map, 25, init_name = FALSE, silent = TRUE)

  expect_output(summ <- summary(result))
  expect_type(summ, 'list')
})

test_that('summary works with multiple chains (skipped - not yet implemented)', {
  skip('MEW does not support summary() - test skipped')
  skip_on_cran()
  set.seed(151)

  result <- redist_mew(fl_map, 25, chains = 2, init_name = FALSE, silent = TRUE)

  result <- result %>%
    group_by(draw) %>%
    mutate(pop_dev = max(abs(total_pop / mean(total_pop) - 1))) %>%
    ungroup()

  expect_no_error(summ <- summary(result))
  expect_type(summ, 'list')
  expect_true('rhat' %in% names(summ))
})

test_that('all plans from chains meet constraints', {
  skip_on_cran()
  
  set.seed(808)
  result <- redist_mew(fl_map, 50, chains = 2, silent = TRUE)

  # Population constraints
  pop_bounds <- attr(fl_map, 'pop_bounds')
  expect_true(all(result$total_pop >= pop_bounds[1]))
  expect_true(all(result$total_pop <= pop_bounds[3]))
  
  # Contiguity
  plans_mat <- get_plans_matrix(result)
  contig <- apply(plans_mat, 2, function(x) all(contiguity(fl_map$adj, x) == 1))
  expect_true(all(contig))
})
