test_that('single chain maintains backward compatibility', {
  skip_on_cran()
  set.seed(123)

  result <- redist_cyclewalk(fl_map, 300, init_name = FALSE)

  expect_s3_class(result, 'redist_plans')
  expect_false('chain' %in% names(result))
  expect_true('mcmc_accept' %in% names(result))
  expect_equal(ncol(get_plans_matrix(result)), 30)
})

test_that('multiple chains work with default ncores', {
  skip_on_cran()
  set.seed(456)

  result <- redist_cyclewalk(fl_map, 300, chains = 2, init_name = FALSE)

  expect_s3_class(result, 'redist_plans')
  expect_true('chain' %in% names(result))
  expect_equal(length(unique(result$chain)), 2)
  expect_equal(sum(result$chain == 1), sum(result$chain == 2))
  expect_equal(ncol(get_plans_matrix(result)), 60)
})

test_that('chains with explicit ncores', {
  set.seed(789)

  result <- redist_cyclewalk(fl_map, 300, chains = 2, ncores = 2, init_name = FALSE)

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

  result <- redist_cyclewalk(fl_map, 20, chains = 2, init_plan = init_vec, init_name = 'shared_init')

  expect_s3_class(result, 'redist_plans')
  # Should have 2 chains (plus NA for shared init)
  expect_equal(length(unique(result$chain[!is.na(result$chain)])), 2)
  expect_true('shared_init' %in% result$draw)

  ndists <- attr(fl_map, 'ndists')
  expect_equal(sum(result$draw == 'shared_init'), ndists)
  # The shared init has NA for chain
  expect_true(all(is.na(result$chain[result$draw == 'shared_init'])))
})

test_that('init_plan as matrix provides one per chain', {
  skip_on_cran()
  set.seed(222)

  init <- redist_smc(fl_map, 2, silent = TRUE)
  init_mat <- get_plans_matrix(init)

  result <- redist_cyclewalk(fl_map, 20, chains = 2, init_plan = init_mat, init_name = 'chain_init')

  expect_s3_class(result, 'redist_plans')
  expect_equal(length(unique(result$chain)), 2)
  expect_true('chain_init 1' %in% result$draw)
  expect_true('chain_init 2' %in% result$draw)
  # Each init should appear only once (ndists rows per plan)
  ndists <- attr(fl_map, 'ndists')
  expect_equal(sum(result$draw == 'chain_init 1'), ndists)
  expect_equal(sum(result$draw == 'chain_init 2'), ndists)
  # Each should be in its respective chain
  expect_equal(sum(result$draw == 'chain_init 1' & result$chain == 1), ndists)
  expect_equal(sum(result$draw == 'chain_init 2' & result$chain == 2), ndists)
})

test_that("init_plan='sample' generates unique inits per chain", {
  skip_on_cran()
  set.seed(333)

  result <- redist_cyclewalk(fl_map, 20, chains = 2, init_plan = 'sample', init_name = 'sampled')

  expect_s3_class(result, 'redist_plans')
  expect_equal(length(unique(result$chain[!is.na(result$chain)])), 2)
  expect_true('sampled 1' %in% result$draw)
  expect_true('sampled 2' %in% result$draw)

  plans_matrix <- get_plans_matrix(result)
  init1_col <- which(colnames(plans_matrix) == 'sampled 1')
  init2_col <- which(colnames(plans_matrix) == 'sampled 2')

  expect_length(init1_col, 1)
  expect_length(init2_col, 1)

  init1 <- plans_matrix[, init1_col]
  init2 <- plans_matrix[, init2_col]

  expect_false(identical(init1, init2))
  expect_true(any(init1 != init2))
})

test_that('return_all=FALSE returns only final plans', {
  set.seed(444)

  result <- redist_cyclewalk(fl_map, 50, chains = 2, return_all = FALSE, init_name = FALSE)

  expect_s3_class(result, 'redist_plans')
  expect_equal(ncol(get_plans_matrix(result)), 2)
  expect_equal(length(unique(result$chain)), 2)
})

test_that('warmup and thin work with chains', {
  skip_on_cran()
  set.seed(555)

  result <- redist_cyclewalk(fl_map, 60, chains = 2, warmup = 20, thin = 2, init_name = FALSE)

  expect_s3_class(result, 'redist_plans')
  expect_true('chain' %in% names(result))
  expected_per_chain <- (60 - 20) / 2
  expect_equal(ncol(get_plans_matrix(result)), expected_per_chain * 2)
})

test_that('edge weights work with chains', {
  skip_on_cran()
  set.seed(666)

  ew <- list(list(edge = c(1, 2), weight = 2.0), list(edge = c(1, 4), weight = 3.0))

  result <- redist_cyclewalk(fl_map, 300, chains = 2, edge_weights = ew, init_name = FALSE)

  expect_s3_class(result, 'redist_plans')
  expect_equal(length(unique(result$chain)), 2)
  # Test that edge weights don't break chain execution
  expect_equal(ncol(get_plans_matrix(result)), 60) # 30 per chain × 2 chains
  expect_true(all(!is.na(result$mcmc_accept)))
})

test_that('diagnostics collected per chain', {
  skip_on_cran()
  set.seed(777)

  result <- redist_cyclewalk(fl_map, 50, chains = 2, init_name = FALSE)

  diag <- attr(result, 'diagnostics')
  expect_type(diag, 'list')
  expect_length(diag, 2) # One diagnostic list per chain
  # Each chain's diagnostics should have the expected fields
  expect_true(all(sapply(diag, function(x) 'accept_prob' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'cycle_length' %in% names(x))))
  expect_true(all(sapply(diag, function(x) 'n_valid_cuts' %in% names(x))))
  # Each chain's diagnostic vectors should have length nsims
  expect_true(all(sapply(diag, function(x) length(x$accept_prob) == 50)))
})

test_that('chains validation rejects invalid values', {
  expect_error(redist_cyclewalk(fl_map, 30, chains = 0), 'chains.*must be positive')
  expect_error(redist_cyclewalk(fl_map, 30, chains = -1), 'chains.*must be positive')
})

test_that('init_plan matrix validation', {
  init <- redist_smc(fl_map, 1, silent = TRUE)
  init_mat <- matrix(get_plans_matrix(init)[, 1], ncol = 1)

  expect_error(
    redist_cyclewalk(fl_map, 30, chains = 2, init_plan = init_mat),
    'init_plan.*matrix must have 2 column'
  )
})

test_that('reference plans added correctly with chains', {
  skip_on_cran()
  set.seed(121)

  init <- redist_smc(fl_map, 1, silent = TRUE)
  init_vec <- get_plans_matrix(init)[, 1]

  result <- redist_cyclewalk(fl_map, 20, chains = 2, init_plan = init_vec, init_name = 'myinit')

  expect_true('myinit' %in% result$draw)
})

test_that('mh_acceptance calculated per chain', {
  skip_on_cran()
  set.seed(131)

  result <- redist_cyclewalk(fl_map, 50, chains = 2)

  mh <- attr(result, 'mh_acceptance')
  expect_type(mh, 'double')
  expect_length(mh, 2)
  expect_true(all(mh >= 0 & mh <= 1))
})

test_that('summary works with single chain', {
  skip_on_cran()
  set.seed(141)

  result <- redist_cyclewalk(fl_map, 50, init_name = FALSE)

  expect_output(summ <- summary(result))
  expect_type(summ, 'list')
})

test_that('summary works with multiple chains', {
  skip_on_cran()
  set.seed(151)

  result <- redist_cyclewalk(fl_map, 50, chains = 2, init_name = FALSE)

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

test_that('summary works for cyclewalk diagnostics', {
  skip_on_cran()
  set.seed(161)

  result <- redist_cyclewalk(fl_map, 100, chains = 2, init_name = FALSE)

  expect_no_error(summ <- summary(result))
  expect_type(summ, 'list')

  expect_true(!is.null(attr(result, 'diagnostics')))
  diagn <- attr(result, 'diagnostics')
  expect_true(is.list(diagn))
  expect_equal(length(diagn), 2) # 2 chains
})
