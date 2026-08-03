test_that('hierarchy levels are normalized into nested integer labels', {
  adj <- list(c(1L), c(0L), c(3L), c(2L))

  levels <- normalize_hierarchy_levels(
    adj,
    list(c('north', 'north', 'north', 'north'), c('a', 'a', 'b', 'b'))
  )

  expect_equal(levels, matrix(c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L), ncol = 2L))
})

test_that('hierarchy enforcement rejects an initial plan with a hierarchy cycle', {
  expect_error(
    redist_mmss(
      fl_map,
      3L,
      init_plan = plans_10[, 1L],
      counties = rep('all', nrow(fl_map)),
      hierarchy_mode = 'strict',
      silent = TRUE
    ),
    'strict hierarchical-plan condition'
  )
})
