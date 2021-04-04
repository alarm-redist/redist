test_that('group percent works', {
  plans <- plans_10[, 1:3]

  out <- redist.group.percent(plans = plans, group_pop = fl25$vap, total_pop = fl25$pop)

  expected <- structure(c(
    0.801126874300292, 0.780898971875266, 0.672846008005056,
    0.724416801454036, 0.780898971875266, 0.742152346737333, 0.712518899818222,
    0.791536703751272, 0.742152346737333
  ), .Dim = c(3L, 3L))

  expect_equal(out, expected)

  # Make sure zero indexing doesn't break it:
  out <- redist.group.percent(plans = plans - 1, group_pop = fl25$vap, total_pop = fl25$pop)

  expected <- structure(c(
    0.801126874300292, 0.780898971875266, 0.672846008005056,
    0.724416801454036, 0.780898971875266, 0.742152346737333, 0.712518899818222,
    0.791536703751272, 0.742152346737333
  ), .Dim = c(3L, 3L))

  expect_equal(out, expected)
})
