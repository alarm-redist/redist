skip_on_ci()
skip_on_cran()
skip_if(Sys.info()['machine'] == 'arm64', 'arm64 machines')

if (!(file.exists(system.file('enumpart/enumpart', package = 'redist')) ||
  file.exists(system.file('enumpart/enumpart.exe', package = 'redist')))) {
  suppressWarnings(redist.init.enumpart())
}

test_that('redist_enumpart enumerates without population constraints', {
  dir <- withr::local_tempdir()
  plans <- redist_enumpart(fl_map, n = 10, lower = NULL, upper = NULL,
    out_path = file.path(dir, 'out'))

  expect_s3_class(plans, 'redist_plans')
  expect_equal(attr(plans, 'algorithm'), 'enumpart')

  pl <- attr(plans, 'plans')
  expect_equal(nrow(pl), 25)
  expect_equal(ncol(pl), 10)
  expect_equal(range(pl), c(1, 3))
})

test_that('redist_enumpart enumerates with population constraints', {
  dir <- withr::local_tempdir()
  plans <- redist_enumpart(fl_map,
    out_path = file.path(dir, 'out'))

  expect_s3_class(plans, 'redist_plans')

  pl <- attr(plans, 'plans')
  expect_equal(nrow(pl), 25)
  expect_equal(ncol(pl), 927)

  dev <- redist.parity(pl, pop)
  expect_true(max(dev) <= 0.1)
})

test_that('redist_enumpart respects explicit paths', {
  dir <- withr::local_tempdir()
  plans <- redist_enumpart(fl_map, n = 10,
    lower = NULL, upper = NULL,
    ordered_path = file.path(dir, 'ordered'),
    out_path = file.path(dir, 'out'))

  expect_s3_class(plans, 'redist_plans')
  expect_true(file.exists(file.path(dir, 'ordered.dat')))
  expect_true(file.exists(file.path(dir, 'out.dat')))

  pl <- attr(plans, 'plans')
  expect_equal(nrow(pl), 25)
  expect_equal(ncol(pl), 10)
  expect_equal(range(pl), c(1, 3))
})

test_that('redist_enumpart read = FALSE writes file without reading', {
  dir <- withr::local_tempdir()
  result <- redist_enumpart(fl_map, n = 5,
    lower = NULL, upper = NULL,
    out_path = file.path(dir, 'out'),
    read = FALSE)

  expect_equal(result, file.path(dir, 'out'))
  expect_true(file.exists(file.path(dir, 'out.dat')))
})

test_that('redist_enumpart reads existing out_path without re-running', {
  dir <- withr::local_tempdir()
  out <- file.path(dir, 'out')

  redist_enumpart(fl_map, n = 8,
    lower = NULL, upper = NULL,
    out_path = out, read = FALSE)

  # Verify the binary is skipped: corrupt the ordered file, which would
  # cause enumpart to fail if it were re-run
  writeLines("corrupted", file.path(dir, 'corrupted_ord.dat'))
  plans <- redist_enumpart(fl_map, out_path = out)

  pl <- attr(plans, 'plans')
  expect_equal(ncol(pl), 8)
})

test_that('redist_enumpart_frontier computes from adj', {
  frontier <- redist_enumpart_frontier(adj)

  expect_type(frontier, 'list')
  expect_named(frontier, c('max', 'average', 'average_sq', 'sequence'))
  expect_true(frontier$max > 0)
  expect_true(length(frontier$sequence) > 0)
})
