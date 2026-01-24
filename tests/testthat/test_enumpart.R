skip_on_ci()
skip_on_cran()
skip_if(Sys.info()['machine'] == 'arm64', 'arm64 machines')

if (!(file.exists(system.file('enumpart/enumpart', package = 'redist')) ||
  file.exists(system.file('enumpart/enumpart.exe', package = 'redist')))) {
  redist.init.enumpart()
}
dir <- withr::local_tempdir()

test_that('enumpart preparation runs correctly', {
  capture.output(
    redist.prep.enumpart(
      adj = adj,
      ordered_path = file.path(dir, 'ordered')
    )
  )
  ord <- scan(file.path(dir, 'ordered.dat'))

  # 102 = 51 edges * 2 vertices
  expect_equal(length(ord), 102)
  expect_equal(min(ord), 1)
  expect_equal(max(ord), 25)
})

test_that('enumpart can sample without constraints', {
  sample_path <- file.path(dir, 'sample.dat')
  if (file.exists(sample_path)) {
    file.remove(sample_path)
  }
  capture.output(
    redist.run.enumpart(
      ordered_path = file.path(dir, 'ordered'),
      out_path = file.path(dir, 'sample'),
      ndists = 3, all = F, n = 10
    )
  )
  m <- matrix(scan(file.path(dir, 'sample.dat')), nrow = 25)

  expect_equal(dim(m), c(25, 10))
  expect_equal(range(m), c(0, 2))
})

test_that('enumpart can read in data', {
  capture.output(
    full <- redist.read.enumpart(out_path = file.path(dir, 'sample'))
  )

  expect_equal(nrow(full), 25)
  expect_equal(ncol(full), 10)
  expect_equal(range(full), c(1, 3))
})

test_that('enumpart can read in partial data', {
  full <- redist.read.enumpart(out_path = file.path(dir, 'sample'))

  samp <- redist.read.enumpart(
    out_path = file.path(dir, 'sample'),
    skip = 5
  )

  expect_equal(nrow(samp), 25)
  expect_equal(ncol(samp), 5)
  expect_equal(range(samp), c(1, 3))
  expect_true(all(full[, 6:10] == samp))

  samp <- redist.read.enumpart(
    out_path = file.path(dir, 'sample'),
    n_max = 8
  )

  expect_equal(nrow(samp), 25)
  expect_equal(ncol(samp), 8)
  expect_equal(range(samp), c(1, 3))
  expect_true(all(full[, 1:8] == samp))
})


test_that('enumpart can sample with unit count constraints', {
  sample_path <- file.path(dir, 'sample.dat')
  if (file.exists(sample_path)) {
    file.remove(sample_path)
  }

  capture.output(
    redist.run.enumpart(file.path(dir, 'ordered'), file.path(dir, 'sample'),
      ndists = 3,
      all = F, n = 100, lower = 4, upper = 16
    )
  )
  m <- matrix(scan(sample_path), nrow = 25)

  expect_equal(dim(m), c(25, 100))
  expect_equal(range(m), c(0, 2))
  # check component sizes
  range_sizes <- range(apply(m, 2, function(x) range(table(x))))
  expect_equal(range_sizes, c(4, 16))
})

test_that('enumpart can sample with population constraints', {
  sample_path <- file.path(dir, 'sample.dat')
  if (file.exists(sample_path)) file.remove(sample_path)

  write(pop, file.path(dir, 'pop.dat'), ncolumns = length(pop))
  target <- sum(pop) / 3
  capture.output(
    redist.run.enumpart(file.path(dir, 'ordered'), file.path(dir, 'sample'),
      ndists = 3, all = T, lower = round(target * 0.9), upper = round(target * 1.1),
      weight_path = file.path(dir, 'pop')
    )
  )
  m <- matrix(scan(sample_path), nrow = 25)

  expect_equal(dim(m), c(25, 927))
  expect_equal(range(m), c(0, 2))
  # check component sizes
  dev <- redist.parity(m, pop)
  expect_true(max(dev) <= 0.1)
})



withr::deferred_clear()
