if (!file.exists(system.file("enumpart/enumpart", package="redist"))) {
    redist.init.enumpart()
}
dir = withr::local_tempdir()

test_that("enumpart preparation runs correctly", {
    expected = c(21, 23, 3, 23, 3, 21, 6, 23, 3, 6, 4, 21, 3, 4, 1, 3, 1, 6, 1,
                 4, 2, 21, 2, 4, 1, 2, 18, 21, 2, 18, 6, 9, 1, 9, 1, 14, 9, 14,
                 2, 16, 16, 18, 1, 15, 14, 15, 1, 13, 13, 15, 2, 17, 16, 17, 1,
                 12, 12, 13, 1, 11, 11, 12, 2, 19, 17, 19, 1, 10, 10, 11, 1, 8,
                 8, 10, 2, 20, 19, 20, 1, 7, 7, 8, 1, 5, 5, 7, 2, 22, 20, 22, 5,
                 22, 5, 25, 22, 25, 5, 24, 22, 24, 24, 25)

    redist.prep.enumpart(g, file.path(dir, "unordered"), file.path(dir, "ordered"))
    expect_equal(scan(file.path(dir, "ordered.dat")), expected)
})

test_that("enumpart can sample without constraints", {
    sample_path = file.path(dir, "sample.dat")
    if (file.exists(sample_path)) file.remove(sample_path)

    redist.run.enumpart(file.path(dir, "ordered"), file.path(dir, "sample"),
                        ndists=3, all=F, n=10)
    m = matrix(scan(file.path(dir, "sample.dat")), nrow=25)

    expect_equal(dim(m), c(25, 10))
    expect_equal(range(m), c(0, 2))
})

test_that("enumpart can sample with unit count constraints", {
    sample_path = file.path(dir, "sample.dat")
    if (file.exists(sample_path)) file.remove(sample_path)

    redist.run.enumpart(file.path(dir, "ordered"), file.path(dir, "sample"),
                        ndists=3, all=F, n=100, lower=4, upper=16)
    m = matrix(scan(sample_path), nrow=25)

    expect_equal(dim(m), c(25, 100))
    expect_equal(range(m), c(0, 2))
    # check component sizes
    range_sizes = range(apply(m, 2, function(x) range(table(x))))
    expect_equal(range_sizes, c(4, 16))
})

test_that("enumpart can sample with population constraints", {
    sample_path = file.path(dir, "sample.dat")
    if (file.exists(sample_path)) file.remove(sample_path)

    write(pop, file.path(dir, "pop.dat"), ncolumns=length(pop))
    target = sum(pop) / 3
    redist.run.enumpart(file.path(dir, "ordered"), file.path(dir, "sample"),
                        ndists=3, all=T, lower=round(target*0.9), upper=round(target*1.1),
                        weight_path=file.path(dir, "pop"))
    m = matrix(scan(sample_path), nrow=25)

    expect_equal(dim(m), c(25, 927))
    expect_equal(range(m), c(0, 2))
    # check component sizes
    dev = redist.parity(m, pop)
    expect_true(max(dev) <= 0.1)
})

withr::deferred_clear()
