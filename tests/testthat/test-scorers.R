# Tests for scorer constructors in R/redist_shortburst.R.
# These check the interface (class, output shape, sign/zero invariants) that
# the redist package layers on top of the redistmetrics numerical kernels.

iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) %>%
    suppressMessages()
plans_iowa <- matrix(iowa$cd_2010, ncol = 1)


test_that("scorer_pop_dev returns max parity deviation, zero only at perfect parity", {
    fn <- scorer_pop_dev(iowa_map)
    expect_s3_class(fn, "redist_scorer")

    out <- fn(plans_iowa)
    expect_type(out, "double")
    expect_length(out, ncol(plans_iowa))
    expect_true(all(out >= 0))

    # 4-precinct grid where the only sensible 2-district plan is perfectly balanced
    grid_map <- redist_map(
        data.frame(pop = rep(1L, 4)),
        ndists = 2,
        pop_tol = 0.5,
        adj = list(1L, c(0L, 2L), c(1L, 3L), 2L)
    )
    perfect <- matrix(c(1L, 1L, 2L, 2L), ncol = 1)
    expect_equal(scorer_pop_dev(grid_map)(perfect), 0)
})

test_that("scorer_splits returns a fraction in [0, 1]", {
    fn <- scorer_splits(iowa_map, counties = region)
    expect_s3_class(fn, "redist_scorer")

    out <- fn(plans_iowa)
    expect_length(out, ncol(plans_iowa))
    expect_true(all(out >= 0 & out <= 1))
})

test_that("scorer_multisplits is bounded above by scorer_splits", {
    fn_s <- scorer_splits(iowa_map, counties = region)
    fn_ms <- scorer_multisplits(iowa_map, counties = region)
    expect_true(all(fn_ms(plans_iowa) <= fn_s(plans_iowa) + 1e-9))
})

test_that("scorer_multisplits returns a fraction in [0, 1]", {
    fn_ms <- scorer_multisplits(iowa_map, counties = region)
    expect_s3_class(fn_ms, "redist_scorer")

    mm <- fn_ms(plans_iowa)
    expect_length(mm, ncol(plans_iowa))
    expect_true(all(mm >= 0 & mm <= 1))
})

test_that("scorer_polsby_popper returns a value in (0, 1] for each plan", {
    skip_on_cran()
    fn <- scorer_polsby_popper(iowa_map)
    expect_s3_class(fn, "redist_scorer")

    pp <- fn(plans_iowa)
    expect_length(pp, ncol(plans_iowa))
    expect_true(all(pp > 0 & pp <= 1))

    # `m` selects the m-th smallest district; larger `m` => >= smaller `m`
    pp2 <- scorer_polsby_popper(iowa_map, m = 2)(plans_iowa)
    expect_true(all(pp2 >= pp - 1e-9))
})

test_that("scorer_status_quo is 1 against the reference plan and < 1 against perturbations", {
    fn <- scorer_status_quo(iowa_map, existing_plan = iowa_map$cd_2010)
    expect_s3_class(fn, "redist_scorer")

    self <- fn(plans_iowa)
    expect_equal(self, 1, tolerance = 1e-9)

    # swap two precincts' districts to create a small perturbation
    perturbed <- plans_iowa
    perturbed[1, 1] <- if (perturbed[1, 1] == 1L) 2L else 1L
    expect_lt(fn(perturbed), 1)
})

test_that("combine_scorers / cbind.redist_scorer returns a matrix-valued scorer", {
    s1 <- scorer_frac_kept(iowa_map)
    s2 <- scorer_pop_dev(iowa_map)
    combined <- cbind(a = s1, b = s2)
    expect_s3_class(combined, "redist_scorer")

    out <- combined(plans_iowa)
    expect_true(is.matrix(out))
    expect_equal(dim(out), c(ncol(plans_iowa), 2))
    expect_equal(colnames(out), c("a", "b"))
    expect_equal(unname(out[, "a"]), s1(plans_iowa))
    expect_equal(unname(out[, "b"]), s2(plans_iowa))

    # combine_scorers is a thin wrapper around cbind
    combined2 <- combine_scorers(a = s1, b = s2)
    expect_equal(combined2(plans_iowa), out)
})

test_that("scorer arithmetic combines underlying scorers correctly", {
    s1 <- scorer_frac_kept(iowa_map)
    s2 <- scorer_pop_dev(iowa_map)

    sum_fn <- s1 + s2
    diff_fn <- s1 - s2
    scaled <- 2 * s1
    expect_s3_class(sum_fn, "redist_scorer")
    expect_s3_class(diff_fn, "redist_scorer")
    expect_s3_class(scaled, "redist_scorer")

    a <- s1(plans_iowa)
    b <- s2(plans_iowa)
    expect_equal(sum_fn(plans_iowa), a + b)
    expect_equal(diff_fn(plans_iowa), a - b)
    expect_equal(scaled(plans_iowa), 2 * a)
})
