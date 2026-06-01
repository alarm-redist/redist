# Tests for confidence-interval helpers in R/confint.R and
# redist_quantile_trunc in R/redist_smc.R.

test_that("redist_quantile_trunc caps values above a self-determined quantile", {
    set.seed(1)
    x <- c(1:20, 1e6)
    trunc <- redist_quantile_trunc(x)
    expect_length(trunc, length(x))
    # the huge outlier is pulled down to the truncation point
    expect_lt(max(trunc), 1e6)
    # smaller values are unchanged
    expect_equal(head(trunc, 10), head(x, 10))
})


test_that("redist_smc_ci returns a 1-row, 3-column tibble bracketing the mean", {
    set.seed(11)
    plans <- redist_smc(fl_map, 200, runs = 2, silent = TRUE) |> suppressWarnings()
    plans <- plans |> mutate(value = total_pop / sum(fl25$pop))

    ci <- redist_smc_ci(plans, value, district = 1)
    expect_s3_class(ci, "tbl_df")
    expect_equal(nrow(ci), 1)
    expect_equal(ncol(ci), 3)
    expect_true(ci$value_lower <= ci$value)
    expect_true(ci$value <= ci$value_upper)
})


test_that("redist_ci dispatches to the algorithm-specific CI", {
    set.seed(11)
    plans <- redist_smc(fl_map, 200, runs = 2, silent = TRUE) |> suppressWarnings()
    plans <- plans |> mutate(value = total_pop / sum(fl25$pop))

    ci_dispatch <- suppressWarnings(redist_ci(plans, value, district = 1))
    ci_direct <- suppressWarnings(redist_smc_ci(plans, value, district = 1))
    expect_equal(ci_dispatch, ci_direct)
})


test_that("redist_ci errors when the algorithm attribute is missing", {
    set.seed(11)
    plans <- redist_smc(fl_map, 50, silent = TRUE)
    attr(plans, "algorithm") <- NULL
    expect_error(redist_ci(plans, total_pop), "algorithm")
})


test_that("redist_mcmc_ci returns a 3-column tibble bracketing the mean", {
    skip_if_not_installed("coda")
    set.seed(2)
    plans <- redist_mergesplit(fl_map, 100, warmup = 20, silent = TRUE) |>
        mutate(value = total_pop / sum(fl25$pop))

    ci <- suppressWarnings(redist_mcmc_ci(plans, value, district = 1))
    expect_s3_class(ci, "tbl_df")
    expect_equal(ncol(ci), 3)
    expect_true(ci$value_lower <= ci$value)
    expect_true(ci$value <= ci$value_upper)
})
