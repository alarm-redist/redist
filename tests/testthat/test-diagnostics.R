# Tests for chain-diagnostic helpers and multi-chain MH acceptance.

test_that("diag_fold and diag_ranknorm behave as documented", {
    set.seed(1)
    x <- c(-2, -1, 0, 1, 2)
    expect_equal(redist:::diag_fold(x), abs(x - median(x)))

    rn <- redist:::diag_ranknorm(x)
    expect_length(rn, length(x))
    expect_true(all(is.finite(rn)))
    # rank-normalized values are increasing in their input rank
    expect_equal(order(rn), order(x))
})

test_that("diag_calc_rhat ~ 1 for iid chains and > 1 when chains drift apart", {
    set.seed(42)
    n_per <- 200; n_chain <- 4
    grp <- rep(seq_len(n_chain), each = n_per)

    # iid: rhat ~ 1
    x_iid <- rnorm(n_per * n_chain)
    expect_lt(redist:::diag_calc_rhat(x_iid, grp), 1.1)

    # drifted: shift each chain mean far apart -> large rhat
    x_drift <- x_iid + rep(c(0, 5, 10, 15), each = n_per)
    expect_gt(redist:::diag_calc_rhat(x_drift, grp), 2)
})

test_that("diag_rhat returns the max over folded and ranknormed Rhat", {
    set.seed(7)
    x <- rnorm(800)
    grp <- rep(1:4, each = 200)
    r <- redist:::diag_rhat(x, grp)
    expect_true(is.finite(r))
    expect_gt(r, 0)

    # `split = TRUE` splits each chain in half
    r_split <- redist:::diag_rhat(x, grp, split = TRUE)
    expect_true(is.finite(r_split))
})

test_that("redist_mergesplit reports a valid MH acceptance rate", {
    set.seed(5118)
    pl <- redist_mergesplit(fl_map, nsims = 30, warmup = 10, silent = TRUE)
    mh <- get_mh_acceptance_rate(pl)
    expect_true(is.numeric(mh))
    expect_true(is.finite(mh))
    expect_true(mh >= 0 && mh <= 1)
})

test_that("redist_mergesplit_parallel runs and yields one diagnostics entry per chain", {
    skip_on_os("windows")
    set.seed(5118)
    pl <- redist_mergesplit_parallel(
        fl_map, nsims = 30, warmup = 10, chains = 2,
        ncores = 2, silent = TRUE
    )
    diag <- attr(pl, "diagnostics")
    expect_length(diag, 2)
})

test_that("redist_smc with runs > 1 yields one diagnostics entry per run", {
    set.seed(5118)
    pl <- redist_smc(fl_map, 50, runs = 2, silent = TRUE)
    diag <- attr(pl, "diagnostics")
    expect_length(diag, 2)
})
