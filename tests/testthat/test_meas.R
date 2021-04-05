test_that("Pairwise distances are computed correctly", {
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
    dist_m = list(
        Hamming=c(0, 2, 4, 2, 0, 2, 4, 2, 0),
        Manhattan=c(0, 4, 6, 4, 0, 2, 6, 2, 0),
        Euclidean=c(0, 2.82842712474619, 3.16227766016838, 2.82842712474619, 0,
                           1.4142135623731, 3.16227766016838, 1.4142135623731, 0),
        VI=c(0, 0.843773942032667, 1.34244720438262, 0.843773942032667, 0,
                    0.550888899101157, 1.34244720438262, 0.550888899101157, 0)
        )
    res = redist.distances(plans_10[, 1:3], measure = "all", total_pop = pop)

    expect_equal(as.numeric(res$Hamming), dist_m$Hamming)
    expect_equal(as.numeric(res$Manhattan), dist_m$Manhattan)
    expect_equal(as.numeric(res$Euclidean), dist_m$Euclidean)
    expect_equal(as.numeric(res$VI), dist_m$VI)
})

test_that("Population parity is computed correctly inside max_dev", {
    dev = c(0.066166599064, 0.036345355141, 0.036345355141,
            0.022645864159, 0.008523619910, 0.052398553498,
            0.095142336454, 0.094822415063, 0.095913575521, 0.075832795370)

    res = max_dev(plans_10[, 1:10], pop, 3)
    expect_equal(res, dev, tolerance=2e-5)
})

test_that("Group percentages are computed correctly", {
    pct = c(0.0664929249178703, 0.252476845951228, 0.192413021989597,
            0.155691917760251, 0.252476845951228, 0.10955590730432)
    res = redist.group.percent(plans_10[, 1:2], fl25$BlackPop, fl25$pop)
    expect_equal(as.numeric(res), pct)
})
