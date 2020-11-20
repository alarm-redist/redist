context("Metrics and Measurements")

test_that("Pairwise distances are computed correctly", {
    dist_m = list(
        Hamming=c(0, 2, 4, 2, 0, 2, 4, 2, 0),
        Manhattan=c(0, 4, 6, 4, 0, 2, 6, 2, 0),
        Euclidean=c(0, 2.82842712474619, 3.16227766016838, 2.82842712474619, 0,
                           1.4142135623731, 3.16227766016838, 1.4142135623731, 0),
        VI=c(0, 0.843773942032667, 1.34244720438262, 0.843773942032667, 0,
                    0.550888899101157, 1.34244720438262, 0.550888899101157, 0)
        )
    res = redist.distances(algdat.p10$cdmat[,1:3], measure="all", pop=pop)

    expect_equal(as.numeric(res$Hamming), dist_m$Hamming)
    expect_equal(as.numeric(res$Manhattan), dist_m$Manhattan)
    expect_equal(as.numeric(res$Euclidean), dist_m$Euclidean)
    expect_equal(as.numeric(res$VI), dist_m$VI)
})

test_that("Population parity is computed correctly", {
    dev = c(0.066155929182306, 0.036334344524997, 0.036334344524997,
            0.022657548802852, 0.00853514319502291, 0.0524105780931325,
            0.0951548494352752, 0.0948120726001336, 0.0959260973143434, 0.0758222359332956)

    res = redist.parity(algdat.p10$cdmat[,1:10], algdat.p10$precinct.data$pop)
    expect_equal(res, dev)
})

test_that("Group percentages are computed correctly", {
    pct = c(0.0664929249178703, 0.252476845951228, 0.192413021989597,
            0.155691917760251, 0.252476845951228, 0.10955590730432)
    res = redist.group.percent(algdat.p10$cdmat[,1:2],
                               algdat.p10$precinct.data$blackpop,
                               algdat.p10$precinct.data$pop)
    expect_equal(as.numeric(res), pct)
})
