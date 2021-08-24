test_that("Population parity is computed correctly inside max_dev", {
    dev = c(0.066166599064, 0.036345355141, 0.036345355141,
            0.022645864159, 0.008523619910, 0.052398553498,
            0.095142336454, 0.094822415063, 0.095913575521, 0.075832795370)

    fl = redist_map(fl25, pop_tol=0.1, ndists=3) %>% suppressMessages()
    plans = redist_plans(plans_10[, 1:10], fl, "enumpart")

    plans = mutate(plans, dev=plan_parity(fl))
    expect_equal(rep(dev, each=3), plans$dev, tolerance=2e-5)
})

test_that("Group percentages are computed correctly", {
    pct = c(0.0664929249178703, 0.252476845951228, 0.192413021989597,
            0.155691917760251, 0.252476845951228, 0.10955590730432)

    fl = redist_map(fl25, pop_tol=0.1, ndists=3) %>% suppressMessages()
    plans = redist_plans(plans_10[, 1:2], fl, "enumpart")

    plans = mutate(plans, dev=group_frac(fl, BlackPop, pop))
    expect_equal(pct, plans$dev, tolerance=2e-5)
})
