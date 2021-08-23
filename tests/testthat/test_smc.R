test_that("SMC runs without errors", {
    N = 20
    res = redist_smc(fl_map, N, silent=T)

    expect_s3_class(res, "redist_plans")
})

test_that("County constraint works", {
    iowa_map = redist_map(iowa, ndists=4, pop_tol=0.05)
    plans = redist_smc(iowa_map, 50, counties=region, silent=TRUE)
    splits = redist.splits(as.matrix(plans), iowa_map$region)
    expect_true(all(splits <= 3L))
    expect_true(all(apply(get_plans_matrix(plans), 2,
                          function(x) all(contiguity(iowa_map$adj, x) == 1))))

    region2 = iowa$region
    region2[25] = NA
    expect_error(redist_smc(iowa_map, 50, counties=region2, silent=TRUE),
                 "missing values")
})

test_that("Single-precinct counties work", {
    bb <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0)))))
    grid <- sf::st_make_grid(bb, n = 4)
    tb <- sf::st_as_sf(grid) %>%
        rename(geometry = x) %>%
        mutate(pop = 1, counties = row_number())
    box_map <- redist_map(tb, ndists = 4, pop_tol = 0.25, total_pop = pop) %>%
        suppressWarnings() %>%
        suppressMessages()

    test <- redist_smc(box_map, 10, counties=counties, silent=TRUE)
    expect_s3_class(test, "redist_plans")
})

test_that("Not egregiously incorrect sampling accuracy (5-prec)", {
    skip_on_cran()
    set.seed(1935)

    g = list(c(1L, 4L), c(0L, 2L, 4L), c(1L, 3L, 4L), c(2L, 4L), c(0L, 1L, 2L, 3L))
    g_pop = c(2, 1, 1, 1, 1)
    map = redist_map(pop=g_pop, ndists=2, pop_tol=0.5, adj=g)
    out = redist_smc(map, 20e3, compactness=0, adapt_k_thresh=1, resample=F, silent=T)
    types = apply(as.matrix(out), 2, function(x) 1L + (x[1] == x[2]))
    wgts = weights(out)
    avg = weighted.mean(types, wgts)
    se = sqrt(sum((types - avg)^2 * (wgts / sum(wgts))^2))
    expect_true(abs(avg - 1.5)/se <= 5)
})

test_that("Not egregiously incorrect sampling accuracy (25-prec)", {
    skip_on_cran()
    set.seed(1935)

    ref_plans = plans_10[, redist.parity(plans_10, pop) <= 0.01]
    log_st_ref = round(log_st_map(adj, ref_plans, rep(1L, 25), 3L), 5)

    out = redist_smc(set_pop_tol(fl_map, 0.01), 6000, compactness=0,
                     adapt_k_thresh=0.99, seq_alpha=0.3, resample=F, silent=T) %>%
        suppressWarnings() # efficiency
    log_st = round(log_st_map(adj, as.matrix(out), rep(1L, 25), 3L), 5)
    types = match(log_st, log_st_ref)

    wgts = weights(out)
    avgs = sapply(seq_along(log_st_ref), function(i) weighted.mean(types==i, wgts))
    ses = sapply(seq_along(log_st_ref), function(i) {
        sqrt(sum(((types==i) - avgs[i])^2 * (wgts / sum(wgts))^2))
    })
    zscores = (avgs - (1/length(log_st_ref))) / ses
    expect_true(all(abs(zscores) <= 6))
})

test_that("Precise population bounds are enforced", {
    map2 = fl_map
    attr(map2, "pop_bounds") = c(52e3, 58e3, 60e3)
    res = redist_smc(map2, 20, silent=T)
    distr_pop = apply(as.matrix(res), 2, function(x) tapply(pop, x, sum))
    expect_true(all(apply(distr_pop, 2, max) <= 60e3))
    expect_true(all(apply(distr_pop, 2, min) >= 52e3))
})

test_that("SMC checks arguments", {
    expect_error(redist_smc(fl_map, 10, compactness=-1), "non-negative")
    expect_error(redist_smc(fl_map, 10, seq_alpha=1.5), "0, 1")
    expect_error(redist_smc(fl_map, 10, adapt_k_thresh=1.5), "0, 1")
    expect_error(redist_smc(fl_map, 0), "positive")
})

