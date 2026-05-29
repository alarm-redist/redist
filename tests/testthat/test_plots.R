test_that("redist.plot.map works", {
    out <- redist.plot.map(shp = iowa, plan = iowa$cd_2010)
    expect_true("ggplot" %in% class(out))

    iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)

    out <- iowa_map %>% redist.plot.map(shp = ., plan = get_existing(.))
    expect_true("ggplot" %in% class(out))

    out <- iowa_map %>% redist.plot.map(shp = ., plan = cd_2010)
    expect_true("ggplot" %in% class(out))

    out <- iowa_map %>% redist.plot.map(shp = ., plan = cd_2010, fill = white)
    expect_true("ggplot" %in% class(out))
})


test_that("redist.plot.adj works", {
    out <- redist.plot.adj(shp = iowa, plan = iowa$cd_2010)
    expect_true("ggplot" %in% class(out))

    iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)

    out <- iowa_map %>% redist.plot.adj(shp = ., plan = get_existing(.))
    expect_true("ggplot" %in% class(out))

    out <- iowa_map %>% redist.plot.map(shp = ., plan = cd_2010)
    expect_true("ggplot" %in% class(out))
})

# Smoke tests for the remaining redist.plot.* functions. These mainly guard
# against silent breakage in plotting helpers — they only assert that a ggplot
# (or patchwork) object is returned, not anything about the visual output.

test_that("redist.plot.majmin returns a ggplot for each supported type", {
    set.seed(1)
    gp <- matrix(runif(4 * 10), nrow = 4)
    expect_s3_class(redist.plot.majmin(gp, type = "hist"), "ggplot")
    expect_s3_class(redist.plot.majmin(gp, type = "box"), "ggplot")
    expect_error(redist.plot.majmin(gp, type = "bogus"))
})

test_that("redist.plot.cores returns a ggplot", {
    iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) %>%
        suppressMessages()
    cores <- suppressWarnings(redist.identify.cores(
        adj = get_adj(iowa_map),
        plan = iowa_map$cd_2010
    ))
    out <- redist.plot.cores(shp = iowa_map, plan = cd_2010, core = cores)
    expect_s3_class(out, "ggplot")
})

test_that("redist.plot.varinfo returns a patchwork ggplot", {
    skip_on_cran()
    iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) %>%
        suppressMessages()
    set.seed(2)
    # Need >= 5 distinct plans for kmeans(centers = 5); jitter the existing one
    plans <- replicate(8, {
        p <- iowa_map$cd_2010
        i <- sample.int(length(p), 2)
        p[i] <- p[rev(i)]
        p
    })
    out <- suppressWarnings(
        redist.plot.varinfo(plans, iowa_map$bvap, iowa_map$vap, iowa_map)
    )
    expect_true(inherits(out, "ggplot") || inherits(out, "patchwork"))
})

test_that("plot.redist_constr handles supported constraint types", {
    iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) %>%
        suppressMessages()
    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_hinge(10, bvap, vap, tgts_group = 0.5) %>%
        add_constr_grp_inv_hinge(5, bvap, vap, tgts_group = 0.3) %>%
        add_constr_grp_pow(1, bvap, vap, tgt_group = 0.55, tgt_other = 0.25)
    expect_s3_class(plot(constr), "ggplot")
    expect_error(plot(constr, type = "bogus"), "group")
})
