test_that("group percent works", {
    plans <- redist_plans(plans_10[, 1:3], fl_map, "enumpart")

    # out <- redist.group.percent(plans = plans, group_pop = fl25$vap, total_pop = fl25$pop)
    out <- group_frac(fl_map, vap, pop, .data=plans)

    expected <- c(
        0.801126874300292, 0.780898971875266, 0.672846008005056,
        0.724416801454036, 0.780898971875266, 0.742152346737333, 0.712518899818222,
        0.791536703751272, 0.742152346737333
    )

    expect_equal(out, expected)

    # check providing raw shp works too
    out <- group_frac(fl25, vap, pop, .data=plans)
    expect_equal(out, expected)
})
