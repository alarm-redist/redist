test_that("constructor works", {
    fl = redist_map(fl25, ndists=3, pop_tol=0.1) %>% suppressMessages()
    x = redist_plans(plans_10, fl, "enumpart")
    expect_s3_class(x, "redist_plans")
    expect_s3_class(x, "tbl")
    expect_true(is.matrix(as.matrix(x)))
    expect_equal(attr(x, "prec_pop"), fl25$pop)
    expect_null(weights(x))
    expect_output(print(x), "drawn using Enumpart", fixed=T)
    expect_output(print(x), "927 sampled plans with 3 districts", fixed=T)
})

test_that("add_reference works", {
    fl = redist_map(fl25, ndists=3, pop_tol=0.1) %>% suppressMessages()
    x = redist_plans(plans_10, fl, "enumpart") %>%
        add_reference(plans_10[,1])
    expect_s3_class(x, "redist_plans")
    expect_equal(as.character(head(x$draw, 3)), rep("plans_10[, 1]", 3))
    expect_equal(colnames(as.matrix(x))[1], "plans_10[, 1]")
    expect_equal(ncol(as.matrix(x)), ncol(plans_10)+1)
    expect_equal(nrow(as.matrix(x)), nrow(plans_10))
})

test_that("subsetting works", {
    fl = redist_map(fl25, ndists=3, pop_tol=0.1) %>% suppressMessages()
    x = redist_plans(plans_10, fl, "enumpart") %>%
        add_reference(plans_10[,1])
    x_samp = subset_sampled(x)
    x_ref = subset_ref(x)
    expect_equal(ncol(as.matrix(x_samp)), ncol(plans_10))
    expect_equal(ncol(as.matrix(x_ref)), 1)
})

test_that("plotting works", {
    fl = redist_map(fl25, ndists=3, pop_tol=0.1) %>% suppressMessages()
    x = redist_plans(plans_10, fl, "enumpart") %>%
        add_reference(plans_10[,1])
    out = plot(x)
    expect_s3_class(out, "ggplot")
    out = hist(x, total_pop)
    expect_s3_class(out, "ggplot")
})

test_that("get_mh_acceptance_rate works", {
    iowa_map <- redist_map(iowa, pop_tol = 0.01, existing_plan = cd_2010)

    # should exist and be number
    capture.output(
    out <- redist_flip(iowa_map, nsims = 5, verbose = FALSE)
    )
    mh <- get_mh_acceptance_rate(out)
    expect_true(is.numeric(mh))
    expect_true(mh >= 0)

    out <- redist_mergesplit(iowa_map, nsims = 5, silent = TRUE)
    mh <- get_mh_acceptance_rate(out)
    expect_true(is.numeric(mh))
    expect_true(mh >= 0)

    # should return NA_real_
    out <- redist_smc(iowa_map, nsims = 5, silent = TRUE)
    mh <- get_mh_acceptance_rate(out)
    expect_true(is.na(mh))
    expect_true(is.numeric(mh))
})
