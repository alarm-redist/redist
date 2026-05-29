# Unit tests for redist_constr constraint constructors.
# These check the structure built into the redist_constr object (strength,
# stored parameters, tidy-evaluation of bare column names) and the input
# validation around it. Numerical correctness of the resulting Gibbs penalties
# is exercised by the SMC/MS/flip distributional tests.

iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05) %>%
    suppressMessages()

test_that("add_constr_status_quo stores reference plan and district count", {
    constr <- redist_constr(iowa_map) %>%
        add_constr_status_quo(strength = 1.5, current = cd_2010)

    expect_s3_class(constr, "redist_constr")
    expect_named(constr, "status_quo")
    el <- constr$status_quo[[1]]
    expect_equal(el$strength, 1.5)
    expect_equal(el$current, iowa_map$cd_2010)
    expect_equal(el$n_current, 4)

    # picks up the map's existing plan when `current` is omitted
    constr2 <- add_constr_status_quo(redist_constr(iowa_map), strength = 1)
    expect_equal(constr2$status_quo[[1]]$current, iowa_map$cd_2010)
})

test_that("add_constr_status_quo errors without an existing plan", {
    map_no_ref <- redist_map(iowa, ndists = 4, pop_tol = 0.05)
    expect_error(
        add_constr_status_quo(redist_constr(map_no_ref), strength = 1),
        "current"
    )
})

test_that("add_constr_grp_pow stores group/total/tgt/pow and uses default pop", {
    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_pow(
            strength = 2, bvap, vap,
            tgt_group = 0.55, tgt_other = 0.25, pow = 1.5
        )

    el <- constr$grp_pow[[1]]
    expect_equal(el$strength, 2)
    expect_equal(el$group_pop, iowa_map$bvap)
    expect_equal(el$total_pop, iowa_map$vap)
    expect_equal(el$tgt_group, 0.55)
    expect_equal(el$tgt_other, 0.25)
    expect_equal(el$pow, 1.5)

    # total_pop default falls back to the map's pop column
    constr2 <- add_constr_grp_pow(redist_constr(iowa_map), 1, bvap)
    expect_equal(constr2$grp_pow[[1]]$total_pop, iowa_map$pop)

    expect_error(
        add_constr_grp_pow(redist_constr(iowa_map), 1, bvap[-1]),
        "length"
    )
})

test_that("add_constr_grp_hinge stores tgts_group", {
    constr <- add_constr_grp_hinge(
        redist_constr(iowa_map), 5, bvap, vap,
        tgts_group = c(0.4, 0.5, 0.6)
    )
    el <- constr$grp_hinge[[1]]
    expect_equal(el$tgts_group, c(0.4, 0.5, 0.6))
    expect_equal(el$group_pop, iowa_map$bvap)
})

test_that("add_constr_grp_inv_hinge stores complement group_pop", {
    constr <- add_constr_grp_inv_hinge(
        redist_constr(iowa_map), 5, bvap, vap, tgts_group = 0.5
    )
    el <- constr$grp_inv_hinge[[1]]
    expect_equal(el$tgts_group, 0.5)
    # the inverse-hinge constructor inverts the group, so group_pop = total - group
    expect_equal(el$group_pop, iowa_map$vap - iowa_map$bvap)
    expect_equal(el$total_pop, iowa_map$vap)
})

test_that("add_constr_compet stores dvote/rvote/pow", {
    constr <- add_constr_compet(
        redist_constr(iowa_map), 1, dem_08, rep_08, pow = 0.75
    )
    el <- constr$compet[[1]]
    expect_equal(el$dvote, iowa_map$dem_08)
    expect_equal(el$rvote, iowa_map$rep_08)
    expect_equal(el$pow, 0.75)

    expect_error(
        add_constr_compet(redist_constr(iowa_map), 1, dem_08[-1], rep_08),
        "length"
    )
})

test_that("add_constr_incumbency stores incumbent indices", {
    constr <- add_constr_incumbency(
        redist_constr(iowa_map), 0.5, incumbents = c(1L, 50L, 99L)
    )
    expect_equal(constr$incumbency[[1]]$incumbents, c(1L, 50L, 99L))
})

test_that("add_constr_splits groups admin and stores counts", {
    constr <- add_constr_splits(redist_constr(iowa_map), 1, admin = region)
    el <- constr$splits[[1]]
    expect_equal(length(el$admin), nrow(iowa_map))
    expect_equal(el$n, length(unique(iowa_map$region)))
    # admin is reduced to consecutive ids 1..n
    expect_equal(sort(unique(el$admin)), seq_len(el$n))

    expect_error(
        add_constr_splits(redist_constr(iowa_map), 1, admin = NULL),
        "NULL"
    )

    bad <- iowa_map$region
    bad[1] <- NA
    expect_error(
        add_constr_splits(redist_constr(iowa_map), 1, admin = bad),
        "NA"
    )
})

test_that("add_constr_multisplits and total_splits mirror splits structure", {
    c1 <- add_constr_multisplits(redist_constr(iowa_map), 1, admin = region)
    c2 <- add_constr_total_splits(redist_constr(iowa_map), 1, admin = region)
    expect_equal(c1$multisplits[[1]]$n, length(unique(iowa_map$region)))
    expect_equal(c2$total_splits[[1]]$n, length(unique(iowa_map$region)))

    bad <- iowa_map$region; bad[1] <- NA
    expect_error(
        add_constr_multisplits(redist_constr(iowa_map), 1, admin = bad),
        "NA"
    )
    expect_error(
        add_constr_total_splits(redist_constr(iowa_map), 1, admin = bad),
        "NA"
    )
})

test_that("add_constr_segregation requires group_pop and stores total_pop", {
    constr <- add_constr_segregation(
        redist_constr(iowa_map), 1, group_pop = bvap
    )
    el <- constr$segregation[[1]]
    expect_equal(el$group_pop, iowa_map$bvap)
    expect_equal(el$total_pop, iowa_map$pop)

    expect_error(
        add_constr_segregation(redist_constr(iowa_map), 1),
        "group_pop"
    )
})

test_that("add_constr_polsby stores area and perimeter info", {
    constr <- add_constr_polsby(redist_constr(iowa_map), 1)
    el <- constr$polsby[[1]]
    expect_equal(length(el$area), nrow(iowa_map))
    expect_true(length(el$from) > 0)
    expect_equal(length(el$from), length(el$to))
    expect_equal(length(el$from), length(el$perimeter))
})

test_that("add_constr_polsby errors on non-sf data", {
    map_no_sf <- redist_map(pop = rep(1, 5), ndists = 2, pop_tol = 0.5,
                            adj = list(1L, c(0L, 2L), c(1L, 3L), c(2L, 4L), 3L))
    expect_error(
        add_constr_polsby(redist_constr(map_no_sf), 1),
        "sf"
    )
})

test_that("add_constr_fry_hold accepts user ssdmat and denominator", {
    ssd <- diag(nrow(iowa_map))
    constr <- add_constr_fry_hold(
        redist_constr(iowa_map), 1,
        ssdmat = ssd, denominator = 10
    )
    el <- constr$fry_hold[[1]]
    expect_equal(el$ssdmat, ssd)
    expect_equal(el$denominator, 10)
    expect_equal(el$total_pop, iowa_map$pop)
})

test_that("add_constr_edges_rem records strength only", {
    constr <- add_constr_edges_rem(redist_constr(iowa_map), 1.25)
    # note: stored under name `edges_removed`, not `edges_rem`
    expect_equal(constr$edges_removed[[1]]$strength, 1.25)
})

test_that("add_constr_qps records strength, cities, and n_cty", {
    cities <- as.integer(seq_len(nrow(iowa_map)) %% 3)
    constr <- suppressMessages(redist:::add_constr_qps(
        redist_constr(iowa_map), 1, cities = cities
    ))
    el <- constr$qps[[1]]
    expect_equal(el$cities, cities)
    expect_equal(el$n_cty, max(cities) + 1)
})

test_that("constraints can be stacked under the same name", {
    constr <- redist_constr(iowa_map) %>%
        add_constr_grp_hinge(1, bvap, vap, tgts_group = 0.5) %>%
        add_constr_grp_hinge(2, hvap, vap, tgts_group = 0.4)
    expect_length(constr$grp_hinge, 2)
    expect_equal(constr$grp_hinge[[1]]$strength, 1)
    expect_equal(constr$grp_hinge[[2]]$strength, 2)
})

test_that("constraint constructors reject non-redist_constr input", {
    expect_error(add_constr_status_quo(list(), 1, 1:3), "redist_constr")
    expect_error(add_constr_grp_pow(list(), 1, 1, 1), "redist_constr")
    expect_error(add_constr_grp_hinge(list(), 1, 1, 1), "redist_constr")
    expect_error(add_constr_grp_inv_hinge(list(), 1, 1, 1), "redist_constr")
    expect_error(add_constr_compet(list(), 1, 1, 1), "redist_constr")
    expect_error(add_constr_incumbency(list(), 1, 1), "redist_constr")
    expect_error(add_constr_splits(list(), 1, 1), "redist_constr")
    expect_error(add_constr_multisplits(list(), 1, 1), "redist_constr")
    expect_error(add_constr_total_splits(list(), 1, 1), "redist_constr")
    expect_error(add_constr_segregation(list(), 1, 1), "redist_constr")
    expect_error(add_constr_polsby(list(), 1), "redist_constr")
    expect_error(add_constr_fry_hold(list(), 1), "redist_constr")
    expect_error(add_constr_edges_rem(list(), 1), "redist_constr")
    expect_error(add_constr_pop_dev(list(), 1), "redist_constr")
})

test_that("nonpositive strength warns for typical constraints", {
    expect_warning(
        add_constr_pop_dev(redist_constr(iowa_map), 0),
        "Nonpositive"
    )
    expect_warning(
        add_constr_status_quo(redist_constr(iowa_map), -1, cd_2010),
        "Nonpositive"
    )
    expect_warning(
        add_constr_splits(redist_constr(iowa_map), 0, admin = region),
        "Nonpositive"
    )
})
