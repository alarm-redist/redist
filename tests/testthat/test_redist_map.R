
test_that("basic constructor works", {
    iowa_map = redist_map(iowa, ndists=4, pop_tol=0.01)
    expect_s3_class(iowa_map, "redist_map")
    expect_equal(attr(iowa_map, "pop_col"), "pop")
    target = sum(iowa$pop) / 4
    expect_equal(attr(iowa_map, "pop_bounds"), c(0.99, 1.0, 1.01)*target)
    expect_type(get_adj(iowa_map), "list")
})

test_that("constructor works with existing plan", {
    iowa_map = redist_map(iowa, existing_plan = cd_2010, pop_tol=0.01)
    expect_equal(attr(iowa_map, "ndists"), 4)

    expect_message(redist_map(iowa, existing_plan = cd_2010))
})

test_that("filtering works", {
    iowa_map = redist_map(iowa, existing_plan = cd_2010, pop_tol=0.01)
    filtered = filter(iowa_map, cd_2010 %in% c(2, 3))
    bounds = attr(iowa_map, "pop_bounds")
    bounds_f = attr(filtered, "pop_bounds")

    expect_equal(bounds[1], bounds_f[1])
    expect_lt(bounds[2], bounds_f[2])
    expect_equal(bounds[3], bounds_f[3])

    expect_equal(attr(filtered, "ndists"), 2)
    expect_equal(nrow(filtered), max(sapply(get_adj(filtered), max))+1)
})

test_that("get_pop_tol works", {
  iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol=0.02)
  tol <- get_pop_tol(iowa_map)

  expect_equal(tol, 0.02)
})

test_that("set_pop_tol works", {
  iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol=0.02)
  iowa_map <- set_pop_tol(iowa_map, 0.05)

  tol <- get_pop_tol(iowa_map)

  expect_equal(tol, 0.05)
})
