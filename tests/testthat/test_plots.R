test_that("redist.plot.map works", {
  out <- redist.plot.map(shp = iowa, plan = iowa$cd_2010)
  expect_true('ggplot' %in% class(out))

  iowa_map = redist_map(iowa, existing_plan = cd_2010, pop_tol=0.01)

  out <- iowa_map  %>% redist.plot.map(shp = ., plan = get_existing(.))
  expect_true('ggplot' %in% class(out))

  out <- iowa_map %>% redist.plot.map(shp = ., plan = cd_2010)
  expect_true('ggplot' %in% class(out))

  out <- iowa_map %>% redist.plot.map(shp = ., plan = cd_2010, fill = white)
  expect_true('ggplot' %in% class(out))
})


test_that("redist.plot.adj works", {
  out <- redist.plot.adj(shp = iowa, plan = iowa$cd_2010)
  expect_true('ggplot' %in% class(out))

  iowa_map = redist_map(iowa, existing_plan = cd_2010, pop_tol=0.01)

  out <- iowa_map %>% redist.plot.adj(shp = ., plan = get_existing(.))
  expect_true('ggplot' %in% class(out))

  out <- iowa_map %>% redist.plot.map(shp = ., plan = cd_2010)
  expect_true('ggplot' %in% class(out))
})
