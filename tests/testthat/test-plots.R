test_that("redist.plot.map works", {
  data(iowa)
  
  out <- redist.plot.map(shp = iowa, plan = iowa$cd_2010)
  
  expect_true('ggplot' %in% class(out))
  
  out <- iowa %>% redist_map(existing_plan = cd_2010) %>% redist.plot.map(shp = ., plan = get_existing(.))
  
  expect_true('ggplot' %in% class(out))

  out <- iowa %>% redist_map(existing_plan = cd_2010) %>% redist.plot.map(shp = ., plan = cd_2010)
  
  expect_true('ggplot' %in% class(out))
  
  out <- iowa %>% redist_map(existing_plan = cd_2010) %>% redist.plot.map(shp = ., plan = cd_2010, fill = white)
  
  expect_true('ggplot' %in% class(out))
})


test_that("redist.plot.adj works", {
  data(iowa)
  
  out <- redist.plot.adj(shp = iowa, plan = iowa$cd_2010)
  
  expect_true('ggplot' %in% class(out))
  
  out <- iowa %>% redist_map(existing_plan = cd_2010) %>% redist.plot.adj(shp = ., plan = get_existing(.))
  
  expect_true('ggplot' %in% class(out))
  
  out <- iowa %>% redist_map(existing_plan = cd_2010) %>% redist.plot.map(shp = ., plan = cd_2010)
  
  expect_true('ggplot' %in% class(out))
  
})
