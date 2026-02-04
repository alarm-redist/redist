devtools::load_all()
library(sf)

# constants ----
n <- 7L
n_chains <- 12L

# targets from Philip ----
target <- structure(list(edges = c(
  28, 29, 30, 31, 32, 33, 34, 35, 36,
  37, 38, 39, 40, 41, 42
), true_prop = c(
  0.020807687768601, 0.0710143952315672,
  0.14977414247499, 0.197857201120775, 0.2074320577516, 0.163717974693393,
  0.103810429703219, 0.0529758100529919, 0.0222509101447097, 0.0076427174466794,
  0.00214991687533369, 0.000477190469004618, 7.99810543944961e-05,
  9.03290674648099e-06, 5.52305514866464e-07
)), class = c(
  'tbl_df',
  'tbl', 'data.frame'
), row.names = c(NA, -15L))

# make shapes ----
shp <- geomander::checkerboard |>
  filter(i < n, j < n) |>
  st_set_crs(3857) |>
  mutate(pop = n)
map <- redist_map(shp, ndists = n, pop_tol = 0.01)

set.seed(n)
init_plans <- redist_smc(map, nsims = 100) |>
  get_plans_matrix()

plans <- redist_cyclewalk(
  map = map,
  nsims = 5e6,
  chains = n_chains,
  thin = 500,
  warmup = 5e5,
  init_plan = init_plans[, seq_len(n_chains)]
)

plans <- plans |>
  mutate(
    edges_rem = comp_edges_rem(pl(), map)
  )
