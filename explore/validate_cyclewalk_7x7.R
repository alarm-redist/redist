devtools::load_all()
library(sf)
library(dplyr)

# constants ----
n <- 7L
n_chains <- 8L

# targets from Philip ----
target <- structure(
  list(
    edges = c(
      28,
      29,
      30,
      31,
      32,
      33,
      34,
      35,
      36,
      37,
      38,
      39,
      40,
      41,
      42
    ),
    true_prop = c(
      0.020807687768601,
      0.0710143952315672,
      0.14977414247499,
      0.197857201120775,
      0.2074320577516,
      0.163717974693393,
      0.103810429703219,
      0.0529758100529919,
      0.0222509101447097,
      0.0076427174466794,
      0.00214991687533369,
      0.000477190469004618,
      7.99810543944961e-05,
      9.03290674648099e-06,
      5.52305514866464e-07
    )
  ),
  class = c(
    'tbl_df',
    'tbl',
    'data.frame'
  ),
  row.names = c(NA, -15L)
)

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
  nsims = 2e6,
  chains = n_chains,
  ncores = n_chains / 2,
  thin = 100,
  warmup = 1e6,
  compactness = 0,
  init_plan = init_plans[, seq_len(n_chains)]
)

plans <- plans |>
  mutate(
    edges_rem = comp_edges_rem(pl(), map)
  )

res <- tibble(
  edges = sort(unique(plans$edges_rem)),
  est_prop = as.vector(
    table(
      plans |> subset_sampled() |> pull(edges_rem)
    ) /
      nrow(subset_sampled(plans))
  )
)

target |>
  full_join(res, by = 'edges') |>
  mutate(
    est_prop = tidyr::replace_na(est_prop, 0),
    diff = est_prop - true_prop
  )


library(ggplot2)
make_valid_plot <- function(pl_cw, target, conf = 0.9, by_chain = TRUE) {
  PAL <- c('#CC79A7', '#E69F00', '#56B4E9', '#20B073', '#0072B2', '#D55E00', '#999999')
  # edge cut compactness
  edge_rg <- 28:42
  d_hist <- purrr::map_dfr(edge_rg, function(k) {
    redist_ci(pl_cw, edges_rem == k, 1, conf = conf, by_chain = by_chain) %>%
      suppressWarnings() %>%
      `colnames<-`(c('est', 'low', 'high')) %>%
      mutate(
        edges = k,
        true = target$true_prop[target$edges == k]
      )
  })

  d_hist |>
    ggplot() +
    geom_col(aes(x = edges, y = true), fill = PAL[1]) +
    geom_pointrange(aes(x = edges, y = est, ymin = low, ymax = high),
      color = '#444444',
      linewidth = 0.7, size = 1.5
    ) +
    scale_y_continuous(str_glue('Proportion of plans'),
      labels = scales::label_percent(),
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = str_glue('{scales::label_comma()(1e6)} samples per run'), x = 'Number of edges removed') +
    theme_bw(base_family = 'Times', base_size = 10)
}

make_valid_plot(plans, target)
