library(tidyverse)
devtools::load_all()

data(fl25)
data(fl25_enum)

# Set up redistricting problem
map <- redist_map(fl25, pop_tol = 0.1, ndists = 3)

set.seed(08450)
plans_sim <- redist_mew(map,
  nsims = 1e5, chains = 8L, thin = 100,
  warmup = 1e4,
  compactness = 1
) |>
  mutate(
    comp_edges = comp_edges_rem(pl(), map),
    dshare = group_frac(map, obama, obama + mccain)
  ) |>
  number_by(dshare) # convert dshare to order statistics

# Reweight enumerated plans to target distribution
# target is proportional to log_st compactness measure
plans_enum <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.1]
wgt_enum <- comp_log_st(plans_enum, map)
wgt_enum <- exp(wgt_enum - max(wgt_enum))
wgt_enum <- 3 * wgt_enum / sum(wgt_enum) # wgts are repeated for each district
plans_enum <- redist_plans(plans_enum, map, algorithm = 'enumpart') |>
  mutate(
    comp_edges = comp_edges_rem(pl(), map),
    dshare = group_frac(map, obama, obama + mccain),
    .wgt = wgt_enum
  ) |>
  number_by(dshare)

# Summarize sampling estimate and ground truth for each summary stat
## Plan compactness
comp_rg <- sort(unique(plans_enum$comp_edges))

# `redist_smc_ci()` checks R-hat diagnostics and computes interval estimates
d_comp_est <- map(comp_rg, function(bin) {
  redist_smc_ci(plans_sim, comp_edges == bin, conf = 0.95) |>
    rename(est = 1, lower = 2, upper = 3) |>
    mutate(comp_edges = bin, .before = est)
}) |>
  bind_rows()
d_comp_enum <- plans_enum |>
  as_tibble() |>
  filter(district == 1L) |> # this stat is the same across districts
  count(comp_edges, wt = .wgt, name = 'true')

## Avg. Democratic vote share by district
d_dem_est <- map(1:3, function(distr) {
  redist_smc_ci(plans_sim, dshare, district = distr, conf = 0.95) |>
    rename(est = 1, lower = 2, upper = 3) |>
    mutate(district = distr, .before = est)
}) |>
  bind_rows()
d_dem_enum <- plans_enum |>
  as_tibble() |>
  summarize(true = sum(dshare * .wgt), .by = district)

# Produce plots
library(patchwork)

p_comp <- left_join(d_comp_est, d_comp_enum, by = 'comp_edges') |>
  ggplot(aes(comp_edges)) +
  geom_col(aes(y = true), fill = '#E4C22B') +
  geom_pointrange(aes(y = est, ymin = lower, ymax = upper), size = 0.2) +
  scale_y_continuous('Frequency', expand = expansion(mult = c(0, 0.05))) +
  labs(x = 'Cut edges (smaller is more compact)') +
  theme_bw(base_size = 9) +
  theme(plot.margin = margin())
p_distr <- d_dem_est |>
  left_join(d_dem_enum, by = 'district') |>
  mutate(across(est:upper, ~ .x - true)) |>
  ggplot(aes(factor(district))) +
  geom_hline(yintercept = 0, lwd = 1.0, color = '#E4C22B') +
  geom_pointrange(aes(y = est, ymin = lower, ymax = upper), size = 0.2) +
  scale_y_continuous(
    'Error in average Democratic vote share',
    limits = c(-0.01, 0.01),
    labels = scales::label_percent(suffix = 'pp')
  ) +
  labs(x = 'District, ordered by Dem. share') +
  theme_bw(base_size = 9) +
  theme(plot.margin = margin(l = 8))
p_comp +
  p_distr +
  plot_layout(widths = c(5, 3)) +
  plot_annotation(theme = theme(plot.margin = margin()))
