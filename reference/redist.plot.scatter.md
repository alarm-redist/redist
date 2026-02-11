# Scatter plot of plan summary statistics

Makes a scatterplot of two quantities of interest across districts or
plans.

## Usage

``` r
redist.plot.scatter(plans, x, y, ..., bigger = TRUE)
```

## Arguments

- plans:

  the `redist_plans` object.

- x:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the quantity to plot on the horizontal axis.

- y:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the quantity to plot on the vertical axis.

- ...:

  passed on to
  [`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

- bigger:

  if TRUE, make the point corresponding to the reference plan larger.

## Value

A ggplot

## Examples

``` r
library(dplyr)
data(iowa)

iowa <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)
plans <- redist_smc(iowa, nsims = 100, silent = TRUE)
plans %>%
    mutate(comp = distr_compactness(iowa)) %>%
    group_by(draw) %>%
    summarize(pop_dev = max(abs(total_pop/mean(total_pop) - 1)),
        comp = comp[1]) %>%
    redist.plot.scatter(pop_dev, comp)
#> Warning: There were 2 warnings in `"draw" %in% names(data)`.
#> The first warning was:
#> ℹ In argument: `comp = distr_compactness(iowa)`.
#> Caused by warning in `distr_compactness()`:
#> ! 'distr_compactness' is deprecated.
#> See help("Deprecated")
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

```
