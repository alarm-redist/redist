# Make a traceplot for a summary statistic

For a statistic in a
[`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
object, make a traceplot showing the evolution of the statistic over
MCMC iterations.

## Usage

``` r
redist.plot.trace(plans, qty, district = 1L, ...)
```

## Arguments

- plans:

  the `redist_plans` object.

- qty:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the statistic.

- district:

  for `redist_plans` objects with multiple districts, which `district`
  to subset to for plotting. Set to `NULL` to perform no subsetting.

- ...:

  passed on to
  [`geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html)

## Value

A ggplot

## Examples

``` r
library(dplyr)
data(iowa)

iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
plans <- redist_mergesplit_parallel(iowa_map, nsims = 200, chains = 2, silent = TRUE) %>%
    mutate(dem = group_frac(iowa_map, dem_08, dem_08 + rep_08)) %>%
    number_by(dem)
redist.plot.trace(plans, dem, district = 1)

```
