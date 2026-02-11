# Plot a district assignment

Plot a district assignment

## Usage

``` r
redist.plot.plans(
  plans,
  draws,
  shp,
  qty = NULL,
  interactive = FALSE,
  ...,
  geom = NULL
)
```

## Arguments

- plans:

  a `redist_plans` object.

- draws:

  the plan(s) to plot. Will match the `draw` column of `x`.

- qty:

  the quantity to plot. Defaults to the district assignment.

- interactive:

  if `TRUE`, show an interactive map in the viewer rather than a static
  map. Only uses the first element of `draws`

- ...:

  additional arguments passed to the plotting functions.

- geom, shp:

  the `redist_map` geometry to use (`geom` is deprecated).

## Value

A ggplot

## Examples

``` r
library(dplyr)
data(iowa)

iowa <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)
plans <- redist_smc(iowa, nsims = 100, silent = TRUE)
redist.plot.plans(plans, c(1, 2, 3, 4), iowa)

```
