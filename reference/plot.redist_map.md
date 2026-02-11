# Plot a `redist_map`

Plot a `redist_map`

## Usage

``` r
# S3 method for class 'redist_map'
plot(x, fill = NULL, by_distr = FALSE, adj = FALSE, ...)
```

## Arguments

- x:

  the `redist_map` object

- fill:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  If present, will be used to color the map units. If using data
  masking, may need to explicitly name argument `fill=...` in
  non-interactive contexts to avoid S3 generic issues.

- by_distr:

  if `TRUE` and `fill` is not missing and, color by district and
  indicate the `fill` variable by shading.

- adj:

  if `TRUE`, force plotting the adjacency graph. Overrides `by_distr`.

- ...:

  passed on to
  [`redist.plot.map`](http://alarm-redist.org/redist/reference/redist.plot.map.md)
  (or
  [`redist.plot.adj`](http://alarm-redist.org/redist/reference/redist.plot.adj.md)
  if `adj=TRUE`). Useful parameters may include `zoom_to`, `boundaries`,
  and `title`.

## Value

ggplot2 object

## Examples

``` r
data(fl25)
d <- redist_map(fl25, ndists = 3, pop_tol = 0.05)
#> Projecting to CRS 3857
plot(d)

plot(d, BlackPop/pop)


data(fl25_enum)
fl25$dist <- fl25_enum$plans[, 5118]
d <- redist_map(fl25, existing_plan = dist)
#> Projecting to CRS 3857
plot(d)

```
