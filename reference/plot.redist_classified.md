# Plot a plan classification

Plot a plan classification

## Usage

``` r
# S3 method for class 'redist_classified'
plot(x, plans, shp, type = "fill", which = NULL, ...)
```

## Arguments

- x:

  a `redist_classified` object, the output of
  [`classify_plans()`](http://alarm-redist.org/redist/reference/classify_plans.md).

- plans:

  a
  [redist_plans](http://alarm-redist.org/redist/reference/redist_plans.md)
  object.

- shp:

  a shapefile or
  [redist_map](http://alarm-redist.org/redist/reference/redist_map.md)
  object.

- type:

  either `"line"` or `"fill"`. Passed on to
  [`compare_plans()`](http://alarm-redist.org/redist/reference/compare_plans.md)
  as `plot`.

- which:

  indices of the splits to plot. Defaults to all

- ...:

  passed on to
  [`compare_plans()`](http://alarm-redist.org/redist/reference/compare_plans.md)

## Value

ggplot comparison plot
