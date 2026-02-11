# Plot Cores

Plot Cores

## Usage

``` r
redist.plot.cores(shp, plan = NULL, core = NULL, lwd = 2)
```

## Arguments

- shp:

  A SpatialPolygonsDataFrame or sf object. Required.

- plan:

  A numeric vector with one entry for each precinct in shp. Used to
  color the districts. Required.

- core:

  Required. integer vector produced by
  [`redist.identify.cores()`](http://alarm-redist.org/redist/reference/redist.identify.cores.md).

- lwd:

  Line width. Defaults to 2.

## Value

ggplot
