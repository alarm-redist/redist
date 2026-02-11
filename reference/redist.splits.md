# Count County Splits

Count County Splits

## Usage

``` r
county_splits(map, counties, .data = cur_plans())

redist.splits(plans, counties)
```

## Arguments

- map:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- counties:

  A vector of county names or county ids.

- .data:

  a
  [`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
  object

- plans:

  A numeric vector (if only one map) or matrix with one row for each
  precinct and one column for each map. Required.

## Value

integer vector with one number for each map
