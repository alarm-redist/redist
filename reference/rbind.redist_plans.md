# Combine multiple sets of redistricting plans

Only works when all the sets are compatible—generated from the same map,
with the same number of districts. Sets of plans will be indexed by the
`chain` column.

## Usage

``` r
# S3 method for class 'redist_plans'
rbind(..., deparse.level = 1)
```

## Arguments

- ...:

  The
  [`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
  objects to combine. If named arguments are provided, the names will be
  used in the `chain` column; otherwise, numbers will be used for the
  `chain` column.

- deparse.level:

  Ignored.

## Value

A new
[`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
object.
