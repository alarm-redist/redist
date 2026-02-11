# Combine scoring functions

`redist_scorer` functions may be combined together to optimize along
multiple dimensions. Rather than linearly combining multiple scorers to
form a single objective as with
[scorer-arith](http://alarm-redist.org/redist/reference/scorer-arith.md),
these functions allow analysts to approximate the Pareto frontier for a
set of scorers.

## Usage

``` r
combine_scorers(...)

# S3 method for class 'redist_scorer'
cbind(..., deparse.level = 1)
```

## Arguments

- ...:

  a numeric or a `redist_scorer` function, from
  [`scorers`](http://alarm-redist.org/redist/reference/scorers.md)

- deparse.level:

  As in [`cbind()`](https://rdrr.io/r/base/cbind.html).

## Value

function of class redist_scorer. Will return a matrix with each column
containing every plan's scores for a particular scoring function.
