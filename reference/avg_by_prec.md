# Average a variable by precinct (Deprecated)

Deprecated in favor of
[`proj_avg()`](http://alarm-redist.org/redist/reference/proj.md). Takes
a column of a `redist_plans` object and averages it across a set of
`draws` for each precinct.

## Usage

``` r
avg_by_prec(plans, x, draws = NA)
```

## Arguments

- plans:

  a `redist_plans` object

- x:

  an expression to average. Tidy-evaluated in `plans`.

- draws:

  which draws to average. `NULL` will average all draws, including
  reference plans. The special value `NA` will average all sampled
  draws. An integer, logical, or character vector indicating specific
  draws may also be provided.

## Value

a vector of length matching the number of precincts, containing the
average.
