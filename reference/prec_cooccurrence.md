# Compute a matrix of precinct co-occurrences

For a map with `n` precincts Returns an `n`-by-`n` matrix, where each
entry measures the fraction of the plans in which the row and column
precincts were in the same district.

## Usage

``` r
prec_cooccurrence(plans, which = NULL, sampled_only = TRUE, ncores = 1)
```

## Arguments

- plans:

  a
  [redist_plans](http://alarm-redist.org/redist/reference/redist_plans.md)
  object.

- which:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  which plans to compute the co-occurrence over. Defaults to all.

- sampled_only:

  if `TRUE`, do not include reference plans.

- ncores:

  the number of parallel cores to use in the computation.

## Value

a symmetric matrix the size of the number of precincts.
