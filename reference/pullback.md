# Pull back plans to unmerged units

Merging map units through
[`merge_by`](http://alarm-redist.org/redist/reference/merge_by.md) or
[dplyr::summarize](https://dplyr.tidyverse.org/reference/summarise.html)
changes the indexing of each unit. Use this function to take a set of
redistricting plans from a `redist` algorithm and re-index them to be
compatible with the original set of units.

## Usage

``` r
pullback(plans, map = NULL)
```

## Arguments

- plans:

  a `redist_plans` object

- map:

  optionally, a `redist_map` object, which will be used to set the new
  population vector

## Value

a new, re-indexed, `redist_plans` object
