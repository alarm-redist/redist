# Freeze Parts of a Map

Freeze Parts of a Map

## Usage

``` r
freeze(freeze_row, plan, .data = cur_map())

redist.freeze(adj, freeze_row, plan = rep(1, length(adj)))
```

## Arguments

- freeze_row:

  Required, logical vector where TRUE freezes and FALSE lets a precinct
  stay free or a vector of indices to freeze

- plan:

  A vector of district assignments, which if provided will create
  separate groups by district. Recommended. In `freeze` defaults to the
  existing plan, if one exists.

- .data:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- adj:

  Required, zero indexed adjacency list.

## Value

integer vector to group by

## Examples

``` r
library(redist)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following object is masked from ‘package:redistmetrics’:
#> 
#>     tally
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
data(fl25)
data(fl25_enum)
data(fl25_adj)
plan <- fl25_enum$plans[, 5118]
freeze_id <- redist.freeze(adj = fl25_adj, freeze_row = (plan == 2),
    plan = plan)

data(iowa)
map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.02)
map <- map %>% merge_by(freeze(cd_2010 == 1, .data = .))
#> Warning: There was 1 warning in `dplyr::summarize()`.
#> ℹ In argument: `dplyr::across(where(is.numeric), sum, na.rm = TRUE)`.
#> ℹ In group 1: `cd_2010 = 1`, `freeze(cd_2010 == 1, .data = .) = 3`.
#> Caused by warning:
#> ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
#> Supply arguments directly to `.fns` through an anonymous function instead.
#> 
#>   # Previously
#>   across(a:b, mean, na.rm = TRUE)
#> 
#>   # Now
#>   across(a:b, \(x) mean(x, na.rm = TRUE))
#> ℹ The deprecated feature was likely used in the redist package.
#>   Please report the issue at <https://github.com/alarm-redist/redist/issues>.
```
