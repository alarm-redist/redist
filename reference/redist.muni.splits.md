# Counts the Number of Municipalities Split Between Districts

Counts the total number of municpalities that are split. Municipalities
in this interpretation do not need to cover the entire state, which
differs from counties.

## Usage

``` r
muni_splits(map, munis, .data = cur_plans())

redist.muni.splits(plans, munis)
```

## Arguments

- map:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- munis:

  A vector of municipality names or ids.

- .data:

  a
  [`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
  object

- plans:

  A numeric vector (if only one map) or matrix with one row for each
  precinct and one column for each map. Required.

## Value

integer vector of length ndist by ncol(plans)

## Examples

``` r
data(iowa)
ia <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop, pop_tol = 0.01)
plans <- redist_smc(ia, 50, silent = TRUE)
ia$region[1:10] <- NA
#old redist.muni.splits(plans, ia$region)
splits_sub_admin(plans, ia, region)
#>   [1] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4
#>  [38] 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 5 5
#>  [75] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
#> [112] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4 5 5 5 5 4 4 4 4 5 5 5 5
#> [149] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
#> [186] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4
```
