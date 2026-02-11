# Counts the Number of Counties Split Between 3 or More Districts

Counts the total number of counties that are split across more than 2
districts.

## Usage

``` r
redist.multisplits(plans, counties)
```

## Arguments

- plans:

  A numeric vector (if only one map) or matrix with one row for each
  precinct and one column for each map. Required.

- counties:

  A vector of county names or county ids.

## Value

integer matrix where each district is a

## Examples

``` r
data(iowa)
ia <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop, pop_tol = 0.01)
plans <- redist_smc(ia, 50, silent = TRUE)
#old redist.multisplits(plans, ia$region)
splits_multi(plans, ia, region)
#>   [1] 2 2 2 2 4 4 4 4 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3
#>  [38] 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3
#>  [75] 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2
#> [112] 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3
#> [149] 3 3 3 3 5 5 5 5 3 3 3 3 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3
#> [186] 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2
```
