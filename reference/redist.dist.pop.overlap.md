# Compare the Population Overlap Across Plans at the District Level

This implements Crespin's 2005 measure of district continuity, as
applied to the geographies represented by a plan, typically precincts or
voting districts. This implementation assumes none of the precincts in
plan_old or plan_new are split.

## Usage

``` r
redist.dist.pop.overlap(plan_old, plan_new, total_pop, normalize_rows = TRUE)
```

## Arguments

- plan_old:

  The reference or original plan to compare against

- plan_new:

  The new plan to compare to the reference plan

- total_pop:

  The total population by precinct This can also take a redist_map
  object and will use the population in that object. If nothing is
  provided, it weights all entries in plan equally.

- normalize_rows:

  Default TRUE. Normalize populations by row. If FALSE, normalizes by
  column. If NULL, does not normalize.

## Value

matrix with length(unique(plan_old)) rows and length(unique(plan_new))
columns

## References

"Using Geographic Information Systems to Measure District Change,
2000-02", Michael Crespin, Political Analysis (2005) 13(3): 253-260

## Examples

``` r
set.seed(5)
data(iowa)
iowa_map <- redist_map(iowa, total_pop = pop, pop_tol = 0.01, ndists = 4)
plans <- redist_smc(iowa_map, 2)
#> SEQUENTIAL MONTE CARLO
#> Sampling 2 99-unit maps with 4 districts and population between 753973 and 769205.
plans_mat <- get_plans_matrix(plans)
ov <- redist.dist.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map)
round(ov, 2)
#>      1    2    3    4
#> 1 0.14 0.00 0.18 0.68
#> 2 0.05 0.17 0.78 0.00
#> 3 0.59 0.06 0.02 0.33
#> 4 0.22 0.76 0.02 0.00

ov_col <- redist.dist.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map, normalize_rows = FALSE)
round(ov_col, 2)
#>      1    2    3    4
#> 1 0.14 0.00 0.18 0.67
#> 2 0.05 0.17 0.78 0.00
#> 3 0.59 0.06 0.02 0.33
#> 4 0.22 0.77 0.02 0.00

ov_un_norm <- redist.dist.pop.overlap(plans_mat[, 1], plans_mat[, 2],
    iowa_map, normalize_rows = NULL)
round(ov_un_norm, 2)
#>        1      2      3      4
#> 1 105215      0 138897 513848
#> 2  40648 131090 592727      0
#> 3 446353  43843  14928 252869
#> 4 165224 582584  18129      0

iowa_map_5 <- iowa_map <- redist_map(iowa, total_pop = pop, pop_tol = 0.01, ndists = 5)
plan_5 <- get_plans_matrix(redist_smc(iowa_map_5, 1))
#> SEQUENTIAL MONTE CARLO
#> Sampling 1 99-unit maps with 5 districts and population between 603178 and 615364.
ov4_5 <- redist.dist.pop.overlap(plans_mat[, 1], plan_5, iowa_map)
round(ov4_5, 2)
#>      1    2    3    4    5
#> 1 0.00 0.00 0.70 0.19 0.11
#> 2 0.00 0.25 0.00 0.37 0.39
#> 3 0.33 0.02 0.11 0.24 0.30
#> 4 0.47 0.53 0.00 0.00 0.00
```
