# Counts the Number of Counties within a District

Counts the total number of counties that are found within a district.
This does not subtract out the number of counties that are found
completely within a district.

## Usage

``` r
redist.district.splits(plans, counties)
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
#old redist.district.splits(plans, ia$region)
splits_count(plans, ia, region)
#>           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> South        4    4    4    3    4    3    3    3    3     3     3     4     4
#> Northeast    2    2    2    1    2    3    3    2    2     2     2     2     2
#> Northwest    2    2    2    2    1    2    2    3    3     3     3     2     2
#> Southeast    3    2    2    3    3    3    3    3    3     2     3     2     2
#> Central      2    2    2    3    2    2    2    2    2     3     2     2     2
#>           [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> South         3     3     3     3     3     3     3     4     3     4     3
#> Northeast     2     2     2     2     4     2     2     2     3     1     2
#> Northwest     3     3     2     3     1     2     2     2     3     2     2
#> Southeast     3     2     2     2     3     3     3     2     3     3     3
#> Central       2     3     3     3     2     2     2     2     2     3     2
#>           [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
#> South         3     3     4     2     2     4     3     3     2     2     2
#> Northeast     2     2     2     2     3     2     1     1     3     2     2
#> Northwest     2     3     2     3     3     2     2     2     2     3     3
#> Southeast     3     3     2     4     3     3     3     3     3     3     3
#> Central       2     2     2     2     2     2     2     2     2     2     2
#>           [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46]
#> South         4     3     3     3     3     3     4     3     3     2     4
#> Northeast     2     2     2     2     2     3     2     2     2     3     2
#> Northwest     2     2     2     2     2     2     2     2     3     3     2
#> Southeast     3     2     3     3     3     3     2     3     3     2     2
#> Central       2     2     2     2     2     2     2     3     2     2     2
#>           [,47] [,48] [,49] [,50] [,51]
#> South         4     4     3     3     4
#> Northeast     2     2     2     3     3
#> Northwest     2     2     3     3     2
#> Southeast     2     2     3     3     2
#> Central       2     2     2     2     2
```
