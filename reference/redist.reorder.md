# Reorders district numbers

Ensures that for each column in the plans object, the first district
listed is 1, the second is 2, up to n districts. Assumes that all
columns have the same number of districts as the first.

## Usage

``` r
redist.reorder(plans)
```

## Arguments

- plans:

  A numeric vector (if only one map) or matrix with one row for each
  precinct and one column for each map.

## Value

integer matrix

## Examples

``` r
cds <- matrix(c(rep(c(4L, 5L, 2L, 1L, 3L), 5),
    rep(c(5L, 4L, 3L, 2L, 1L), 2), rep(c(4L, 5L, 2L, 1L, 3L), 3)), nrow = 25)
redist.reorder(cds)
#>       [,1] [,2]
#>  [1,]    1    1
#>  [2,]    2    2
#>  [3,]    3    3
#>  [4,]    4    4
#>  [5,]    5    5
#>  [6,]    1    1
#>  [7,]    2    2
#>  [8,]    3    3
#>  [9,]    4    4
#> [10,]    5    5
#> [11,]    1    2
#> [12,]    2    1
#> [13,]    3    4
#> [14,]    4    5
#> [15,]    5    3
#> [16,]    1    2
#> [17,]    2    1
#> [18,]    3    4
#> [19,]    4    5
#> [20,]    5    3
#> [21,]    1    2
#> [22,]    2    1
#> [23,]    3    4
#> [24,]    4    5
#> [25,]    5    3
```
