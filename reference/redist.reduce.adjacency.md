# Reduce Adjacency List

Tool to help reduce adjacency lists for analyzing subsets of maps.

## Usage

``` r
redist.reduce.adjacency(adj, keep_rows)
```

## Arguments

- adj:

  A zero-indexed adjacency list. Required.

- keep_rows:

  row numbers of precincts to keep

## Value

zero indexed adjacency list with max value length(keep_rows) - 1

## Examples

``` r
data(fl25_adj)
redist.reduce.adjacency(fl25_adj, c(2, 3, 4, 6, 21))
#> [[1]]
#> [1] 2 4
#> 
#> [[2]]
#> [1] 2 3 4
#> 
#> [[3]]
#> [1] 0 1 4
#> 
#> [[4]]
#> [1] 1
#> 
#> [[5]]
#> [1] 0 1 2
#> 
```
