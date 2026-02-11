# Relabel Discontinuous Counties

Relabel Discontinuous Counties

## Usage

``` r
redist.county.relabel(adj, counties, simplify = TRUE)
```

## Arguments

- adj:

  adjacency list

- counties:

  character vector of county names

- simplify:

  boolean - TRUE returns a numeric vector of ids, while FALSE appends a
  number when there are multiple connected components.

## Value

character vector of county names

## Examples

``` r
set.seed(2)
data(fl25)
data(fl25_adj)
counties <- sample(c(rep("a", 20), rep("b", 5)))
redist.county.relabel(fl25_adj, counties)
#>  [1] 2 1 1 2 1 1 1 1 1 1 1 1 1 2 1 1 1 1 3 1 1 1 1 4 1
```
