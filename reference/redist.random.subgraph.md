# Return a random subgraph of a shape

`random.subgraph` returns a random subset of the shp provided

## Usage

``` r
redist.random.subgraph(shp, n, adj = NULL)
```

## Arguments

- shp:

  sf object or SpatialPolygonsDataFrame

- n:

  number of edges to sample. n must be a positive integer.

- adj:

  Optional. zero indexed adjacency list.

## Value

sf dataframe with n rows

## Details

Snowball sampling with backtracking
