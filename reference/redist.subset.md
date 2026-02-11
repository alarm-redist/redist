# Subset a shp

Subsets a shp object along with its adjacency. Useful for running
smaller analyses on pairs of districts. Provide population, ndists,
pop_tol, and sub_ndists to get proper population parity constraints on
subsets.

## Usage

``` r
redist.subset(shp, adj, keep_rows, total_pop, ndists, pop_tol, sub_ndists)
```

## Arguments

- shp:

  An sf object

- adj:

  A zero-indexed adjacency list. Created with `redist.adjacency` if not
  supplied.

- keep_rows:

  row numbers of precincts to keep. Random submap selected if not
  supplied.

- total_pop:

  numeric vector with one entry for the population of each precinct.

- ndists:

  integer, number of districts in whole map

- pop_tol:

  The strength of the hard population constraint.

- sub_ndists:

  integer, number of districts in subset map

## Value

a list containing the following components:

- shp:

  The subsetted shp object

- adj:

  The subsetted adjacency list for shp

- keep_rows:

  The indices of the rows kept.

- sub_ndists:

  The number of districts in the subset.

- sub_pop_tol:

  The new parity constraint for a subset.
