# Prepares a run of the enumpart algorithm by ordering edges

Prepares a run of the enumpart algorithm by ordering edges

## Usage

``` r
redist.prep.enumpart(
  adj,
  unordered_path,
  ordered_path,
  weight_path = NULL,
  total_pop = NULL
)
```

## Arguments

- adj:

  zero indexed adjacency list

- unordered_path:

  valid path to output the unordered adjacency map to

- ordered_path:

  valid path to output the ordered adjacency map to

- weight_path:

  A path (not including ".dat") to store a space-delimited file
  containing a vector of vertex weights. Only supply with total_pop.

- total_pop:

  the vector of precinct populations. Only supply with weight_path

## Value

0 on success

## References

Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
"The Essential Role of Empirical Validation in Legislative Redistricting
Simulation." Forthcoming, Statistics and Public Policy.

## Examples

``` r
if (FALSE) { # \dontrun{
temp <- tempdir()
data(fl25)
adj <- redist.adjacency(fl25)
redist.prep.enumpart(adj = adj, unordered_path = paste0(temp, "/unordered"),
    ordered_path = paste0(temp, "/ordered"))
} # }
```
