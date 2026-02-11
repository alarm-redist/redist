# Runs the enumpart algorithm

Runs the enumpart algorithm

## Usage

``` r
redist.run.enumpart(
  ordered_path,
  out_path,
  ndists = 2,
  all = TRUE,
  n = NULL,
  weight_path = NULL,
  lower = NULL,
  upper = NULL,
  options = NULL
)
```

## Arguments

- ordered_path:

  Path used in redist.prep.enumpart (not including ".dat")

- out_path:

  Valid path to output the enumerated districts

- ndists:

  number of districts to enumerate

- all:

  boolean. TRUE outputs all districts. FALSE samples n districts.

- n:

  integer. Number of districts to output if all is FALSE. Returns
  districts selected from uniform random distribution.

- weight_path:

  A path (not including ".dat") to a space-delimited file containing a
  vector of vertex weights, to be used along with `lower` and `upper`.

- lower:

  A lower bound on each partition's total weight, implemented by
  rejection sampling.

- upper:

  An upper bound on each partition's total weight.

- options:

  Additional enumpart arguments. Not recommended for use.

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
redist.run.enumpart(ordered_path = paste0(temp, "/ordered"),
    out_path = paste0(temp, "/enumerated"))
} # }
```
