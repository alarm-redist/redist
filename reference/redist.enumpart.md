# Enumerate All Parititions (Fifield et al. 2020)

Single function for standard enumeration analysis, using ZDD methodology
(Fifield, Imai, Kawahara, and Kenny 2020).

## Usage

``` r
redist.enumpart(
  adj,
  unordered_path,
  ordered_path,
  out_path,
  ndists = 2,
  all = TRUE,
  n = NULL,
  weight_path = NULL,
  lower = NULL,
  upper = NULL,
  init = FALSE,
  read = TRUE,
  total_pop = NULL
)
```

## Arguments

- adj:

  zero indexed adjacency list.

- unordered_path:

  valid path to output the unordered adjacency map to

- ordered_path:

  valid path to output the ordered adjacency map to

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

- init:

  Runs redist.init.enumpart. Defaults to false. Should be run on first
  use.

- read:

  boolean. Defaults to TRUE. reads

- total_pop:

  the vector of precinct populations

## Value

List with entries district_membership and parity.

## References

Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020). The
essential role of empirical validation in legislative redistricting
simulation. *Statistics and Public Policy*, 7(1), 52-68.
