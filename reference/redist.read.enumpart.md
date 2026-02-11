# Read Results from enumpart

Read Results from enumpart

## Usage

``` r
redist.read.enumpart(out_path, skip = 0, n_max = -1L)
```

## Arguments

- out_path:

  out_path specified in redist.run.enumpart

- skip:

  number of lines to skip

- n_max:

  max number of lines to read

## Value

district_membership matrix

## References

Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
"The Essential Role of Empirical Validation in Legislative Redistricting
Simulation." Forthcoming, Statistics and Public Policy.

## Examples

``` r
if (FALSE) { # \dontrun{
temp <- tempdir()
cds <- redist.read.enumpart(out_path = paste0(temp, "/enumerated"))
} # }
```
