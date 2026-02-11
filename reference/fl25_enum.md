# All Partitions of 25 Precincts into 3 Congressional Districts (No Population Constraint)

This data set contains demographic and geographic information about 25
contiguous precincts in the state of Florida. The data lists all
possible partitions of the 25 precincts into three contiguous
congressional districts. The 25-precinct shapefile may be found in
[`fl25`](http://alarm-redist.org/redist/reference/fl25.md)

## Usage

``` r
data("fl25_enum")
```

## Format

A list with two entries:

- `plans`:

  A matrix containing every partition of the 25 precincts into three
  contiguous congressional districts, with no population constraint.

- `pop_dev`:

  A vector containing the maximum population deviation across the three
  districts for each plan.

## References

Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
(2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
Carlo." Working Paper. Available at
<http://imai.princeton.edu/research/files/redist.pdf>.

Massey, Douglas and Nancy Denton. (1987) "The Dimensions of Social
Segregation". Social Forces.

## Examples

``` r
data(fl25_enum)
```
