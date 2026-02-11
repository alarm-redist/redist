# Florida 25 Precinct File

This data set contains the 25-precinct shapefile and related data for
each precinct. All possible partitions of the 25 precincts into three
contiguous congressional districts are stored in
[`fl25_enum`](http://alarm-redist.org/redist/reference/fl25_enum.md),
and the corresponding adjacency graph is stored in `fl25_adj`.

## Format

A list storing the adjacency graph for the 25-precinct subset of
Florida.

## References

Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
(2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
Carlo." Working Paper. Available at
<http://imai.princeton.edu/research/files/redist.pdf>.

## Examples

``` r
data(fl25_adj)
```
