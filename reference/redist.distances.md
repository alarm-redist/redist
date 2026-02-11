# Compute Distance between Partitions

Compute Distance between Partitions

## Usage

``` r
plan_distances(plans, measure = "variation of information", ncores = 1)

redist.distances(plans, measure = "Hamming", ncores = 1, total_pop = NULL)
```

## Arguments

- plans:

  A matrix with one row for each precinct and one column for each map.
  Required.

- measure:

  String vector indicating which distances to compute. Implemented
  currently are "Hamming", "Manhattan", "Euclidean", and "variation of
  information", Use "all" to return all implemented measures. Not case
  sensitive, and any unique substring is enough, e.g. "ham" for Hamming,
  or "info" for variation of information.

- ncores:

  Number of cores to use for parallel computing. Default is 1.

- total_pop:

  The vector of precinct populations. Used only if computing variation
  of information. If not provided, equal population of precincts will be
  assumed, i.e. the VI will be computed with respect to the precincts
  themselves, and not the population.

## Value

`distance_matrix` returns a numeric distance matrix for the chosen
metric.

a named list of distance matrices, one for each distance measure
selected.

## Details

Hamming distance measures the number of different precinct assignments
between plans. Manhattan and Euclidean distances are the 1- and 2-norms
for the assignment vectors. All three of the Hamming, Manhattan, and
Euclidean distances implemented here are not invariant to permutations
of the district labels; permuting will cause large changes in measured
distance, and maps which are identical up to a permutation may be
computed to be maximally distant.

Variation of Information is a metric on population partitions (i.e.,
districtings) which is invariant to permutations of the district labels,
and arises out of information theory. It is calculated as \$\$ VI(\xi,
\xi') = -\sum\_{i=1}^n\sum\_{j=1}^n pop(\xi_i \cap \xi'\_j)/P
(2log(pop(\xi_i \cap \xi'\_j)) - log(pop(\xi_i)) - log(pop(\xi'\_j)))
\$\$ where \\\xi,\xi'\\ are the partitions, \\\xi_i,\xi_j\\ the
individual districts, \\pop(\cdot)\\ is the population, and \\P\\ the
total population of the state. VI is also expressible as the difference
between the joint entropy and the mutual information (see references).

## References

Cover, T. M. and Thomas, J. A. (2006). *Elements of information theory.*
John Wiley & Sons, 2 edition.

## Examples

``` r
data(fl25)
data(fl25_enum)

plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
distances <- redist.distances(plans_05)
distances$Hamming[1:5, 1:5]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    2    3    2    8
#> [2,]    2    0    4    3    9
#> [3,]    3    4    0    1    5
#> [4,]    2    3    1    0    6
#> [5,]    8    9    5    6    0
```
