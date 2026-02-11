# Redistricting via Compact Random Seed and Grow Algorithm

`redist.crsg` generates redistricting plans using a random seed a grow
algorithm. This is the compact districting algorithm described in Chen
and Rodden (2013).

## Usage

``` r
redist.crsg(
  adj,
  total_pop,
  shp,
  ndists,
  pop_tol,
  verbose = TRUE,
  maxiter = 5000
)
```

## Arguments

- adj:

  List of length N, where N is the number of precincts. Each list
  element is an integer vector indicating which precincts that precinct
  is adjacent to. It is assumed that precinct numbers start at 0.

- total_pop:

  numeric vector of length N, where N is the number of precincts. Each
  element lists the population total of the corresponding precinct, and
  is used to enforce pop_tol constraints.

- shp:

  An sf dataframe to compute area and centroids with.

- ndists:

  integer, the number of districts we want to partition the precincts
  into.

- pop_tol:

  numeric, indicating how close district population targets have to be
  to the target population before algorithm converges. pop_tol=0.05 for
  example means that all districts must be between 0.95 and 1.05 times
  the size of target.pop in population size.

- verbose:

  boolean, indicating whether the time to run the algorithm is printed.

- maxiter:

  integer, indicating maximum number of iterations to attempt before
  convergence to population constraint fails. If it fails once, it will
  use a different set of start values and try again. If it fails again,
  redist.rsg() returns an object of all NAs, indicating that use of more
  iterations may be advised. Default is 5000.

## Value

list, containing three objects containing the completed redistricting
plan.

- `plan`: A vector of length N, indicating the district membership of
  each precinct.

- `district_list` A list of length Ndistrict. Each list contains a
  vector of the precincts in the respective district.

- `district_pop` A vector of length Ndistrict, containing the population
  totals of the respective districts.

## References

Jowei Chen and Jonathan Rodden (2013) “Unintentional Gerrymandering:
Political Geography and Electoral Bias in Legislatures.” Quarterly
Journal of Political Science. 8(3): 239-269.

## Examples

``` r
data("fl25")
adj <- redist.adjacency(fl25)
redist.crsg(adj = adj, total_pop = fl25$pop, shp = fl25, ndists = 2, pop_tol = .1)
#> 
#> ==================== 
#> redist.crsg(): Automated Redistricting Starts
#> 
#> 
#>  2 districts built using 25 precincts in 0 seconds...
#> 
#> $plan
#>  [1] 1 2 2 2 2 1 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2
#> 
#> $district_list
#> $district_list[[1]]
#> [1] 14 13  0  8  5 12 11 10  9
#> 
#> $district_list[[2]]
#>  [1]  7 20 22 18 21 23 24  4 17 19  6 15  2 16  3  1
#> 
#> 
#> $district_pop
#> [1] 81127 93916
#> 
```
