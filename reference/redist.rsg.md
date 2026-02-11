# Redistricting via Random Seed and Grow Algorithm

`redist.rsg` generates redistricting plans using a random seed a grow
algorithm. This is the non-compact districting algorithm described in
Chen and Rodden (2013). The algorithm can provide start values for the
other redistricting routines in this package.

## Usage

``` r
redist.rsg(adj, total_pop, ndists, pop_tol, verbose = TRUE, maxiter = 5000)
```

## Arguments

- adj:

  List of length N, where N is the number of precincts. Each list
  element is an integer vector indicating which precincts that precinct
  is adjacent to. It is assumed that precinct numbers start at 0.

- total_pop:

  numeric vector of length N, where N is the number of precincts. Each
  element lists the population total of the corresponding precinct, and
  is used to enforce population constraints.

- ndists:

  integer, the number of districts we want to partition the precincts
  into.

- pop_tol:

  numeric, indicating how close district population targets have to be
  to the target population before algorithm converges. thresh=0.05 for
  example means that all districts must be between 0.95 and 1.05 times
  the size of target.pop in population size.

- verbose:

  boolean, indicating whether the time to run the algorithm is printed.

- maxiter:

  integer, indicating maximum number of iterations to attempt before
  convergence to population constraint fails. If it fails once, it will
  use a different set of start values and try again. If it fails again,
  redist.rsg() returns an object of all NAs, indicating that use of more
  iterations may be advised.

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

## Author

Benjamin Fifield, Department of Politics, Princeton University
<benfifield@gmail.com>, <https://www.benfifield.com/>

Michael Higgins, Department of Statistics, Kansas State University
<mikehiggins@k-state.edu>,
<https://www.k-state.edu/stats/about/people/HigginsMichael.html>

Kosuke Imai, Department of Politics, Princeton University
<imai@harvard.edu>, <https://imai.fas.harvard.edu>

James Lo, <jameslo@princeton.edu>

Alexander Tarr, Department of Electrical Engineering, Princeton
University <atarr@princeton.edu>

## Examples

``` r
### Real data example from test set
data(fl25)
data(fl25_adj)

res <- redist.rsg(adj = fl25_adj, total_pop = fl25$pop,
    ndists = 3, pop_tol = 0.05)
#> 
#> ==================== 
#> redist.rsg(): Automated Redistricting Starts
#> 
#> 
#>  3 districts built using 25 precincts in 0 seconds...
#> 
```
