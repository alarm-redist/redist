# A set of redistricting plans

A `redist_plans` object is essentially a data frame of summary
information on each district and each plan, along with the matrix of
district assignments and information about the simulation process used
to generate the plans.

## Usage

``` r
redist_plans(plans, map, algorithm, wgt = NULL, ...)
```

## Arguments

- plans:

  a matrix with `n_precinct` columns and `n_sims` rows, or a single
  vector of precinct assignments.

- map:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- algorithm:

  the algorithm used to generate the plans (usually "smc" or "mcmc")

- wgt:

  the weights to use, if any.

- ...:

  Other named attributes to set

## Value

a new `redist_plans` object.

## Details

The first two columns of the data frame will be `draw`, a factor
indexing the simulation draw, and `district`, an integer indexing the
districts within a plan. The data frame will therefore have
`n_sims*ndists` rows. As a data frame, the usual `dplyr` methods will
work.

Other useful methods for `redist_plans` objects:

- [`summary.redist_plans`](http://alarm-redist.org/redist/reference/summary.redist_plans.md)

- [`add_reference`](http://alarm-redist.org/redist/reference/add_reference.md)

- [`subset_sampled`](http://alarm-redist.org/redist/reference/subset_sampled.md)

- [`subset_ref`](http://alarm-redist.org/redist/reference/subset_sampled.md)

- [`pullback`](http://alarm-redist.org/redist/reference/pullback.md)

- [`number_by`](http://alarm-redist.org/redist/reference/number_by.md)

- [`match_numbers`](http://alarm-redist.org/redist/reference/match_numbers.md)

- [`is_county_split`](http://alarm-redist.org/redist/reference/is_county_split.md)

- [`prec_assignment`](http://alarm-redist.org/redist/reference/prec_assignment.md)

- [`plan_distances`](http://alarm-redist.org/redist/reference/redist.distances.md)

- [`get_plans_matrix`](http://alarm-redist.org/redist/reference/get_plans_matrix.md)

- [`get_plans_weights`](http://alarm-redist.org/redist/reference/get_plans_weights.md)

- [`get_sampling_info`](http://alarm-redist.org/redist/reference/get_sampling_info.md)

- [`as.matrix.redist_plans`](http://alarm-redist.org/redist/reference/get_plans_matrix.md)

- [`plot.redist_plans`](http://alarm-redist.org/redist/reference/plot.redist_plans.md)

## Examples

``` r
data(iowa)

iowa <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)
rsg_plan <- redist.rsg(iowa$adj, iowa$pop, ndists = 4, pop_tol = 0.05)$plan
#> 
#> ==================== 
#> redist.rsg(): Automated Redistricting Starts
#> 
#> 
#>  4 districts built using 99 precincts in 0.05 seconds...
#> 
redist_plans(rsg_plan, iowa, "rsg")
#> A <redist_plans> containing 1 sampled plan
#> Plans have 4 districts from a 99-unit map, and were drawn using random
#> seed-and-grow.
#> Plans matrix: int [1:99, 1] 1 1 2 3 4 2 2 1 2 2 ...
#> # A tibble: 4 × 3
#>   draw  district total_pop
#> * <fct>    <int>     <dbl>
#> 1 1            1    735276
#> 2 1            2    777906
#> 3 1            3    763285
#> 4 1            4    769888
```
