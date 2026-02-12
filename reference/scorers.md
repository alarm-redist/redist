# Scoring functions for `redist_shortburst`

The output of these functions may be passed into
[`redist_shortburst()`](http://alarm-redist.org/redist/reference/redist_shortburst.md)
as `score_fn`. Scoring functions have type `redist_scorer` and may be
combined together using basic arithmetic operations.

## Usage

``` r
scorer_group_pct(map, group_pop, total_pop, k = 1)

scorer_pop_dev(map)

scorer_splits(map, counties)

scorer_multisplits(map, counties)

scorer_frac_kept(map)

scorer_polsby_popper(map, perim_df = NULL, areas = NULL, m = 1)

scorer_status_quo(map, existing_plan = get_existing(map))
```

## Arguments

- map:

  A
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object.

- group_pop:

  A numeric vector with the population of the group for every precinct.

- total_pop:

  A numeric vector with the population for every precinct.

- k:

  the k-th from the top group fraction to return as the score.

- counties:

  A numeric vector with an integer from 1:n_counties

- perim_df:

  perimeter distance dataframe from
  [`redistmetrics::prep_perims()`](http://alarm-redist.org/redistmetrics/reference/prep_perims.md)

- areas:

  area of each precinct (ie `st_area(map)`)

- m:

  the m-th from the bottom Polsby Popper to return as the score.
  Defaults to 1, the minimum Polsby Popper score

- existing_plan:

  A vector containing the current plan.

## Value

A scoring function of class `redist_scorer` which returns a single
numeric value per plan. Larger values are generally better for
`frac_kept`, `group_pct`, and `polsby_popper` and smaller values are
better for `splits` and `pop_dev`.

## Details

Function details:

- `scorer_group_pct` returns the `k`-th top group percentage across
  districts. For example, if the group is Democratic voters and `k=3`,
  then the function returns the 3rd-highest fraction of Democratic
  voters across all districts. Can be used to target `k` VRA districts
  or partisan gerrymanders.

- `scorer_pop_dev` returns the maximum population deviation within a
  plan. Smaller values are closer to population parity, so use
  `maximize=FALSE` with this scorer.

- `scorer_splits` returns the fraction of counties that are split within
  a plan. Higher values have more county splits, so use `maximize=FALSE`
  with this scorer.

- `scorer_frac_kept` returns the fraction of edges kept in each
  district. Higher values mean more compactness.

- `scorer_polsby_popper` returns the `m`-th Polsby Popper score within a
  plan. Higher scores correspond to more compact districts. Use
  `m=ndists/2` to target the median compactness, `m=1` to target the
  minimum compactness.

- `scorer_status_quo` returns 1 - the rescaled variation of information
  distance between the plan and the `existing_plan`. Larger values
  indicate the plan is closer to the existing plan.

## Examples

``` r
# \donttest{
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)

scorer_frac_kept(iowa_map)
#> function (plans) 
#> {
#>     (edges - n_removed(adj, plans, ndists))/edges
#> }
#> <bytecode: 0x55a8a0893230>
#> <environment: 0x55a8b0208e48>
#> attr(,"class")
#> [1] "redist_scorer" "function"     
scorer_status_quo(iowa_map)
#> function (plans) 
#> {
#>     1 - 0.5 * var_info_vec(plans, existing_plan, pop)/log(ndists)
#> }
#> <bytecode: 0x55a8b0c0ce98>
#> <environment: 0x55a8b0c0d6b0>
#> attr(,"class")
#> [1] "redist_scorer" "function"     
scorer_group_pct(iowa_map, dem_08, tot_08, k = 2)
#> function (plans) 
#> {
#>     group_pct_top_k(plans, group_pop, total_pop, k, ndists)
#> }
#> <bytecode: 0x55a8aa28d250>
#> <environment: 0x55a8aa28db10>
#> attr(,"class")
#> [1] "redist_scorer" "function"     
1.5*scorer_frac_kept(iowa_map) + 0.4*scorer_status_quo(iowa_map)
#> function (plans) 
#> {
#>     fn1(plans) + fn2(plans)
#> }
#> <bytecode: 0x55a8aec2fd88>
#> <environment: 0x55a8aec305a0>
#> attr(,"class")
#> [1] "redist_scorer" "function"     
1.5*scorer_frac_kept(iowa_map) + scorer_frac_kept(iowa_map)*scorer_status_quo(iowa_map)
#> function (plans) 
#> {
#>     fn1(plans) + fn2(plans)
#> }
#> <bytecode: 0x55a8aec2fd88>
#> <environment: 0x55a8ab1c9740>
#> attr(,"class")
#> [1] "redist_scorer" "function"     
cbind(
    comp = scorer_frac_kept(iowa_map),
    sq = scorer_status_quo(iowa_map)
)
#> function (plans) 
#> {
#>     do.call(cbind, c(lapply(fns, function(fn) {
#>         fn(plans)
#>     }), list(deparse.level = deparse.level)))
#> }
#> <bytecode: 0x55a8a88ea520>
#> <environment: 0x55a8a953ffc0>
#> attr(,"class")
#> [1] "redist_scorer" "function"     
# }
```
