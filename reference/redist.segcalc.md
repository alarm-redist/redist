# Segregation index calculation for MCMC redistricting.

`redist.segcalc` calculates the dissimilarity index of segregation (see
Massey & Denton 1987 for more details) for a specified subgroup under
any redistricting plan.

## Usage

``` r
segregation_index(
  map,
  group_pop,
  total_pop = map[[attr(map, "pop_col")]],
  .data = cur_plans()
)

redist.segcalc(plans, group_pop, total_pop)
```

## Arguments

- map:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- group_pop:

  A vector of populations for some subgroup of interest.

- total_pop:

  A vector containing the populations of each geographic unit.

- .data:

  a
  [`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
  object

- plans:

  A matrix of congressional district assignments or a redist object.

## Value

`redist.segcalc` returns a vector where each entry is the dissimilarity
index of segregation (Massey & Denton 1987) for each redistricting plan
in `algout`.

## References

Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
(2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
Carlo." Working Paper. Available at
<http://imai.princeton.edu/research/files/redist.pdf>.

Massey, Douglas and Nancy Denton. (1987) "The Dimensions of Social
Segregation". Social Forces.

## Examples

``` r
# \donttest{
data(fl25)
data(fl25_enum)
data(fl25_adj)

## Get an initial partition
init_plan <- fl25_enum$plans[, 5118]
fl25$init_plan <- init_plan

## 25 precinct, three districts - no pop constraint ##
fl_map <- redist_map(fl25, existing_plan = 'init_plan', adj = fl25_adj)
#> Projecting to CRS 3857
alg_253 <- redist_flip(fl_map, nsims = 10000)
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> Starting chain 1
#> ℹ Starting swMH().
#> ■                                  0% | ETA: 9s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.78
#> 


## Get Republican Dissimilarity Index from simulations
# old: rep_dmi_253 <- redist.segcalc(alg_253, fl25$mccain, fl25$pop)
rep_dmi_253 <- seg_dissim(alg_253, fl25, mccain, pop)  |>
    redistmetrics::by_plan(ndists = 3)
# }
```
