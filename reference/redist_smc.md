# SMC Redistricting Sampler (McCartan and Imai 2023)

`redist_smc` uses a Sequential Monte Carlo algorithm (McCartan and Imai
2023) to generate representative samples of congressional or legislative
redistricting plans according to contiguity, population, compactness,
and administrative boundary constraints.

## Usage

``` r
redist_smc(
  map,
  nsims,
  counties = NULL,
  compactness = 1,
  constraints = list(),
  resample = TRUE,
  runs = 1L,
  ncores = 0L,
  init_particles = NULL,
  n_steps = NULL,
  adapt_k_thresh = 0.99,
  seq_alpha = 0.5,
  truncate = (compactness != 1),
  trunc_fn = redist_quantile_trunc,
  pop_temper = 0,
  final_infl = 1,
  ref_name = NULL,
  verbose = FALSE,
  silent = FALSE
)
```

## Arguments

- map:

  A
  [`redist_map()`](http://alarm-redist.org/redist/reference/redist_map.md)
  object.

- nsims:

  The number of samples to draw.

- counties:

  A vector containing county (or other administrative or geographic
  unit) labels for each unit, which may be integers ranging from 1 to
  the number of counties, or a factor or character vector. If provided,
  the algorithm will only generate maps which split up to `ndists-1`
  counties. Even there are fewer counties than `ndists - 1`, the
  spanning trees will change the results of the simulations. There is no
  strength parameter associated with this constraint. To adjust the
  number of county splits further, or to constrain a second type of
  administrative split, consider using
  [`add_constr_splits()`](http://alarm-redist.org/redist/reference/constraints.md),
  [`add_constr_multisplits()`](http://alarm-redist.org/redist/reference/constraints.md),
  and
  [`add_constr_total_splits()`](http://alarm-redist.org/redist/reference/constraints.md).

- compactness:

  Controls the compactness of the generated districts, with higher
  values preferring more compact districts. Must be nonnegative. See the
  'Details' section for more information, and computational
  considerations.

- constraints:

  A
  [`redist_constr()`](http://alarm-redist.org/redist/reference/redist_constr.md)
  object or a list containing information on sampling constraints. See
  [constraints](http://alarm-redist.org/redist/reference/constraints.md)
  for more information.

- resample:

  Whether to perform a final resampling step so that the generated plans
  can be used immediately. Set this to `FALSE` to perform direct
  importance sampling estimates, or to adjust the weights manually.

- runs:

  How many independent parallel runs to conduct. Each run will have
  `nsims` simulations. Multiple runs allows for estimation of simulation
  standard errors. Output will only be shown for the first run. For
  compatibility with MCMC methods, runs are identified with the `chain`
  column in the output.

- ncores:

  How many cores to use to parallelize plan generation within each run.
  The default, 0, will use the number of available cores on the machine
  as long as `nsims` and the number of units is large enough. If
  `runs>1` you will need to set this manually. If more than one core is
  used, the sampler output will not be fully reproducible with
  [`set.seed()`](https://rdrr.io/r/base/Random.html). If full
  reproducibility is desired, set `ncores=1`.

- init_particles:

  A matrix of partial plans to begin sampling from. For advanced use
  only. The matrix must have `nsims` columns and a row for every
  precinct. It is important to ensure that the existing districts meet
  contiguity and population constraints, or there may be major issues
  when sampling.

- n_steps:

  How many steps to run the SMC algorithm for. Each step splits off a
  new district. Defaults to all remaining districts. If fewer than the
  number of remaining splits, reference plans are disabled.

- adapt_k_thresh:

  The threshold value used in the heuristic to select a value `k_i` for
  each splitting iteration. Higher values are more accurate but may
  require more computation. Set to 1 for the most conservative sampling.
  Must be between 0 and 1.

- seq_alpha:

  The amount to adjust the weights by at each resampling step; higher
  values prefer exploitation, while lower values prefer exploration.
  Must be between 0 and 1.

- truncate:

  Whether to truncate the importance sampling weights at the final step
  by `trunc_fn`. Recommended if `compactness` is not 1. Truncation only
  applied if `resample=TRUE`.

- trunc_fn:

  A function which takes in a vector of weights and returns a truncated
  vector. If the [loo](https://mc-stan.org/loo/reference/loo.html)
  package is installed (strongly recommended), will default to
  Pareto-smoothed Importance Sampling (PSIS) rather than naive
  truncation.

- pop_temper:

  The strength of the automatic population tempering. Try values of
  0.01-0.05 to start if the algorithm gets stuck on the final few
  splits.

- final_infl:

  A multiplier for the population constraint on the final iteration.
  Used to loosen the constraint when the sampler is getting stuck on the
  final split. `pop_temper` should be tried first, since using
  `final_infl` will actually change the target distribution.

- ref_name:

  a name for the existing plan, which will be added as a reference plan,
  or `FALSE` to not include the initial plan in the output. Defaults to
  the column name of the existing plan.

- verbose:

  Whether to print out intermediate information while sampling.
  Recommended.

- silent:

  Whether to suppress all diagnostic information.

## Value

`redist_smc` returns a
[redist_plans](http://alarm-redist.org/redist/reference/redist_plans.md)
object containing the simulated plans.

## Details

This function draws samples from a specific target measure controlled by
the `map`, `compactness`, and `constraints` parameters.

Key to ensuring good performance is monitoring the efficiency of the
resampling process at each SMC stage. Unless `silent=FALSE`, this
function will print out the effective sample size of each resampling
step to allow the user to monitor the efficiency. If `verbose=TRUE` the
function will also print out information on the \\k_i\\ values
automatically chosen and the acceptance rate (based on the population
constraint) at each step. Users should also check diagnostics of the
sample by running
[`summary.redist_plans()`](http://alarm-redist.org/redist/reference/summary.redist_plans.md).

Higher values of `compactness` sample more compact districts; setting
this parameter to 1 is computationally efficient and generates nicely
compact districts. Values of other than 1 may lead to highly variable
importance sampling weights. In these cases, these weights are by
default truncated using
[`redist_quantile_trunc()`](http://alarm-redist.org/redist/reference/redist_quantile_trunc.md)
to stabilize the resulting estimates, but if truncation is used, a
specific truncation function should probably be chosen by the user.

## References

McCartan, C., & Imai, K. (2023). Sequential Monte Carlo for Sampling
Balanced and Compact Redistricting Plans. *Annals of Applied Statistics*
17(4). Available at
[doi:10.1214/23-AOAS1763](https://doi.org/10.1214/23-AOAS1763) .

## Examples

``` r
# \donttest{
data(fl25)

fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#> Projecting to CRS 3857

sampled_basic <- redist_smc(fl_map, 5000)
#> SEQUENTIAL MONTE CARLO
#> Sampling 5000 25-unit maps with 3 districts and population between 52513 and 64182.

constr <- redist_constr(fl_map)
constr <- add_constr_incumbency(constr, strength = 100, incumbents = c(3, 6, 25))
sampled_constr <- redist_smc(fl_map, 5000, constraints = constr)
#> SEQUENTIAL MONTE CARLO
#> Sampling 5000 25-unit maps with 3 districts and population between 52513 and 64182.

# Multiple parallel independent runs
redist_smc(fl_map, 1000, runs = 2)
#> A <redist_plans> containing 2,000 sampled plans
#> Plans have 3 districts from a 25-unit map, and were drawn using Sequential
#> Monte Carlo.
#> With plans resampled from weights
#> Plans matrix: int [1:25, 1:2000] 2 1 1 1 2 1 2 3 2 3 ...
#> # A tibble: 6,000 × 4
#>    draw  chain district total_pop
#>    <fct> <int>    <int>     <dbl>
#>  1 1         1        1     52653
#>  2 1         1        2     60679
#>  3 1         1        3     61711
#>  4 2         1        1     62770
#>  5 2         1        2     56227
#>  6 2         1        3     56046
#>  7 3         1        1     57349
#>  8 3         1        2     55024
#>  9 3         1        3     62670
#> 10 4         1        1     54084
#> # ℹ 5,990 more rows

# One run with multiple cores
redist_smc(fl_map, 1000, ncores = 2)
#> SEQUENTIAL MONTE CARLO
#> Sampling 1000 25-unit maps with 3 districts and population between 52513 and 64182.
#> Split [0/2] ■                                | ETA?
#> Split [2/2] ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | ETA 0s
#> 
#> A <redist_plans> containing 1,000 sampled plans
#> Plans have 3 districts from a 25-unit map, and were drawn using Sequential
#> Monte Carlo.
#> With plans resampled from weights
#> Plans matrix: int [1:25, 1:1000] 3 2 2 2 2 3 1 1 3 1 ...
#> # A tibble: 3,000 × 3
#>    draw  district total_pop
#>  * <fct>    <int>     <dbl>
#>  1 1            1     55024
#>  2 1            2     58845
#>  3 1            3     61174
#>  4 2            1     57208
#>  5 2            2     54128
#>  6 2            3     63707
#>  7 3            1     57892
#>  8 3            2     55321
#>  9 3            3     61830
#> 10 4            1     56832
#> # ℹ 2,990 more rows
# }
```
