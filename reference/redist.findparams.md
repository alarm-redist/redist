# Run parameter testing for `redist_flip`

`redist.findparams` is used to find optimal parameter values of
`redist_flip` for a given map.

## Usage

``` r
redist.findparams(
  map,
  nsims,
  init_plan = NULL,
  adapt_lambda = FALSE,
  adapt_eprob = FALSE,
  params,
  ssdmat = NULL,
  group_pop = NULL,
  counties = NULL,
  nstartval_store = 1,
  maxdist_startval = 100,
  maxiterrsg = 5000,
  report_all = TRUE,
  parallel = FALSE,
  ncores = NULL,
  log = FALSE,
  verbose = TRUE
)
```

## Arguments

- map:

  A
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object.

- nsims:

  The number of simulations run before a save point.

- init_plan:

  A vector containing the congressional district labels of each
  geographic unit. The default is `NULL`. If not provided, random and
  contiguous congressional district assignments will be generated using
  `redist.rsg`.

- adapt_lambda:

  Whether to adaptively tune the lambda parameter so that the
  Metropolis-Hastings acceptance probability falls between 20% and 40%.
  Default is FALSE.

- adapt_eprob:

  Whether to adaptively tune the edgecut probability parameter so that
  the Metropolis-Hastings acceptance probability falls between 20% and
  40%. Default is FALSE.

- params:

  A matrix of parameter values to test, such as the output of
  `expand.grid`. Parameters accepted for `params` include `eprob`,
  `lambda`, `pop_tol`, `beta`, and `constraint`.

- ssdmat:

  A matrix of squared distances between geographic units. The default is
  `NULL`.

- group_pop:

  A vector of populations for some sub-group of interest. The default is
  `NULL`.

- counties:

  A vector of county membership assignments. The default is `NULL`.

- nstartval_store:

  The number of maps to sample from the preprocessing chain for use as
  starting values in future simulations. Default is 1.

- maxdist_startval:

  The maximum distance from the starting map that sampled maps should
  be. Default is 100 (no restriction).

- maxiterrsg:

  Maximum number of iterations for random seed-and-grow algorithm to
  generate starting values. Default is 5000.

- report_all:

  Whether to report all summary statistics for each set of parameter
  values. Default is `TRUE`.

- parallel:

  Whether to run separate parameter settings in parallel. Default is
  `FALSE`.

- ncores:

  Number of parallel tasks to run, declared outside of the function.
  Default is `NULL`.

- log:

  Whether to open a log to track progress for each parameter combination
  being tested. Default is FALSE.

- verbose:

  Whether to print additional information about the tests. Default is
  `TRUE`.

## Value

`redist.findparams` returns a print-out of summary statistics about each
parameter setting.

## Details

This function allows users to test multiple parameter settings of
`redist_flip` in preparation for a longer run for analysis.

## References

Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
(2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
Carlo." Working Paper. Available at
<http://imai.princeton.edu/research/files/redist.pdf>.

## Examples

``` r
# \donttest{
data(fl25)
data(fl25_enum)
data(fl25_adj)

## Get an initial partition
init_plan <- fl25_enum$plans[, 5118]

params <- expand.grid(eprob = c(.01, .05, .1))

# Make map
map_fl <- redist_map(fl25, ndists = 3, pop_tol = 0.2)
#> Projecting to CRS 3857
## Run the algorithm
redist.findparams(map_fl,
    init_plan = init_plan, nsims = 10000, params = params)
#> ## ------------------------------
#>  ## redist.findparams(): Parameter tuning for redist_flip()
#>  ## Searching over 3 parameter combinations
#>  ## ------------------------------
#> 
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> Starting chain 1
#> ℹ Starting swMH().
#> ■                                  0% | ETA: 8s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.99
#> 
#> No maps available under parameter set 1.
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> Starting chain 1
#> ℹ Starting swMH().
#> ■■■■■■■■■■■■■                     40% | ETA:  1s | MH Acceptance: 0.96
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.96
#> 
#> No maps available under parameter set 2.
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> Starting chain 1
#> ℹ Starting swMH().
#> ■■■■■■■■■■■■■■■■■■■■■■■■■         81% | ETA:  0s | MH Acceptance: 0.93
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.93
#> 
#> No maps available under parameter set 3.
#>  ## -------------------------------------
#>  ## -------------------------------------
#>  ## Parameter Values for Simulation 1 
#>  ## No valid maps generated under these parameters
#>  ## -------------------------------------
#>  ## -------------------------------------
#> 
#>  ## -------------------------------------
#>  ## -------------------------------------
#>  ## Parameter Values for Simulation 2 
#>  ## No valid maps generated under these parameters
#>  ## -------------------------------------
#>  ## -------------------------------------
#> 
#>  ## -------------------------------------
#>  ## -------------------------------------
#>  ## Parameter Values for Simulation 3 
#>  ## No valid maps generated under these parameters
#>  ## -------------------------------------
#>  ## -------------------------------------
#> 
#> $diagnostics
#> [1] " ## -------------------------------------\n ## -------------------------------------\n ## Parameter Values for Simulation 1 \n ## No valid maps generated under these parameters\n ## -------------------------------------\n ## -------------------------------------\n\n ## -------------------------------------\n ## -------------------------------------\n ## Parameter Values for Simulation 2 \n ## No valid maps generated under these parameters\n ## -------------------------------------\n ## -------------------------------------\n\n ## -------------------------------------\n ## -------------------------------------\n ## Parameter Values for Simulation 3 \n ## No valid maps generated under these parameters\n ## -------------------------------------\n ## -------------------------------------\n\n"
#> 
#> $startvals
#> $startvals[[1]]
#> NULL
#> 
#> 
# }
```
