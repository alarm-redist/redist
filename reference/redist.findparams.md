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
#> ℹ Starting swMH().
#> ■                                  0% | ETA: 8s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.99
#> 
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> ℹ Starting swMH().
#> ■■                                 4% | ETA: 2s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.96
#> 
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> ℹ Starting swMH().
#> ■■■■■■■■■■■■■■■                   48% | ETA:  1s | MH Acceptance: 0.93
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.92
#> 
#>  ## -------------------------------------
#>  ## -------------------------------------
#>  ## Parameter Values for Simulation 1 
#> ## Edgecut probability = 0.01
#> ## Lambda = 0
#> ## No hard population constraint applied
#> ## -------------------------------------
#> ## Diagnostics:
#> ## Metropolis-Hastings Acceptance Ratio = 0.993
#> ## Mean population parity distance = 0.851
#> ## Median population parity distance = 0.828
#> ## Population parity range = 0.029 1.93
#> ## MCMC Iteration quantiles of population parity median = 0.876 0.807 0.876 0.768
#> 
#> ## Mean share of geographies equal to initial assignment = 0.352
#> ## Median share of geographies equal to initial assignment = 0.36
#> ## Range of share of geographies equal to initial assignment = 0 0.96
#> ## MCMC Iteration quantiles of geography distance to initial assignment = 0.4 0.4 0.32 0.24
#> ## -------------------------------------
#> ## -------------------------------------
#> 
#>  ## -------------------------------------
#>  ## -------------------------------------
#>  ## Parameter Values for Simulation 1 
#> ## Edgecut probability = 0.05
#> ## Lambda = 0
#> ## No hard population constraint applied
#> ## -------------------------------------
#> ## Diagnostics:
#> ## Metropolis-Hastings Acceptance Ratio = 0.96
#> ## Mean population parity distance = 0.837
#> ## Median population parity distance = 0.835
#> ## Population parity range = 0.008 1.868
#> ## MCMC Iteration quantiles of population parity median = 0.717 0.963 0.808 0.586
#> 
#> ## Mean share of geographies equal to initial assignment = 0.325
#> ## Median share of geographies equal to initial assignment = 0.32
#> ## Range of share of geographies equal to initial assignment = 0 0.96
#> ## MCMC Iteration quantiles of geography distance to initial assignment = 0.36 0.36 0.28 0.28
#> ## -------------------------------------
#> ## -------------------------------------
#> 
#>  ## -------------------------------------
#>  ## -------------------------------------
#>  ## Parameter Values for Simulation 1 
#> ## Edgecut probability = 0.1
#> ## Lambda = 0
#> ## No hard population constraint applied
#> ## -------------------------------------
#> ## Diagnostics:
#> ## Metropolis-Hastings Acceptance Ratio = 0.923
#> ## Mean population parity distance = 0.899
#> ## Median population parity distance = 0.893
#> ## Population parity range = 0.009 1.946
#> ## MCMC Iteration quantiles of population parity median = 0.875 0.891 0.912 0.889
#> 
#> ## Mean share of geographies equal to initial assignment = 0.313
#> ## Median share of geographies equal to initial assignment = 0.28
#> ## Range of share of geographies equal to initial assignment = 0 0.96
#> ## MCMC Iteration quantiles of geography distance to initial assignment = 0.28 0.24 0.32 0.36
#> ## -------------------------------------
#> ## -------------------------------------
#> 
#> $diagnostics
#> [1] " ## -------------------------------------\n ## -------------------------------------\n ## Parameter Values for Simulation 1 \n## Edgecut probability = 0.01\n## Lambda = 0\n## No hard population constraint applied\n## -------------------------------------\n## Diagnostics:\n## Metropolis-Hastings Acceptance Ratio = 0.993\n## Mean population parity distance = 0.851\n## Median population parity distance = 0.828\n## Population parity range = 0.029 1.93\n## MCMC Iteration quantiles of population parity median = 0.876 0.807 0.876 0.768\n\n## Mean share of geographies equal to initial assignment = 0.352\n## Median share of geographies equal to initial assignment = 0.36\n## Range of share of geographies equal to initial assignment = 0 0.96\n## MCMC Iteration quantiles of geography distance to initial assignment = 0.4 0.4 0.32 0.24\n## -------------------------------------\n## -------------------------------------\n\n ## -------------------------------------\n ## -------------------------------------\n ## Parameter Values for Simulation 1 \n## Edgecut probability = 0.05\n## Lambda = 0\n## No hard population constraint applied\n## -------------------------------------\n## Diagnostics:\n## Metropolis-Hastings Acceptance Ratio = 0.96\n## Mean population parity distance = 0.837\n## Median population parity distance = 0.835\n## Population parity range = 0.008 1.868\n## MCMC Iteration quantiles of population parity median = 0.717 0.963 0.808 0.586\n\n## Mean share of geographies equal to initial assignment = 0.325\n## Median share of geographies equal to initial assignment = 0.32\n## Range of share of geographies equal to initial assignment = 0 0.96\n## MCMC Iteration quantiles of geography distance to initial assignment = 0.36 0.36 0.28 0.28\n## -------------------------------------\n## -------------------------------------\n\n ## -------------------------------------\n ## -------------------------------------\n ## Parameter Values for Simulation 1 \n## Edgecut probability = 0.1\n## Lambda = 0\n## No hard population constraint applied\n## -------------------------------------\n## Diagnostics:\n## Metropolis-Hastings Acceptance Ratio = 0.923\n## Mean population parity distance = 0.899\n## Median population parity distance = 0.893\n## Population parity range = 0.009 1.946\n## MCMC Iteration quantiles of population parity median = 0.875 0.891 0.912 0.889\n\n## Mean share of geographies equal to initial assignment = 0.313\n## Median share of geographies equal to initial assignment = 0.28\n## Range of share of geographies equal to initial assignment = 0 0.96\n## MCMC Iteration quantiles of geography distance to initial assignment = 0.28 0.24 0.32 0.36\n## -------------------------------------\n## -------------------------------------\n\n"
#> 
#> $startvals
#> $startvals[[1]]
#>       [,1]
#>  [1,]    3
#>  [2,]    1
#>  [3,]    2
#>  [4,]    1
#>  [5,]    1
#>  [6,]    2
#>  [7,]    1
#>  [8,]    1
#>  [9,]    2
#> [10,]    3
#> [11,]    2
#> [12,]    2
#> [13,]    2
#> [14,]    2
#> [15,]    2
#> [16,]    2
#> [17,]    2
#> [18,]    2
#> [19,]    1
#> [20,]    1
#> [21,]    2
#> [22,]    1
#> [23,]    2
#> [24,]    1
#> [25,]    1
#> 
#> $startvals[[2]]
#>       [,1]
#>  [1,]    2
#>  [2,]    2
#>  [3,]    2
#>  [4,]    2
#>  [5,]    2
#>  [6,]    3
#>  [7,]    2
#>  [8,]    2
#>  [9,]    1
#> [10,]    2
#> [11,]    2
#> [12,]    2
#> [13,]    2
#> [14,]    2
#> [15,]    2
#> [16,]    2
#> [17,]    2
#> [18,]    2
#> [19,]    2
#> [20,]    2
#> [21,]    2
#> [22,]    2
#> [23,]    3
#> [24,]    2
#> [25,]    2
#> 
#> $startvals[[3]]
#>       [,1]
#>  [1,]    1
#>  [2,]    3
#>  [3,]    3
#>  [4,]    3
#>  [5,]    1
#>  [6,]    1
#>  [7,]    1
#>  [8,]    1
#>  [9,]    1
#> [10,]    2
#> [11,]    2
#> [12,]    2
#> [13,]    2
#> [14,]    1
#> [15,]    2
#> [16,]    3
#> [17,]    3
#> [18,]    1
#> [19,]    3
#> [20,]    3
#> [21,]    1
#> [22,]    3
#> [23,]    1
#> [24,]    3
#> [25,]    3
#> 
#> 
# }
```
