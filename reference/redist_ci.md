# Confidence Intervals for SMC and MCMC Estimates

Builds a confidence interval for a quantity of interest. If multiple
runs are available, uses the between-run variation to estimate the
standard error. If only one run is available, uses information on the
SMC particle/plan genealogy to estimate the standard error, using a
variant of the method of Olson & Douc (2019). The multiple-run estimator
is more reliable, especially for situations with many districts, and
should be used when parallelism is available. All reference plans are
ignored.

## Usage

``` r
redist_ci(plans, x, district = 1L, conf = 0.9, by_chain = FALSE)

redist_smc_ci(plans, x, district = 1L, conf = 0.9, by_chain = FALSE)

redist_mcmc_ci(plans, x, district = 1L, conf = 0.9, by_chain = FALSE)
```

## Arguments

- plans:

  a
  [redist_plans](http://alarm-redist.org/redist/reference/redist_plans.md)
  object.

- x:

  the quantity to build an interval for. Tidy-evaluated within `plans`.

- district:

  for
  [redist_plans](http://alarm-redist.org/redist/reference/redist_plans.md)
  objects with multiple districts, which `district` to subset to. Set to
  `NULL` to perform no subsetting.

- conf:

  the desired confidence level.

- by_chain:

  Whether the confidence interval should indicate overall sampling
  uncertainty (`FALSE`) or per-chain sampling uncertainty (`TRUE`). In
  the latter case the intervals will be wider by a factor of
  `sqrt(runs)`.

## Value

A tibble with three columns: `X`, `X_lower`, and `X_upper`, where `X` is
the name of the vector of interest, containing the mean and confidence
interval. When used inside
[`summarize()`](https://dplyr.tidyverse.org/reference/summarise.html)
this will create three columns in the output data.

## Functions

- `redist_smc_ci()`: Compute confidence intervals for SMC output.

- `redist_mcmc_ci()`: Compute confidence intervals for MCMC output.

## References

Lee, A., & Whiteley, N. (2018). Variance estimation in the particle
filter. Biometrika, 105(3), 609-625.

Olsson, J., & Douc, R. (2019). Numerically stable online estimation of
variance in particle filters. Bernoulli, 25(2), 1504-1535.

H. P. Chan and T. L. Lai. A general theory of particle filters in hidden
Markov models and some applications. Ann. Statist., 41(6):2877–2904,
2013.

## Examples

``` r
library(dplyr)
data(iowa)

iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
plans <- redist_mergesplit_parallel(iowa_map, nsims = 200, chains = 2, silent = TRUE) %>%
    mutate(dem = group_frac(iowa_map, dem_08, dem_08 + rep_08)) %>%
    number_by(dem)
redist_smc_ci(plans, dem)
#> Warning: Runs have not converged for this statistic.
#> ℹ R-hat is 2.292
#> → Increase the number of samples.
#> # A tibble: 1 × 3
#>     dem dem_lower dem_upper
#>   <dbl>     <dbl>     <dbl>
#> 1 0.478     0.467     0.490
```
