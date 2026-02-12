# Diagnostic information on sampled plans

Prints diagnostic information, which varies by algorithm. All algorithms
compute the
[`plans_diversity()`](http://alarm-redist.org/redist/reference/plans_diversity.md)
of the samples.

## Usage

``` r
# S3 method for class 'redist_plans'
summary(object, district = 1L, all_runs = TRUE, vi_max = 100, ...)
```

## Arguments

- object:

  a
  [redist_plans](http://alarm-redist.org/redist/reference/redist_plans.md)
  object

- district:

  For R-hat values, which district to use for district-level summary
  statistics. We strongly recommend calling
  [`match_numbers()`](http://alarm-redist.org/redist/reference/match_numbers.md)
  or
  [`number_by()`](http://alarm-redist.org/redist/reference/number_by.md)
  before examining these district-level statistics.

- all_runs:

  When there are multiple SMC runs, show detailed summary statistics for
  all runs (the default), or only the first run?

- vi_max:

  The maximum number of plans to sample in computing the pairwise
  variation of information distance (sample diversity).

- ...:

  additional arguments (ignored)

## Value

A data frame containing diagnostic information, invisibly.

## Details

For SMC and MCMC, if there are multiple runs/chains, R-hat values will
be computed for each summary statistic. These values should be close
to 1. If they are not, then there is too much between-chain variation,
indicating that there are not enough samples. R-hat values are
calculated after rank-normalization and folding. MCMC chains are split
in half before R-hat is computed. For summary statistics that vary
across districts, R-hat is calculated for the first district only.

For SMC, diagnostics statistics include:

- **Effective samples**: the effective sample size at each iteration,
  computed using the SMC weights. Larger is better. The percentage in
  parentheses is the ratio of the effective samples to the total
  samples.

- **Acceptance rate**: the fraction of drawn spanning trees which yield
  a valid redistricting plan within the population tolerance. Very small
  values (\< 1%) can indicate a bottleneck and may lead to a lack of
  diversity.

- **Standard deviation of the log weights**: More variable weights
  (larger s.d.) indicate less efficient sampling. Values greater than 3
  are likely problematic.

- **Maximum unique plans:** an upper bound on the number of unique
  redistricting plans that survive each stage. The percentage in
  parentheses is the ratio of this number to the total number of
  samples. Small values (\< 100) indicate a bottleneck, which leads to a
  loss of sample diversity and a higher variance.

- **Estimated `k` parameter**: How many spanning tree edges were
  considered for cutting at each split. Mostly informational, though
  large jumps may indicate a need to increase `adapt_k_thresh`.

- **Bottleneck**: An asterisk will appear in the right column if a
  bottleneck appears likely, based on the values of the other
  statistics.

In the event of problematic diagnostics, the function will provide
suggestions for improvement.

## Examples

``` r
data(iowa)
iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.1)
plans <- redist_smc(iowa_map, 100)
#> SEQUENTIAL MONTE CARLO
#> Sampling 100 99-unit maps with 4 districts and population between 685430 and 837748.
summary(plans)
#> SMC: 100 sampled plans of 4 districts on 99 units
#> `adapt_k_thresh`=0.99 • `seq_alpha`=0.5
#> `pop_temper`=0
#> 
#> Plan diversity 80% range: 0.46 to 0.80
#> 
#> Sampling diagnostics for SMC run 1 of 1 (100 samples)
#>          Eff. samples (%) Acc. rate Log wgt. sd  Max. unique Est. k 
#> Split 1        98 (98.1%)     14.2%        0.27    64 (101%)     17 
#> Split 2        97 (97.0%)     31.8%        0.35    56 ( 89%)     10 
#> Split 3        94 (94.2%)     14.5%        0.47    55 ( 87%)      6 
#> Resample       77 (77.4%)       NA%        0.47    78 (123%)     NA 
#> 
#> •  Watch out for low effective samples, very low acceptance rates (less than
#> 1%), large std. devs. of the log weights (more than 3 or so), and low numbers
#> of unique plans. R-hat values for summary statistics should be between 1 and
#> 1.05.
```
