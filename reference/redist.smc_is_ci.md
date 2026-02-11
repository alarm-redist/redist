# (Deprecated) Confidence Intervals for Importance Sampling Estimates

Builds a confidence interval for a quantity of interest, given
importance sampling weights.

## Usage

``` r
redist.smc_is_ci(x, wgt, conf = 0.99)
```

## Arguments

- x:

  A numeric vector containing the quantity of interest

- wgt:

  A numeric vector containing the nonnegative importance weights. Will
  be normalized automatically.

- conf:

  The confidence level for the interval.

## Value

A two-element vector of the form `[lower, upper]` containing the
importance sampling confidence interval.
