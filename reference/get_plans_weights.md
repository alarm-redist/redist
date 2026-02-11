# Extract the sampling weights from a redistricting simulation.

May be `NULL` if no weights exist (MCMC or optimization methods).

## Usage

``` r
get_plans_weights(plans)

# S3 method for class 'redist_plans'
weights(object, ...)
```

## Arguments

- plans, object:

  the `redist_plans` object

- ...:

  Ignored.

## Value

A numeric vector of weights, with an additional attribute `resampled`
indicating whether the plans have been resampled according to these
weights. If weights have been resampled, this returns the weights before
resampling (i.e., they do not correspond to the resampled plans).

numeric vector
