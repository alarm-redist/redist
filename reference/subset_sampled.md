# Subset to sampled or reference draws

Subset to sampled or reference draws

## Usage

``` r
subset_sampled(plans, matrix = TRUE)

subset_ref(plans, matrix = TRUE)
```

## Arguments

- plans:

  the `redist_plans` object

- matrix:

  if `TRUE`, the default, also subset the plans matrix. If the plans
  matrix is not needed, turning this off may save some time.

## Value

a `redist_plans` object, with only rows corresponding to simulated (or
reference) draws remaining.
