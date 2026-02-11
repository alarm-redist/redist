# Tally a variable by district

Tally a variable by district

## Usage

``` r
tally_var(map, x, .data = pl())
```

## Arguments

- map:

  a `redist_map` object

- x:

  a variable to tally. Tidy-evaluated.

- .data:

  a `redist_plans` object or matrix of plans

## Value

a vector containing the tallied values by district and plan
(column-major)
