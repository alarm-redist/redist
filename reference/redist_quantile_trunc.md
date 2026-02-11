# Helper function to truncate importance weights

Defined as `pmin(x, quantile(x, 1 - length(x)^(-0.5)))`

## Usage

``` r
redist_quantile_trunc(x)
```

## Arguments

- x:

  the weights

## Value

numeric vector

## Examples

``` r
redist_quantile_trunc(c(1, 2, 3, 4))
#> [1] 1.0 2.0 2.5 2.5
```
