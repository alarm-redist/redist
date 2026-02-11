# Create County IDs

Create County IDs

## Usage

``` r
redist.county.id(counties)
```

## Arguments

- counties:

  vector of counties, required.

## Value

A vector with an ID that corresponds from 1:n counties

## Examples

``` r
set.seed(2)
counties <- sample(c(rep("a", 20), rep("b", 5)))
redist.county.id(counties)
#>  [1] 2 1 1 2 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 1 1 1 1 2 1
```
