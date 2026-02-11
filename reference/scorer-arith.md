# Scoring function arithmetic

`redist_scorer` functions may be multiplied by constants and/or added
together to form linear combinations.

## Usage

``` r
# S3 method for class 'redist_scorer'
x * fn2

# S3 method for class 'redist_scorer'
fn1 + fn2

# S3 method for class 'redist_scorer'
fn1 - fn2
```

## Arguments

- x:

  a numeric or a `redist_scorer` function, from
  [`scorers`](http://alarm-redist.org/redist/reference/scorers.md)

- fn2:

  a `redist_scorer` function, from
  [`scorers`](http://alarm-redist.org/redist/reference/scorers.md)

- fn1:

  a `redist_scorer` function, from
  [`scorers`](http://alarm-redist.org/redist/reference/scorers.md)

## Value

function of class redist_scorer
