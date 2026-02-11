# Summary plots for `\link{redist_plans}`

If no arguments are passed, defaults to plotting the sampling weights
for the
[`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
object. If no weights exist, plots district populations.

## Usage

``` r
# S3 method for class 'redist_plans'
plot(x, ..., type = "distr_qtys")
```

## Arguments

- x:

  the `redist_plans` object.

- ...:

  passed on to the underlying function

- type:

  the name of the plotting function to use. Will have `redist.plot.`,
  prepended to it; e.g., use `type="plans"` to call
  [`redist.plot.plans`](http://alarm-redist.org/redist/reference/redist.plot.plans.md).
