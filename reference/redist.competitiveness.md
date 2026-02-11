# Compute Competitiveness

Currently only implements the competitiveness function in equation (5)
of Cho & Liu 2016.

## Usage

``` r
competitiveness(map, rvote, dvote, .data = cur_plans())

redist.competitiveness(plans, rvote, dvote, alpha = 1, beta = 1)
```

## Arguments

- map:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- rvote:

  A numeric vector with the Republican vote for each precinct.

- dvote:

  A numeric vector with the Democratic vote for each precinct.

- .data:

  a
  [`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
  object

- plans:

  A numeric vector (if only one map) or matrix with one row for each
  precinct and one column for each map. Required.

- alpha:

  A numeric value for the alpha parameter for the talisman metric

- beta:

  A numeric value for the beta parameter for the talisman metric

## Value

Numeric vector with competitiveness scores

## Examples

``` r
data(fl25)
data(fl25_enum)

plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
# old: comp <- redist.competitiveness(plans_05, fl25$mccain, fl25$obama)
comp <- compet_talisman(plans_05, fl25, mccain, obama)
```
