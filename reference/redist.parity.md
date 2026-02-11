# Calculates Maximum Deviation from Population Parity

Computes the deviation from population parity from a plan. Higher values
indicate that (at least) a single district in the map deviates from
population parity. See Details.

## Usage

``` r
redist.parity(plans, total_pop)

plan_parity(map, .data = pl(), ...)
```

## Arguments

- plans:

  A matrix with one row for each precinct and one column for each map.
  Required.

- total_pop:

  A numeric vector with the population for every precinct.

- map:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- .data:

  a
  [`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
  object

- ...:

  passed on to `redist.parity`

## Value

numeric vector with the population parity for each column

## Details

With a map with `pop` representing the populations of each district, the
deviation from population parity is given as
`max(abs(pop - parity) / parity)` where `parity = sum(pop)/length(pop)`
is the population size for the average district. Therefore, the metric
can be thought of as the maximum percent deviation from equal
population. For example, a value of 0.03 in this metric indicates that
all districts are within 3 percent of population parity.
