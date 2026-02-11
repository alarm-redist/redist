# Create Constraints for SMC

Create Constraints for SMC

## Usage

``` r
redist.constraint.helper(
  constraints = "vra",
  tgt_min = 0.55,
  group_pop,
  total_pop,
  ndists,
  nmmd,
  strength_vra = 2500,
  pow_vra = 1.5
)
```

## Arguments

- constraints:

  Vector of constraints to include. Currently only 'vra' implemented.

- tgt_min:

  Defaults to 0.55. If 'vra' included, the minority percent to encourage
  in each district.

- group_pop:

  A vector of populations for some subgroup of interest.

- total_pop:

  A vector containing the populations of each geographic unit.

- ndists:

  The total number of districts.

- nmmd:

  The number of majority minority districts to target for 'vra'
  constraint

- strength_vra:

  The strength of the 'vra' constraint. Defaults to 2500.

- pow_vra:

  The exponent for the 'vra' constraint. Defaults to 1.5.

## Value

list of lists for each constraint selected
