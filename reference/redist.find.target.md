# Find Majority Minority Remainder

Given a percent goal for majority minority districts, this computes the
average value of minority in non-majority minority districts. This value
is "tgt_other" in `redist_flip` and `redist_smc`.

## Usage

``` r
redist.find.target(tgt_min, group_pop, total_pop, ndists, nmmd)
```

## Arguments

- tgt_min:

  target group population for majority minority district

- group_pop:

  A vector of populations for some subgroup of interest.

- total_pop:

  A vector containing the populations of each geographic unit.

- ndists:

  The number of congressional districts.

- nmmd:

  The number of majority minority districts.

## Value

numeric value to target
