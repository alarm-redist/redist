# Add a reference plan to a set of plans

This function facilitates comparing an existing (i.e., non-simulated)
redistricting plan to a set of simulated plans.

## Usage

``` r
add_reference(plans, ref_plan, name = NULL)
```

## Arguments

- plans:

  a `redist_plans` object

- ref_plan:

  an integer vector containing the reference plan. It will be renumbered
  to 1..`ndists`.

- name:

  a human-readable name for the reference plan. Defaults to the name of
  `ref_plan`.

## Value

a modified `redist_plans` object containing the reference plan
