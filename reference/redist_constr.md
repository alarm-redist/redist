# Set up constraints for sampling

`redist_constr` objects are used to specify constraints when sampling
redistricting plans with
[`redist_smc()`](http://alarm-redist.org/redist/reference/redist_smc.md)
and
[`redist_mergesplit()`](http://alarm-redist.org/redist/reference/redist_mergesplit.md).
Each constraint is specified as a function which scores a given plan.
Higher scores are penalized and sampled less frequently.

## Usage

``` r
redist_constr(map = tibble())
```

## Arguments

- map:

  a
  [`redist_map()`](http://alarm-redist.org/redist/reference/redist_map.md)
  object; the map that will be used in sampling

## Value

a `redist_constr` object, which is just a list with a certain nested
structure.

## Details

The `redist_constr` object keeps track of sampling constraints in a
nested list. You can view the exact structure of this list by calling
[`str()`](https://rdrr.io/r/utils/str.html). Constraints may be added by
using one of the following functions:

- [`add_constr_compet()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_custom()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_edges_rem()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_fry_hold()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_grp_hinge()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_grp_inv_hinge()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_grp_pow()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_incumbency()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_log_st()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_multisplits()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_polsby()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_pop_dev()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_segregation()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_splits()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_status_quo()`](http://alarm-redist.org/redist/reference/constraints.md)

- [`add_constr_total_splits()`](http://alarm-redist.org/redist/reference/constraints.md)

More information about each constraint can be found on the relevant
constraint page.

## Examples

``` r
data(iowa)
map_ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
constr <- redist_constr(map_ia)
constr <- add_constr_splits(constr, strength = 1.5, admin = region)
print(constr)
#> A <redist_constr> with 1 constraint
#> • A splits constraint of strength 1.5
```
