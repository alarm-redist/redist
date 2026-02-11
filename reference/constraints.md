# Sampling constraints

The
[`redist_smc()`](http://alarm-redist.org/redist/reference/redist_smc.md)
and
[`redist_mergesplit()`](http://alarm-redist.org/redist/reference/redist_mergesplit.md)
algorithms in this package allow for additional constraints on the
redistricting process to be encoded in the target distribution for
sampling. These functions are provided to specify these constraints. All
arguments are quoted and evaluated in the context of the data frame
provided to
[`redist_constr()`](http://alarm-redist.org/redist/reference/redist_constr.md).

## Usage

``` r
add_constr_status_quo(constr, strength, current)

add_constr_grp_pow(
  constr,
  strength,
  group_pop,
  total_pop = NULL,
  tgt_group = 0.5,
  tgt_other = 0.5,
  pow = 1
)

add_constr_grp_hinge(
  constr,
  strength,
  group_pop,
  total_pop = NULL,
  tgts_group = c(0.55)
)

add_constr_grp_inv_hinge(
  constr,
  strength,
  group_pop,
  total_pop = NULL,
  tgts_group = c(0.55)
)

add_constr_compet(constr, strength, dvote, rvote, pow = 0.5)

add_constr_incumbency(constr, strength, incumbents)

add_constr_splits(constr, strength, admin)

add_constr_multisplits(constr, strength, admin)

add_constr_total_splits(constr, strength, admin)

add_constr_pop_dev(constr, strength)

add_constr_segregation(constr, strength, group_pop, total_pop = NULL)

add_constr_polsby(constr, strength, perim_df = NULL)

add_constr_fry_hold(
  constr,
  strength,
  total_pop = NULL,
  ssdmat = NULL,
  denominator = 1
)

add_constr_log_st(constr, strength, admin = NULL)

add_constr_edges_rem(constr, strength)

add_constr_custom(constr, strength, fn)
```

## Arguments

- constr:

  A
  [`redist_constr()`](http://alarm-redist.org/redist/reference/redist_constr.md)
  object

- strength:

  The strength of the constraint. Higher values mean a more restrictive
  constraint.

- current:

  The reference map for the status quo constraint.

- group_pop:

  A vector of group population

- total_pop:

  A vector of total population. Defaults to the population vector used
  for sampling.

- tgt_group, tgt_other:

  Target group shares for the power-type constraint.

- pow:

  The exponent for the power-type constraint.

- tgts_group:

  A vector of target group shares for the hinge-type constraint.

- dvote, rvote:

  A vector of Democratic or Republican vote counts

- incumbents:

  A vector of unit indices for incumbents. For example, if three
  incumbents live in the precincts that correspond to rows 1, 2, and 100
  of your
  [redist_map](http://alarm-redist.org/redist/reference/redist_map.md),
  entering incumbents = c(1, 2, 100) would avoid having two or more
  incumbents be in the same district.

- admin:

  A vector indicating administrative unit membership

- perim_df:

  A dataframe output from
  [`redistmetrics::prep_perims`](http://alarm-redist.org/redistmetrics/reference/prep_perims.md)

- ssdmat:

  Squared distance matrix for Fryer Holden constraint

- denominator:

  Fryer Holden minimum value to normalize by. Default is 1 (no
  normalization).

- fn:

  A function

## Details

All constraints are fed into a Gibbs measure, with coefficients on each
constraint set by the corresponding `strength` parameter. The strength
can be any real number, with zero corresponding to no constraint. Higher
and higher `strength` values will eventually cause the algorithm's
accuracy and efficiency to suffer. Whenever you use constraints, be sure
to check all sampling diagnostics.

The `status_quo` constraint adds a term measuring the variation of
information distance between the plan and the reference, rescaled to
\[0, 1\].

The `grp_hinge` constraint takes a list of target group percentages. It
matches each district to its nearest target percentage, and then applies
a penalty of the form \\\sqrt{max(0, tgt - grouppct)}\\, summing across
districts. This penalizes districts which are below their target
percentage. Use
[`plot.redist_constr()`](http://alarm-redist.org/redist/reference/plot.redist_constr.md)
to visualize the effect of this constraint and calibrate `strength`
appropriately.

The `grp_inv_hinge` constraint takes a list of target group percentages.
It matches each district to its nearest target percentage, and then
applies a penalty of the form \\\sqrt{max(0, grouppct - tgt)}\\, summing
across districts. This penalizes districts which are above their target
percentage. Use
[`plot.redist_constr()`](http://alarm-redist.org/redist/reference/plot.redist_constr.md)
to visualize the effect of this constraint and calibrate `strength`
appropriately.

The `grp_pow` constraint (for expert use) adds a term of the form
\\(\|tgtgroup-grouppct\|\|tgtother-grouppct\|)^{pow})\\, which
encourages districts to have group shares near either `tgt_group` or
`tgt_other`. Values of `strength` depend heavily on the values of these
parameters and especially the `pow` parameter. Use
[`plot.redist_constr()`](http://alarm-redist.org/redist/reference/plot.redist_constr.md)
to visualize the effect of this constraint and calibrate `strength`
appropriately.

The `compet` constraint encourages competitiveness by applying the
`grp_pow` constraint with target percentages set to 50%. For
convenience, it is specified with Democratic and Republican vote shares.

The `incumbency` constraint adds a term counting the number of districts
containing paired-up incumbents. Values of `strength` should generally
be small, given that the underlying values are counts.

The `splits` constraint adds a term counting the number of counties
which are split once or more. Values of `strength` should generally be
small, given that the underlying values are counts.

The `multisplits` constraint adds a term counting the number of counties
which are split twice or more. Values of `strength` should generally be
small, given that the underlying values are counts.

The `total_splits` constraint adds a term counting the total number of
times each county is split, summed across counties (i.e., counting the
number of excess district-county pairs). Values of `strength` should
generally be small, given that the underlying values are counts.

The `edges_rem` constraint adds a term counting the number of edges
removed from the adjacency graph. This is only usable with
[`redist_flip()`](http://alarm-redist.org/redist/reference/redist_flip.md),
as other algorithms implicitly use this via the `compactness` parameter.
Values of `strength` should generally be small, given that the
underlying values are counts.

The `log_st` constraint constraint adds a term counting the log number
of spanning trees. This is only usable with
[`redist_flip()`](http://alarm-redist.org/redist/reference/redist_flip.md),
as other algorithms implicitly use this via the `compactness` parameter.

The `polsby` constraint adds a term encouraging compactness as defined
by the Polsby Popper metric. Values of `strength` may be of moderate
size.

The `fry_hold` constraint adds a term encouraging compactness as defined
by the Fryer Holden metric. Values of `strength` should be extremely
small, as the underlying values are massive when the true minimum Fryer
Holden denominator is not known.

The `segregation` constraint adds a term encouraging segregation among
minority groups, as measured by the dissimilarity index.

The `pop_dev` constraint adds a term encouraging plans to have smaller
population deviations from the target population.

The `custom` constraint allows the user to specify their own constraint
using a function which evaluates districts one at a time. The provided
function `fn` should take two arguments: a vector describing the current
plan assignment for each unit as its first argument, and an integer
describing the district which to evaluate in the second argument.
`which([plans == distr])` would give the indices of the units that are
assigned to a district `distr` in any iteration. The function must
return a single scalar for each plan - district combination, where a
value of 0 indicates no penalty is applied. If users want to penalize an
entire plan, they can have the penalty function return a scalar that
does not depend on the district. It is important that `fn` not use
information from precincts not included in `distr`, since in the case of
SMC these precincts may not be assigned any district at all (`plan` will
take the value of 0 for these precincts). The flexibility of this
constraint comes with an additional computational cost, since the other
constraints are written in C++ and so are more performant.

## Examples

``` r
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
constr <- redist_constr(iowa_map)
constr <- add_constr_splits(constr, strength = 1.5, admin = name)
constr <- add_constr_grp_hinge(constr, strength = 100,
    dem_08, tot_08, tgts_group = c(0.5, 0.6))
# encourage districts to have the same number of counties
constr <- add_constr_custom(constr, strength = 1000, fn = function(plan, distr) {
    # notice that we only use information on precincts in `distr`
    abs(sum(plan == distr) - 99/4)
})
print(constr)
#> A <redist_constr> with 3 constraints
#> • A splits constraint of strength 1.5
#> • A (hinge-type) group share constraint of strength 100
#>     group_pop : num [1:99] 1924 1118 3971 2970 1739 ...
#>     total_pop : num [1:99] 4053 2206 7059 6176 3435 ...
#>     tgts_group: num [1:2] 0.5 0.6
#> • A custom constraint of strength 1000
#>     fn:function (plan, distr)  
```
