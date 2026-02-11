# 'Flip' Markov Chain Monte Carlo Redistricting Simulation (Fifield et al. 2020)

This function allows users to simulate redistricting plans using a
Markov Chain Monte Carlo algorithm (Fifield, Higgins, Imai, and Tarr
2020). Several constraints corresponding to substantive requirements in
the redistricting process are implemented, including population parity
and geographic compactness. In addition, the function includes
multiple-swap and simulated tempering functionality to improve the
mixing of the Markov Chain.

## Usage

``` r
redist_flip(
  map,
  nsims,
  warmup = 0,
  init_plan,
  constraints = add_constr_edges_rem(redist_constr(map), 0.4),
  thin = 1,
  eprob = 0.05,
  lambda = 0,
  temper = FALSE,
  betaseq = "powerlaw",
  betaseqlength = 10,
  betaweights = NULL,
  adapt_lambda = FALSE,
  adapt_eprob = FALSE,
  exact_mh = FALSE,
  adjswaps = TRUE,
  init_name = NULL,
  verbose = TRUE,
  nthin
)
```

## Arguments

- map:

  A
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object.

- nsims:

  The number of samples to draw, not including warmup.

- warmup:

  The number of warmup samples to discard.

- init_plan:

  A vector containing the congressional district labels of each
  geographic unit. The default is `NULL`. If not provided, a random
  initial plan will be generated using `redist_smc`. You can also
  request to initialize using `redist.rsg` by supplying 'rsg', though
  this is not recommended behavior.

- constraints:

  A `redist_constr` object.

- thin:

  The amount by which to thin the Markov Chain. The default is `1`.

- eprob:

  The probability of keeping an edge connected. The default is `0.05`.

- lambda:

  lambda The parameter determining the number of swaps to attempt each
  iteration of the algorithm. The number of swaps each iteration is
  equal to Pois(`lambda`) + 1. The default is `0`.

- temper:

  Whether to use simulated tempering algorithm. Default is FALSE.

- betaseq:

  Sequence of beta values for tempering. The default is `powerlaw` (see
  Fifield et. al (2020) for details).

- betaseqlength:

  Length of beta sequence desired for tempering. The default is `10`.

- betaweights:

  betaweights Sequence of weights for different values of beta. Allows
  the user to upweight certain values of beta over others. The default
  is `NULL` (equal weighting).

- adapt_lambda:

  adapt_lambda Whether to adaptively tune the lambda parameter so that
  the Metropolis-Hastings acceptance probability falls between 20% and
  40%. Default is FALSE.

- adapt_eprob:

  eprob Whether to adaptively tune the edgecut probability parameter so
  that the Metropolis-Hastings acceptance probability falls between 20%
  and 40%. Default is FALSE.

- exact_mh:

  Whether to use the approximate (FALSE) or exact (TRUE)
  Metropolis-Hastings ratio calculation for accept-reject rule. Default
  is FALSE.

- adjswaps:

  Flag to restrict swaps of beta so that only values adjacent to current
  constraint are proposed. The default is `TRUE`.

- init_name:

  a name for the initial plan, or `FALSE` to not include the initial
  plan in the output. Defaults to the column name of the existing plan,
  or "`<init>`" if the initial plan is sampled.

- verbose:

  Whether to print initialization statement. Default is `TRUE`.

- nthin:

  Deprecated. Use `thin`.

## Value

A
[`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
object containing the simulated plans.

## Details

`redist_flip` allows for Gibbs constraints to be supplied via a list
object passed to `constraints`. `redist_flip` uses a small compactness
constraint by default, as this improves the realism of the maps greatly
and also leads to large speed improvements. (One of the most time
consuming aspects of the flip MCMC backend is checking for district
shattering, which is slowed down even further by non-compact districts.
As such, it is recommended that all flip simulations use at least a
minimal compactness constraint, even if you weaken it from the default
settings.) The default is a `compact` constraint using the
`edges-removed` metric with a weight of 0.6. For very small maps (\< 100
precincts), you will likely want to weaken (lower) this constraint,
while for very large maps (\> 5000 precincts), you will likely want to
strengthen (increase) this constraint. Otherwise, for most maps, the
default constraint should be a good starting place.

`redist_flip` samples from a known target distribution which can be
described using the `constraints`. The following describes the
constraints available. The general advice is to set weights in a way
that gets between 20% and 40% acceptance on average, though more tuning
advice is available in the vignette on using MCMC methods.Having too
small of an acceptance rate indicates that the weights within
`constraints` are too large and will impact sampling efficiency. If the
Metropolis Hastings acceptance rate is too large, this may impact the
target distribution, but may be fine for general exploration of possible
maps.

There are currently 9 implemented constraint types, though \``compact`
and `partisan` have sub-types which are specified via a character
`metric` within their respective list objects. The constraints are as
follows:

- `compact` - biases the algorithm towards drawing more compact
  districts.

- weight - the coefficient to put on the Gibbs constraint

- metric - which metric to use. Must be one of `edges-removed` (the
  default), `polsby-popper`, `fryer-holden`, or `log-st`. Using Polsby
  Popper is generally not recommended, as `edges-removed` is faster and
  highly correlated. `log-st` can be used to match the target
  distribution of `redist_smc` or `redist_mergesplit`.

- areas - Only used with `polsby-popper` - A vector of precinct areas.

- borderlength_mat - Only used with `polsby-popper` - A matrix of
  precinct border lengths.

- ssdmat - Only used with `fryer-holden` - A matrix of squared distances
  between precinct centroids.

- ssd_denom - Only used with `fryer-holden` - a positive integer to use
  as the normalizing constant for the Relative Proximity Index.

- `population` - A Gibbs constraint to complement the hard population
  constraint set by `pop_tol`. This penalizes moves which move away from
  smaller population parity deviations. It is very useful when an
  `init_plan` sits outside of the desired `pop_tol` but there are
  substantive reasons to use that plan. This constraint uses the input
  to `total_pop`.

- weight - the coefficient to put on the Gibbs constraint

- `countysplit` This is a Gibbs constraint to minimize county splits.
  Unlike SMC's county constraint, this allows for more than `ndists - 1`
  splits and does not require that counties are contiguous.

- weight - the coefficient to put on the Gibbs constraint

- `hinge` This uses the proportion of a group in a district and matches
  to the nearest target proportion, and then creates a penalty of
  \\\sqrt{max(0, nearest.target - group.pct)}\\.

- weight - the coefficient to put on the Gibbs constraint

- minorityprop - A numeric vector of minority proportions (between 0
  and 1) which districts should aim to have

- `vra` This takes two target proportions of the presence of a minority
  group within a district. \\(\|target.min - group.pct\|\|target.other -
  group.pct\|)^{1.5})\\

- weight - the coefficient to put on the Gibbs constraint

- target_min - the target minority percentage. Often, this is set to
  0.55 to encourage minority majority districts.

- target_other - the target minority percentage for non majority
  minority districts.

- `minority` This constraint sorts the districts by the proportion of a
  group in a district and compares the highest districts to the entries
  of minorityprop. This takes the form \\\sum\_{i=1}^{n}
  \sqrt{\|group.pct(i) - minorityprop(i)\| }\\ where n is the length of
  minorityprop input.

- weight - the coefficient to put on the Gibbs constraint

- minorityprop - A numeric vector of minority proportions (between 0
  and 1) which districts should aim to have

- `similarity` This is a status-quo constraint which penalizes plans
  which are very different from the starting place. It is useful for
  local exploration.

- weight - the coefficient to put on the Gibbs constraint

- `partisan` This is a constraint which minimizes partisan bias, either
  as measured as the difference from proportional representation or as
  the magnitude of the efficiency gap.

- weight - the coefficient to put on the Gibbs constraint

- rvote - An integer vector of votes for Republicans or other party

- dvote - An integer vector of votes for Democrats or other party

- metric - which metric to use. Must be one of
  `proportional-representation` or `efficiency-gap`.

- `segregation` This constraint attempts to minimize the degree of
  dissimilarity between districts by group population.

- weight - the coefficient to put on the Gibbs constraint

## References

Fifield, B., Higgins, M., Imai, K., & Tarr, A. (2020). Automated
redistricting simulation using Markov chain Monte Carlo. *Journal of
Computational and Graphical Statistics*, 29(4), 715-728.

## Examples

``` r
data(iowa)
iowa_map <- redist_map(iowa, ndists = 4, existing_plan = cd_2010, total_pop = pop,
    pop_tol = 0.05)
sims <- redist_flip(map = iowa_map, nsims = 100)
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> ℹ Starting swMH().
#> ■                                  1% | ETA: 0s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.68
#> 
```
