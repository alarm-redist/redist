# Parallel Merge-Split/Recombination MCMC Redistricting Sampler

`redist_mergesplit_parallel()` runs
[`redist_mergesplit()`](http://alarm-redist.org/redist/reference/redist_mergesplit.md)
on several chains in parallel.

## Usage

``` r
redist_mergesplit_parallel(
  map,
  nsims,
  chains = 1,
  warmup = if (is.null(init_plan)) 10 else max(100, nsims%/%5),
  thin = 1L,
  init_plan = NULL,
  counties = NULL,
  compactness = 1,
  constraints = list(),
  constraint_fn = function(m) rep(0, ncol(m)),
  adapt_k_thresh = 0.99,
  k = NULL,
  ncores = NULL,
  cl_type = "PSOCK",
  return_all = TRUE,
  init_name = NULL,
  silly_adj_fix = FALSE,
  verbose = FALSE,
  silent = FALSE
)
```

## Arguments

- map:

  A
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object.

- nsims:

  The number of samples to draw, including warmup.

- chains:

  the number of parallel chains to run. Each chain will have `nsims`
  draws. If `init_plan` is sampled, each chain will be initialized with
  its own sampled plan.

- warmup:

  The number of warmup samples to discard. Recommended to be at least
  the first 20% of samples, and in any case no less than around 100
  samples, unless initializing from a random plan.

- thin:

  Save every `thin`-th sample. Defaults to no thinning (1).

- init_plan:

  The initial state of the map, provided as a single vector to be shared
  across all chains, or a matrix with `chains` columns. If not provided,
  will default to the reference map of the map object, or if none
  exists, will sample a random initial state using redist_smc. You can
  also request a random initial state for each chain by setting
  init_plan="sample".

- counties:

  A vector containing county (or other administrative or geographic
  unit) labels for each unit, which may be integers ranging from 1 to
  the number of counties, or a factor or character vector. If provided,
  the algorithm will generate maps tend to follow county lines. There is
  no strength parameter associated with this constraint. To adjust the
  number of county splits further, or to constrain a second type of
  administrative split, consider using
  [`add_constr_splits()`](http://alarm-redist.org/redist/reference/constraints.md),
  [`add_constr_multisplits()`](http://alarm-redist.org/redist/reference/constraints.md),
  and
  [`add_constr_total_splits()`](http://alarm-redist.org/redist/reference/constraints.md).

- compactness:

  Controls the compactness of the generated districts, with higher
  values preferring more compact districts. Must be nonnegative. See the
  'Details' section for more information, and computational
  considerations.

- constraints:

  A list containing information on constraints to implement. See the
  'Details' section for more information.

- constraint_fn:

  A function which takes in a matrix where each column is a
  redistricting plan and outputs a vector of log-weights, which will be
  added the the final weights.

- adapt_k_thresh:

  The threshold value used in the heuristic to select a value `k_i` for
  each splitting iteration. Set to 0.9999 or 1 if the algorithm does not
  appear to be sampling from the target distribution. Must be between 0
  and 1.

- k:

  The number of edges to consider cutting after drawing a spanning tree.
  Should be selected automatically in nearly all cases.

- ncores:

  the number of parallel processes to run. Defaults to the maximum
  available.

- cl_type:

  the cluster type (see `makeCluster()`). Safest is `"PSOCK"`, but
  `"FORK"` may be appropriate in some settings.

- return_all:

  if `TRUE` return all sampled plans; otherwise, just return the final
  plan from each chain.

- init_name:

  a name for the initial plan, or `FALSE` to not include the initial
  plan in the output. Defaults to the column name of the existing plan,
  or "`<init>`" if the initial plan is sampled.

- silly_adj_fix:

  Heuristic for fixing weird inputs.

- verbose:

  Whether to print out intermediate information while sampling.
  Recommended.

- silent:

  Whether to suppress all diagnostic information.

## Value

A
[`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
object with all of the simulated plans, and an additional `chain` column
indicating the chain the plan was drawn from.

## Details

This function draws samples from a specific target measure, controlled
by the `map`, `compactness`, and `constraints` parameters.

Key to ensuring good performance is monitoring the acceptance rate,
which is reported at the sample level in the output. Users should also
check diagnostics of the sample by running
[`summary.redist_plans()`](http://alarm-redist.org/redist/reference/summary.redist_plans.md).

Higher values of `compactness` sample more compact districts; setting
this parameter to 1 is computationally efficient and generates nicely
compact districts.

## References

Carter, D., Herschlag, G., Hunter, Z., and Mattingly, J. (2019). A
merge-split proposal for reversible Monte Carlo Markov chain sampling of
redistricting plans. arXiv preprint arXiv:1911.01503.

McCartan, C., & Imai, K. (2023). Sequential Monte Carlo for Sampling
Balanced and Compact Redistricting Plans. *Annals of Applied Statistics*
17(4). Available at
[doi:10.1214/23-AOAS1763](https://doi.org/10.1214/23-AOAS1763) .

DeFord, D., Duchin, M., and Solomon, J. (2019). Recombination: A family
of Markov chains for redistricting. arXiv preprint arXiv:1911.05725.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fl25)
fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
sampled <- redist_mergesplit_parallel(fl_map, nsims = 100, chains = 100)
} # }
```
