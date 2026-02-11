# Redistricting Optimization through Short Bursts

This function uses
[`redist_mergesplit()`](http://alarm-redist.org/redist/reference/redist_mergesplit.md)
or
[`redist_flip()`](http://alarm-redist.org/redist/reference/redist_flip.md)
to optimize a redistrict plan according to a user-provided criteria. It
does so by running the Markov chain for "short bursts" of usually 10
iterations, and then starting the chain anew from the best plan in the
burst, according to the criteria. This implements the ideas in the
below-referenced paper, "Voting Rights, Markov Chains, and Optimization
by Short Bursts."

## Usage

``` r
redist_shortburst(
  map,
  score_fn = NULL,
  stop_at = NULL,
  burst_size = ifelse(backend == "mergesplit", 10L, 50L),
  max_bursts = 500L,
  maximize = TRUE,
  init_plan = NULL,
  counties = NULL,
  constraints = redist_constr(map),
  compactness = 1,
  adapt_k_thresh = 0.95,
  reversible = TRUE,
  fixed_k = NULL,
  return_all = TRUE,
  thin = 1L,
  backend = "mergesplit",
  flip_lambda = 0,
  flip_eprob = 0.05,
  verbose = TRUE
)
```

## Arguments

- map:

  A [redist_map](http://alarm-redist.org/redist/reference/redist_map.md)
  object.

- score_fn:

  A function which takes a matrix of plans and returns a score (or,
  generally, a row vector) for each plan. Can also be a purrr-style
  anonymous function. See
  [`?scorers`](http://alarm-redist.org/redist/reference/scorers.md) for
  some function factories for common scoring rules.

- stop_at:

  A threshold to stop optimization at. When `score_fn` returns a row
  vector per plan, `maximize` can be an equal-length vector specifying a
  threshold for each dimension, which must all be met for the algorithm
  to stop.

- burst_size:

  The size of each burst. 10 is recommended for the `mergesplit` backend
  and 50 for the `flip` backend. Can also provide burst schedule
  function which takes the current iteration (an integer) and returns
  the desired burst size. This can be a random function.

- max_bursts:

  The maximum number of bursts to run before returning.

- maximize:

  If `TRUE`, try to maximize the score; otherwise, try to minimize it.
  When `score_fn` returns a row vector per plan, `maximize` can be an
  equal-length vector specifying whether each dimension should be
  maximized or minimized.

- init_plan:

  The initial state of the map. If not provided, will default to the
  reference map of the `map` object, or if none exists, will sample a
  random initial state using
  [`redist_smc()`](http://alarm-redist.org/redist/reference/redist_smc.md).
  You can also request a random initial state by setting
  `init_plan="sample"`.

- counties:

  A vector containing county (or other administrative or geographic
  unit) labels for each unit, which may be integers ranging from 1 to
  the number of counties, or a factor or character vector. If provided,
  the algorithm will only generate maps which split up to `ndists-1`
  counties. If no county-split constraint is desired, this parameter
  should be left blank.

- constraints:

  A `redist_constr` with Gibbs constraints.

- compactness:

  Controls the compactness of the generated districts, with higher
  values preferring more compact districts. Must be non-negative. See
  [`redist_mergesplit`](http://alarm-redist.org/redist/reference/redist_mergesplit.md)
  for more information.

- adapt_k_thresh:

  The threshold value used in the heuristic to select a value `k_i` for
  each splitting iteration.

- reversible:

  If `FALSE` and `backend="mergesplit"`, the Markov chain used will not
  be reversible. This may speed up optimization.

- fixed_k:

  If not `NULL`, will be used to set the `k` parameter for the
  `mergesplit` backend. If e.g. `k=1` then the best edge in each
  spanning tree will be used. Lower values may speed up optimization at
  the cost of the Markov chain no longer targeting a known distribution.
  Recommended only in conjunction with `reversible=FALSE`.

- return_all:

  Whether to return all the burst results or just the best one
  (generally, the Pareto frontier). Recommended for monitoring purposes.

- thin:

  Save every `thin`-th sample. Defaults to no thinning (1). Ignored if
  `return_all=TRUE`.

- backend:

  the MCMC algorithm to use within each burst, either "mergesplit" or
  "flip".

- flip_lambda:

  The parameter determining the number of swaps to attempt each
  iteration of flip mcmc. The number of swaps each iteration is equal to
  Pois(lambda) + 1. The default is 0.

- flip_eprob:

  The probability of keeping an edge connected in flip mcmc. The default
  is 0.05.

- verbose:

  Whether to print out intermediate information while sampling.
  Recommended for monitoring purposes.

## Value

a `redist_plans` object containing the final best plan (or the best
plans after each burst, if `return_all=TRUE`.

## References

Cannon, S., Goldbloom-Helzner, A., Gupta, V., Matthews, J. N., & Suwal,
B. (2020). Voting Rights, Markov Chains, and Optimization by Short
Bursts. arXiv preprint arXiv:2011.02288.

## Examples

``` r
# \donttest{
data(iowa)

iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
redist_shortburst(iowa_map, scorer_frac_kept(iowa_map), max_bursts = 50)
#> MERGE-SPLIT SHORT BURSTS
#> Sampling up to 50 bursts of 10 iterations each.
#> Burst  Improve?    score 
#>     5            0.788288
#>    10            0.788288
#>    15     🥳     0.810811
#>    20            0.810811
#>    25            0.810811
#>    30            0.810811
#>    35            0.810811
#>    36     🙂     0.819820
#>    40            0.819820
#>    44     😀     0.824324
#>    45            0.824324
#>    46     🎇     0.828829
#>    50            0.828829
#> A <redist_plans> containing 50 sampled plans and 1 reference plan
#> Plans have 4 districts from a 99-unit map, and were drawn using short bursts.
#> Plans matrix: int [1:99, 1:51] 1 1 2 3 4 2 2 4 2 2 ...
#> # A tibble: 204 × 5
#>    draw   district total_pop score burst_size
#>    <fct>     <int>     <dbl> <dbl>      <int>
#>  1 <init>        1    761612 0.788         NA
#>  2 <init>        2    761548 0.788         NA
#>  3 <init>        3    761624 0.788         NA
#>  4 <init>        4    761571 0.788         NA
#>  5 1             1    761612 0.788         10
#>  6 1             2    761548 0.788         10
#>  7 1             3    761624 0.788         10
#>  8 1             4    761571 0.788         10
#>  9 2             1    761612 0.788         10
#> 10 2             2    761548 0.788         10
#> # ℹ 194 more rows
redist_shortburst(iowa_map, ~ 1 - scorer_frac_kept(iowa_map)(.), max_bursts = 50)
#> MERGE-SPLIT SHORT BURSTS
#> Sampling up to 50 bursts of 10 iterations each.
#> Burst  Improve?    score 
#>     5            0.211712
#>     7     🎆     0.216216
#>    10            0.216216
#>    15            0.216216
#>    19     🌈     0.265766
#>    20            0.265766
#>    23     🙂     0.274775
#>    25            0.274775
#>    30            0.274775
#>    35            0.274775
#>    37     😎     0.279279
#>    38     🎉     0.288288
#>    40            0.288288
#>    45            0.288288
#>    50            0.288288
#> A <redist_plans> containing 50 sampled plans and 1 reference plan
#> Plans have 4 districts from a 99-unit map, and were drawn using short bursts.
#> Plans matrix: int [1:99, 1:51] 1 1 2 3 4 2 2 4 2 2 ...
#> # A tibble: 204 × 5
#>    draw   district total_pop score burst_size
#>    <fct>     <int>     <dbl> <dbl>      <int>
#>  1 <init>        1    761612 0.212         NA
#>  2 <init>        2    761548 0.212         NA
#>  3 <init>        3    761624 0.212         NA
#>  4 <init>        4    761571 0.212         NA
#>  5 1             1    761612 0.212         10
#>  6 1             2    761548 0.212         10
#>  7 1             3    761624 0.212         10
#>  8 1             4    761571 0.212         10
#>  9 2             1    761612 0.212         10
#> 10 2             2    761548 0.212         10
#> # ℹ 194 more rows
# }
```
