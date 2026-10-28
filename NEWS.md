# redist 5.0.0
* Allows for sampling multimember districting plans from both SMC and 
mergesplit MCMC. For SMC the district seat sizes (how many 
legislators a district can have) must be a range of values e.g. (3, 4, 5, 6).
MCMC supports non-continuous ranges. 
    * `redist_map` now takes `nseats` and `seats_range` arguments for multi-member maps
    * `redist_plans` now store a `seats` column for multi-member and partial plans. 
* Introduces new methods (different `sampling_space` arguments) for sampling plans 
in both `redist_smc` and `redist_mergesplit`. 
The final output is still plans from the same distribution as before but new 
sampling spaces and their associated splitting methods tend to perform better. 
The `adapt_k_thresh` parameter is not needed for these new sampling
spaces.
* Improves SMC and Mergesplit MCMC performance by pre-allocating and reusing as 
much memory as possible while drawing spanning trees.
* Replaces old SMC weights with generally lower variance optimal weights.
* Adds the option to add Mergesplit MCMC steps at any point during an SMC run 
using the new `ms_params` argument.  
Adding mergesplit steps can help boost convergence for plans with a larger 
number of districts without increasing the sample size. Maps with 10 or more districts
will add mergesplit steps by default.
* Introduces a new method for splitting in SMC - generalized region splits.
Instead of splitting off one district at a time this allows for splitting into
two arbitrary sized regions. This tends to reduce the runtime for a fixed
sample size since, on average, it draws spanning trees on smaller subgraphs then 
single district splits. Currently only supported for single-member maps.
* When counties are used `redist_mergesplit` now samples from the same target 
distribution as `redist_smc` (it guarantees no more than the number of districts
minus 1 splits). For more details 
see O'Sullivan, McCartan, and Imai (2026) <doi:10.48550/arXiv.2603.22188>. 
* `redist_mergesplit_parallel` is deprecated. `redist_mergesplit` can be run in
parallel using the `ncores` argument. 
* `redist_mergesplit` inputs now work differently. 
    * `nsims` is now the number of plans saved. 
    * `warmup` is the number of steps to run the chain for before collecting any samples.
    * `thin` means we will run the chain for `thin - 1` steps between saving plans
    * Overall the chain will be run for `warmup + nsims * thin` and return `nsims` plans.
    * `adapt_k_thresh` must now be passed in as an argument in `split_params` 
    and it is only needed when `sampling_space = graph_plan`.
* `redist_smc` inputs now work differently
    * `ncores` now only controls the number of threads spawned in c++ for a single
    R process and defaults to single-threading. The number of independent R 
    processes can be controlled by the `nproc` argument in `control`.
    * `adapt_k_thresh` must now be passed in as an argument in `split_params` 
    and it is only needed when `sampling_space = graph_plan`.
    * `truncate`, `trunc_fn`, and `final_infl` are no longer accepted as arguments.
    * Partial sampling starting from an `init_particles` input now expects 1-indexed 
    values instead of 0-indexed.
    * `seq_alpha` cannot be changed from 1 when mergesplit steps are used.
* Normalizing constants can now be estimated from `redist_smc` runs with the
`est_norm_unbiased` and `est_norm_biased` functions. This is only supported
when `seq_alpha` is set to 1.
* The diagnostic-checking function `summary.redist_plans()` has new warnings 
and recommendations for how to fix samples from `redist_smc` along with faster
rhat computation. 
* New functionality for `redist_constr` has been added:
    * All constraints now allow for rejection sampling in `redist_smc` and 
    `redist_mergesplit`. All constraints now take an optional threshold 
    argument `thresh` where for a newly split plan if the score is
    greater than or equal to `thresh` then the plan will be automatically rejected. 
    This amounts to giving plans where any region has a score above `thresh` a 
    probability of 0. 
    * Constraints now score all regions, not just districts by default. 
    What sizes can be scored is controlled by the new `only_nseats` argument.
    * Adds the following new plan-wide constraints which are a function of 
    the entire plan as opposed to individual districts:
      * `add_constr_custom_plan`
      * `add_constr_min_group_frac`
      * `add_constr_plan_incumbency`
      * `add_constr_plan_splits`
      * `add_constr_total_plan_splits`
    * When in the splitting process constraints are applied can be adjusted with
    the new `only_nregions` argument. 
    * `add_constr_edges_rem` is now supported for `redist_smc` and `redist_mergesplit`
* The mergesplit backend for `redist_shortburst` now uses uniform edge sampling 
with forest space for the backend instead of sampling with graph space and all
`k` related parameters have been removed. 

# redist 4.4.0
* Old `enumpart` functions are deprecated in favor of `redist_enumpart()` which provides a more consistent syntax.

# redist 4.3.2
* Allows for parallel flip with `chains` argument in `redist_flip()`.
* Fixes URL issues causing a note on CRAN. We have pointed the links to our website directly.

# redist 4.3.0
* Improves SMC performance by pre-allocating some memory while drawing spanning trees.
* Replaces SMC label-counting adjustments (exact and importance-sampling-based) with a new backward kernel that eliminates approximation error and requires far less computation
* 4.2.0 introduced some regressions in `redist_shortburst()` along with the new features. The following issues are fixed: 
  * the returned `redist_plans` object would store the wrong score for the '<init>' input. This issue only impacted the returned object and the correct score was used during the shortburst process. (#180)
  * the function would return immediately if `stop_at` was specified and `minimize = FALSE`. (#181) 
  
 * Add `summary()` support for plans sampled with the `flip` algorithm. This does not replace the full flip diagnostic suite, but provides an easy way to compute r-hats.

# redist 4.2.0
* Deprecate functionality that is provided by `redistmetrics` package.
* Improve contiguity checking speed drastically.
* Support for multiple independent scoring functions in `redist_shortburst()`.
With multiple scorers, the algorithm will stochastically explore to try to 
find the largest Pareto frontier for the scores. The frontier can be accessed with
`attr(<plans obj>, "pareto_score")`.
* Removes the MPI vignette which relied on older implementations of `redist.mcmc()`, which was replaced by `redist.flip()` a few years ago, and finally `redist_flip()`.

# redist 4.1.1
* Resolves a sanitizer error for CRAN

# redist 4.1.0
* Improved diagnostic output
* New `redist_ci` interface for confidence interval calculation
* Improved plotting options with `redist.plot.distr_qtys()` for custom geometry types.
* Improved resampling efficiency at the final SMC stage
* Faster implementation of loop-erased random walk in C++
* Faster random number generation in C++
* Updated citation information

# redist 4.0.0
* A new constraint interface that is more flexible, user friendly, and consistent 
across algorithms (see `redist_constr()` and `?constraints`). For the first time,
user-defined custom constraints are supported and integrated within all three 
algorithms.
* New diagnostic-checking function, `summary.redist_plans()`
* Summary statistics have been broken out into a new `redistmetrics` package
This will speed up compilation time and also provides a cleaner, more extensible 
interface for the implementation of additional metrics.
* Parallel computing support for the SMC algorithm, both within and across sampling runs
* Reproducible across-run parallelism throughout the package, via `doRNG`
* Much faster `match_numbers()` using the Hungarian method
* `min_move_parity()` calculates how much population needs to be moved between 
districts in order to completely balance a redistricting plan.
* Support for partial SMC simulations, where fewer districts are drawn than the 
total number. Allows advanced users to manually combine partial runs to 
form complete maps.
* Improved algorithm reporting, including new progress bars and `cli` errors and 
warnings throughout the package
* Update the SMC algorithm to include a missing correction factor for the number
of ways to sequentially label districts. This factor should not have an effect
on substantive conclusions and summary statistics.
* Remove deprecated functions
* Many bug fixes (see https://github.com/alarm-redist/redist/issues)

# redist 3.1.6
* Utilities for using municipalities as well as counties in split calculations

# redist 3.1.5
* skip SMC test on Linux

# redist 3.1.4
* skip SMC test on Solaris

# redist 3.1.2
* Fixes crash caused by `redist.splits()`

# redist 3.1.1
* Fixes printing bug in `color_graph()`

# redist 3.1.0
* Removes prior deprecated functions and arguments
* Fix bugs (#78, #81, #86)
* Introduces `redist_mergesplit_parallel()`
* Adds `rbind()` generic for `redist_plans` objects
* Improves sampling speed for SMC and Merge-split with county constraint
* Adds county split measures.
* Adds population overlap measures for plan comparisons.
* Deprecates `redist.smc()` in favor of `redist_smc()` and `redist.mergesplit()` in favor of `redist_mergesplit()`.
# redist 3.0.2
* Fix bugs (#60, #61, #62, #70, #71, #72), including s2 compatibility, Solaris fixes, and improved dplyr verb robustness.

# redist 3.0.1

* New tidy interface, including new `redist_map` and `redist_plans` objects
* Merge-split MCMC now available in `redist_mergesplit()`
* Short burst MCMC optimization now available in `redist_shortburst()` along
  with scoring functions (`?scorers`)
* Improved Flip MCMC interface and performance improvements
* New support for larger simulation size limits
* Functions to freeze parts of a map and extract district cores
* New VRA constraint
* Many new plotting functions
* Consistent function and argument names
* New partisanship and compactness metrics
* Performance improvements to compactness calculations
* Plan comparison and classification in `compare_plans()` and `classify_plans()`
* New `iowa` dataset and cleaned-up package data
* New vignettes for redistricting analysis and workflows
* Various bug fixes


# redist 2.0.4

* New `redist.subset` allows for easy subsetting of an adjacency graph
* Added a `NEWS.md` file to track changes to the package
