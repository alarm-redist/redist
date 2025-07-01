# 4.3.0
* Improves SMC performance by pre-allocating some memory while drawing spanning trees.
* Replaces SMC label-counting adjustments (exact and importance-sampling-based) with a new backward kernel that eliminates approximation error and requires far less computation
* 4.2.0 introduced some regressions in `redist_shortburst()` along with the new features. The following issues are fixed: 
  * the returned `redist_plans` object would store the wrong score for the '<init>' input. This issue only impacted the returned object and the correct score was used during the shortburst process. (#180)
  * the function would return immediately if `stop_at` was specified and `minimize = FALSE`. (#181) 
  
 * Add `summary()` support for plans sampled with the `flip` algorithm. This does not replace the full flip diagnostic suite, but provides an easy way to compute r-hats.

# 4.2.0
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
