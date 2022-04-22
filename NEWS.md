# redist 4.0.0
* A new constraint interface that is more flexible and user friendly (see `redist_constr()` and `?constraints`)
* New diagnostic-checking function, `summary.redist_plans()`
* Many bug fixes (see https://github.com/alarm-redist/redist/issues)
* Remove deprecated functions
* Update the SMC algorithm to include a missing correction factor for the number
of ways to sequentially label districts. This factor should not have an effect
on substantive conclusions and summary statistics.

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
