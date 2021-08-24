# redist 3.1.0

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
