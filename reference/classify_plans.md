# Hierarchically classify a set of redistricting plans

Applies hierarchical clustering to a distance matrix computed from a set
of plans and takes the first `k` splits.

## Usage

``` r
classify_plans(dist_mat, k = 8, method = "complete")
```

## Arguments

- dist_mat:

  a distance matrix, the output of
  [`plan_distances()`](http://alarm-redist.org/redist/reference/redist.distances.md)

- k:

  the number of groupings to create

- method:

  the clustering method to use. See
  [`hclust()`](https://rdrr.io/r/stats/hclust.html) for options.

## Value

An object of class `redist_classified`, which is a list with two
elements:

- groups:

  A character vector of group labels of the form `"I.A.1.a.i"`, one for
  each plan.

- splits:

  A list of splits in the hierarchical clustering. Each list element is
  a list of two mutually exclusive vectors of plan indices, labeled by
  their group classification, indicating the plans on each side of the
  split.

Use
[`plot.redist_classified()`](http://alarm-redist.org/redist/reference/plot.redist_classified.md)
for a visual summary.
