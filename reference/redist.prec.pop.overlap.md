# Compare the Population Overlap Across Plans at the Precinct Level

Compare the Population Overlap Across Plans at the Precinct Level

## Usage

``` r
redist.prec.pop.overlap(
  plan_old,
  plan_new,
  total_pop,
  weighting = "s",
  normalize = TRUE,
  index_only = FALSE,
  return_mat = FALSE
)
```

## Arguments

- plan_old:

  The reference plan to compare against

- plan_new:

  The new plan to compare to the reference plan

- total_pop:

  The total population by precinct This can also take a redist_map
  object and will use the population in that object. If nothing is
  provided, it weights all entries in plan equally.

- weighting:

  Should weighting be done by sum of populations `'s'`, mean of
  populations `'m'`, geometric mean of populations `'g'`, or none `'n'`

- normalize:

  Should entries be normalized by the total population

- index_only:

  Default is FALSE. TRUE returns only one numeric index, the mean of the
  upper triangle of the matrix, under the weighting and normalization
  chosen.

- return_mat:

  Defaults to FALSE, where it returns the summary by row. If TRUE
  returns matrix with length(plan_old) rows and columns. Ignored if
  index_only = TRUE.

## Value

numeric vector with length(plan_old) entries

## Examples

``` r
set.seed(5)
data(iowa)
iowa_map <- redist_map(iowa, total_pop = pop, pop_tol = 0.01, ndists = 4)
plans <- redist_smc(iowa_map, 2, silent = TRUE)
plans_mat <- get_plans_matrix(plans)
ov_vec <- redist.prec.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map)
redist.prec.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map,  weighting = "s",
    normalize = FALSE, index_only = TRUE)
#> [1] 14020.52
```
