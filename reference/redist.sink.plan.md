# Sink Plans to 1:ndists

Takes a plan and renumbers it to be from 1:ndists

## Usage

``` r
redist.sink.plan(plan)
```

## Arguments

- plan:

  vector of assignments, required.

## Value

A vector with an ID that corresponds from 1:ndists, and attribute `n`
indicating the number of districts.

## Examples

``` r
data(fl25_enum)
plan <- fl25_enum$plans[, 5118]
# Subset based on something:
plan <- plan[plan != 2]
plan <- vctrs::vec_group_id(plan)
# Now plan can be used with redist_flip()
plan
#>  [1] 1 2 1 2 2 1 2 2 2 2 2 2
#> attr(,"n")
#> [1] 2
```
