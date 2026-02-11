# Renumber districts to match an existing plan

District numbers in simulated plans are by and large random. This
function attempts to renumber the districts across all simulated plans
to match the numbers in a provided plan, using the Hungarian algorithm.

## Usage

``` r
match_numbers(
  data,
  plan,
  total_pop = attr(data, "prec_pop"),
  col = "pop_overlap"
)
```

## Arguments

- data:

  a `redist_plans` object.

- plan:

  a character vector giving the name of the plan to match to (e.g., for
  a reference plan), or an integer vector containing the plan itself.

- total_pop:

  a vector of population counts. Should not be needed for most
  `redist_plans` objects.

- col:

  the name of a new column to store the vector of population overlap
  with the reference plan: the fraction of the total population who are
  in the same district under each plan and the reference plan. Set to
  `NULL` if no column should be created. renumbering options in any
  plan.

## Value

a modified `redist_plans` object. New district numbers will be stored as
an ordered factor variable in the `district` column. The district
numbers in the plan matrix will match the levels of this factor.

## Examples

``` r
data(iowa)

iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
plans <- redist_smc(iowa_map, 100, silent = TRUE)
match_numbers(plans, "cd_2010")
#> A <redist_plans> containing 100 sampled plans and 1 reference plan
#> Plans have 4 districts from a 99-unit map, and were drawn using Sequential
#> Monte Carlo.
#> With plans resampled from weights
#> Plans matrix: int [1:99, 1:101] 1 1 2 3 4 2 2 4 2 2 ...
#> # A tibble: 404 × 4
#>    draw    district total_pop pop_overlap
#>    <fct>   <ord>        <dbl>       <dbl>
#>  1 cd_2010 1           761612       1    
#>  2 cd_2010 2           761548       1    
#>  3 cd_2010 3           761624       1    
#>  4 cd_2010 4           761571       1    
#>  5 1       1           758881       0.657
#>  6 1       2           741379       0.657
#>  7 1       3           796636       0.657
#>  8 1       4           749459       0.657
#>  9 2       1           790031       0.824
#> 10 2       2           724485       0.824
#> # ℹ 394 more rows
```
