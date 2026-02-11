# Calculates Sparse Population Moves to Minimize Population Deviation

This function computes a minimal set of population moves (e.g., 5 people
from district 1 to district 3) to maximally balance the population
between districts. The moves are only allowed between districts that
share the territory of a county, so that any boundary adjustments are
guaranteed to preserve all unbroken county boundaries.

## Usage

``` r
min_move_parity(map, plan, counties = NULL, penalty = 0.2)
```

## Arguments

- map:

  a [redist_map](http://alarm-redist.org/redist/reference/redist_map.md)

- plan:

  an integer vector containing the plan to be balanced. Tidy-evaluated.

- counties:

  an optional vector of counties, whose boundaries will be preserved.
  Tidy-evaluated.

- penalty:

  the larger this value, the more to encourage sparsity.

## Value

a list with components:

- `moves`:

  A tibble describing the population moves

- `pop_old`:

  The current district populations

- `pop_new`:

  The district populations after the moves

## Examples

``` r
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
min_move_parity(iowa_map, cd_2010)
#> $moves
#> # A tibble: 3 × 3
#>    from    to  move
#>   <int> <int> <dbl>
#> 1     2     1    35
#> 2     4     1     5
#> 3     3     4    23
#> 
#> $pop_old
#> [1] 761548 761624 761612 761571
#> 
#> $pop_new
#> [1] 761588 761589 761589 761589
#> 
```
