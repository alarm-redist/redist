# Calculate Frontier Size

Calculate Frontier Size

## Usage

``` r
redist.calc.frontier.size(ordered_path)
```

## Arguments

- ordered_path:

  path to ordered path created by redist.prep.enumpart

## Value

List, four objects

- `max` numeric, maximum frontier size

- `average` numeric, average frontier size

- `average_sq` numeric, average((frontier size)^2)

- `sequence` numeric vector, lists out all sizes for every frontier

## Examples

``` r
if (FALSE) { # \dontrun{
data(fl25)
adj <- redist.adjacency(fl25)
redist.prep.enumpart(adj, "unordered", "ordered")
redist.calc.frontier.size("ordered")
} # }
```
