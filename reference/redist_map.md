# Create a `redist_map` object.

Sets up a redistricting problem.

## Usage

``` r
redist_map(
  ...,
  existing_plan = NULL,
  pop_tol = NULL,
  total_pop = c("pop", "population", "total_pop", "POP100"),
  ndists = NULL,
  pop_bounds = NULL,
  adj = NULL,
  adj_col = "adj",
  planarize = 3857
)

as_redist_map(x)
```

## Arguments

- ...:

  column elements to be bound into a `redist_map` object or a single
  `list` or `data.frame`. These will be passed on to the
  [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
  constructor.

- existing_plan:

  [`<tidy-select>`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)
  the existing district assignment. Must be numeric or convertible to
  numeric.

- pop_tol:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the population tolerance. The percentage deviation from the average
  population will be constrained to be no more than this number. If
  `existing_plan` is provided, defaults to the parity of that plan;
  otherwise, defaults to 0.01.

- total_pop:

  [`<tidy-select>`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)
  the vector of precinct populations. Defaults to the `pop`,
  `population`, or `total_pop` columns, if one exists.

- ndists:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the integer number of districts to partition the map into. Must be
  specified if `existing_plan` is not supplied.

- pop_bounds:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  more specific population bounds, in the form of
  `c(lower, target, upper)`.

- adj:

  the adjacency graph for the object. Defaults to being computed from
  the data if it is coercible to a shapefile.

- adj_col:

  the name of the adjacency graph column

- planarize:

  a number, indicating the CRS to project the shapefile to if it is
  latitude-longitude based. Set to NULL or FALSE to avoid planarizing.

- x:

  an object to be coerced

## Value

A redist_map object

## Details

A `redist_map` object is a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
which contains an adjacency list and additional information about the
number of districts and population bounds. It supports all of the
`dplyr` generics, and will adjust the adjacency list and attributes
according to these functions; i.e., if we `filter` to a subset of units,
the graph will change to subset to these units, and the population
bounds will adjust accordingly. If an existing map is also attached to
the object, the number of districts will also adjust. Subsetting with
`` `[` `` and `` `[[` `` does not recompute graphs or attributes.

Other useful methods for `redist_map` objects:

- [`merge_by`](http://alarm-redist.org/redist/reference/merge_by.md)

- [`get_adj`](http://alarm-redist.org/redist/reference/get_adj.md)

- [`plot.redist_map`](http://alarm-redist.org/redist/reference/plot.redist_map.md)

## Examples

``` r
data(fl25)
d <- redist_map(fl25, ndists = 3, pop_tol = 0.05, total_pop = pop)
#> Projecting to CRS 3857
dplyr::filter(d, pop >= 10e3)
#> Warning: Your subset was not based on districts.
#> → Please use `set_pop_tol()` to update your <redist_map> or create a new
#>   <redist_map> with the correct number of districts.
#> A <redist_map> with 5 units and 13 fields
#> To be partitioned into 3 districts with population between 58,347.67 - 5.0% and 58,347.67 + 5.0%
#> With geometry:
#>     bbox:           xmin: -9107951 ymin: 3047341 xmax: -9065127 ymax: 3085829
#>     projected CRS:  WGS 84 / Pseudo-Mercator
#> # A tibble: 5 × 13
#>   geoid10    pop   vap obama mccain TotPop BlackPop HispPop   VAP BlackVAP
#> * <chr>    <dbl> <dbl> <dbl>  <dbl>  <dbl>    <dbl>   <dbl> <dbl>    <dbl>
#> 1 2519.0_0 15993 12379   647   1522  15993      863    1980 12379      582
#> 2 2396.0_0 12853  8965   640   1027  12853     2368    3860  8965     1468
#> 3 2515.0_0 22218 14513   506    641  22218     5911    8800 14513     3417
#> 4 2525.0_0 16734 14720   826   1414  16734      197     791 14720      155
#> 5 1246.0_0 12289  7734   366    271  12289     1693    9773  7734     1062
#> # ℹ 3 more variables: HispVAP <dbl>, geometry <POLYGON [m]>, adj <list>
```
