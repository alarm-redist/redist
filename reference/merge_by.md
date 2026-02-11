# Merge map units

In performing a county-level or cores-based analysis it is often
necessary to merge several units together into a larger unit. This
function performs this operation, modifying the adjacency graph as
needed and attempting to properly aggregate other data columns.

## Usage

``` r
merge_by(.data, ..., by_existing = TRUE, drop_geom = TRUE, collapse_chr = TRUE)
```

## Arguments

- .data:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- ...:

  [`<tidy-select>`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)
  the column(s) to merge by

- by_existing:

  if an existing assignment is present, whether to also group by it

- drop_geom:

  whether to drop the geometry column. Recommended, as otherwise a
  costly geometric merge is required.

- collapse_chr:

  if `TRUE`, preserve character columns by collapsing their values. For
  example, a county name column in Iowa might be merged and have entries
  such as "Cedar~Clinton~Des Moines". Set to `FALSE` to drop character
  columns instead.

## Value

A merged
[`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
object
