# Renumber districts to match a quantity of interest

District numbers in simulated plans are by and large random. This
function will renumber the districts across all simulated plans in order
of a provided quantity of interest.

## Usage

``` r
number_by(data, x, desc = FALSE)
```

## Arguments

- data:

  a `redist_plans` object

- x:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the quantity of interest.

- desc:

  `TRUE` if district should be sorted in descending order.

## Value

a modified `redist_plans` object. New district numbers will be stored as
an ordered factor variable in the `district` column. The district
numbers in the plan matrix will match the levels of this factor.
