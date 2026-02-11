# Plot a histogram of a summary statistic

Plots a histogram of a statistic of a
[`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
object, with a reference line for each reference plan, if applicable.

## Usage

``` r
redist.plot.hist(plans, qty, bins = NULL, ...)

# S3 method for class 'redist_plans'
hist(x, qty, ...)
```

## Arguments

- plans:

  the `redist_plans` object.

- qty:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the statistic.

- bins:

  the number of bins to use in the histogram. Defaults to
  Freedman-Diaconis rule.

- ...:

  passed on to
  [`geom_histogram`](https://ggplot2.tidyverse.org/reference/geom_histogram.html)

- x:

  [`<data-masking>`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)
  the statistic.

## Value

A ggplot

## Examples

``` r
library(dplyr)
data(iowa)

iowa <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
plans <- redist_smc(iowa, nsims = 100, silent = TRUE)
group_by(plans, draw) %>%
    summarize(pop_dev = max(abs(total_pop/mean(total_pop) - 1))) %>%
    redist.plot.hist(pop_dev)

```
