# Calculate gerrymandering metrics for a set of plans

`redist.metrics` is used to compute different gerrymandering metrics for
a set of maps.

## Usage

``` r
partisan_metrics(map, measure, rvote, dvote, ..., .data = cur_plans())

redist.metrics(
  plans,
  measure = "DSeats",
  rvote,
  dvote,
  tau = 1,
  biasV = 0.5,
  respV = 0.5,
  bandwidth = 0.01,
  draw = 1
)
```

## Arguments

- map:

  a
  [`redist_map`](http://alarm-redist.org/redist/reference/redist_map.md)
  object

- measure:

  A vector with a string for each measure desired from list "DSeats",
  "DVS", "EffGap", "EffGapEqPop", "TauGap", "MeanMedian", "Bias",
  "BiasV", "Declination", "Responsiveness", "LopsidedWins",
  "RankedMarginal", and "SmoothedSeat". Use "all" to get all metrics.
  "DSeats" and "DVS" are always computed, so it is recommended to always
  return those values.

- rvote:

  A numeric vector with the Republican vote for each precinct.

- dvote:

  A numeric vector with the Democratic vote for each precinct.

- ...:

  passed on to `redist.metrics`

- .data:

  a
  [`redist_plans`](http://alarm-redist.org/redist/reference/redist_plans.md)
  object

- plans:

  A numeric vector (if only one map) or matrix with one row for each
  precinct and one column for each map. Required.

- tau:

  A non-negative number for calculating Tau Gap. Only used with option
  "TauGap". Defaults to 1.

- biasV:

  A value between 0 and 1 to compute bias at. Only used with option
  "BiasV". Defaults to 0.5.

- respV:

  A value between 0 and 1 to compute responsiveness at. Only used with
  option "Responsiveness". Defaults to 0.5.

- bandwidth:

  A value between 0 and 1 for computing responsiveness. Only used with
  option "Responsiveness." Defaults to 0.01.

- draw:

  A numeric to specify draw number. Defaults to 1 if only one map
  provided and the column number if multiple maps given. Can also take a
  factor input, which will become the draw column in the output if its
  length matches the number of entries in plans. If the `plans` input is
  a `redist_plans` object, it extracts the `draw` identifier.

## Value

A tibble with a column for each specified measure and a column that
specifies the map number.

## Details

This function computes specified compactness scores for a map. If there
is more than one precinct specified for a map, it aggregates to the
district level and computes one score.

- DSeats is computed as the expected number of Democratic seats with no
  change in votes.

- DVS is the Democratic Vote Share, which is the two party vote share
  with Democratic votes as the numerator.

- EffGap is the Efficiency Gap, calculated with votes directly.

- EffGapEqPop is the Efficiency Gap under an Equal Population
  assumption, calculated with the DVS.

- TauGap is the Tau Gap, computed with the Equal Population assumption.

- MeanMedian is the Mean Median difference.

- Bias is the Partisan Bias computed at 0.5.

- BiasV is the Partisan Bias computed at value V.

- Declination is the value of declination at 0.5.

- Responsiveness is the responsiveness at the user-supplied value with
  the user-supplied bandwidth.

- LopsidedWins computed the Lopsided Outcomes value, but does not
  produce a test statistic.

- RankedMarginal computes the Ranked Marginal Deviation (0-1, smaller is
  better). This is also known as the "Gerrymandering Index" and is
  sometimes presented as this value divided by 10000.

- SmoothedSeat computes the Smoothed Seat Count Deviation (0-1, smaller
  is R Bias, bigger is D Bias).

## References

Jonathan N. Katz, Gary King, and Elizabeth Rosenblatt. 2020. Theoretical
Foundations and Empirical Evaluations of Partisan Fairness in
District-Based Democracies. American Political Science Review, 114, 1,
Pp. 164-178.

Gregory S. Warrington. 2018. "Quantifying Gerrymandering Using the Vote
Distribution." Election Law Journal: Rules, Politics, and Policy. Pp.
39-57.http://doi.org/10.1089/elj.2017.0447

Samuel S.-H. Wang. 2016. "Three Tests for Practical Evaluation of
Partisan Gerrymandering." Stanford Law Review, 68, Pp. 1263 - 1321.

Gregory Herschlag, Han Sung Kang, Justin Luo, Christy Vaughn Graves,
Sachet Bangia, Robert Ravier & Jonathan C. Mattingly (2020) Quantifying
Gerrymandering in North Carolina, Statistics and Public Policy, 7:1,
30-38, DOI: 10.1080/2330443X.2020.1796400

## Examples

``` r
data(fl25)
data(fl25_enum)
plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
# old: redist.metrics(plans_05, measure = "DSeats", rvote = fl25$mccain, dvote = fl25$obama)
part_dseats(plans_05, fl25, mccain, obama)
#>   [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 3 3 3 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 2
#>  [38] 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2
#>  [75] 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3
#> [112] 2 2 2 2 2 2 3 3 3 3 3 3 2 2 2 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2
#> [149] 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2
#> [186] 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2
#> [223] 3 3 3 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [260] 2 2 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 3
#> [297] 3 3 3 3 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2
#> [334] 2 2 2 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 3 3 3 2 2 2 2
#> [371] 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 3 3
#> [408] 3 3 3 3 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3
#> [445] 3 3 3 2 2 2 3 3 3 2 2 2 2 2 2 3 3 3 2 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 3
#> [482] 3 3 2 2 2 3 3 3 3 3 3 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3 2 2 2 2 2
#> [519] 2 3 3 3 3 3 3 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3
#> [556] 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3
```
