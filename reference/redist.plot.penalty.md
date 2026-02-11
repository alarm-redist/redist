# (Deprecated) Visualize Group Power Penalty

Plots the shape of the
[`add_constr_grp_pow()`](http://alarm-redist.org/redist/reference/constraints.md)
penalty.

## Usage

``` r
redist.plot.penalty(
  tgt_min = 0.55,
  tgt_other = 0.25,
  strength_vra = 2500,
  pow_vra = 1.5,
  limits = TRUE
)
```

## Arguments

- tgt_min:

  double, defaults to 0.55. The minority target percent.

- tgt_other:

  double, defaults to 0.25. The other group target percent.

- strength_vra:

  double, strength of the VRA constraint.

- pow_vra:

  double, exponent of the VRA constraint.

- limits:

  Whether to limit y axis to 0,500. Default is TRUE for comparability
  across values.

## Value

ggplot

## Details

This function allows you to plot the un-exponentiated penalty
implemented as
[`add_constr_grp_pow()`](http://alarm-redist.org/redist/reference/constraints.md).
The function takes two key inputs, `tgt_min` and `tgt_other` which
center the minimum penalty spots. A higher y-value indicates a higher
penalty and incentivizes moving towards a spot with a lower y-value. The
x-axis indicates the group population proportion in a given district.
