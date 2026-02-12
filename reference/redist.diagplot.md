# Diagnostic plotting functionality for MCMC redistricting.

`redist.diagplot` generates several common MCMC diagnostic plots.

## Usage

``` r
redist.diagplot(sumstat,
plot = c("trace", "autocorr", "densplot", "mean", "gelmanrubin"),
logit = FALSE, savename = NULL)
```

## Arguments

- sumstat:

  A vector, list, `mcmc` or `mcmc.list` object containing a summary
  statistic of choice.

- plot:

  The type of diagnostic plot to generate: one of "trace", "autocorr",
  "densplot", "mean", "gelmanrubin". If `plot = "gelmanrubin"`, the
  input `sumstat` must be of class `mcmc.list` or `list`.

- logit:

  Flag for whether to apply the logistic transformation for the summary
  statistic. The default is `FALSE`.

- savename:

  Filename to save the plot. Default is `NULL`.

## Value

Returns a plot of file type `.pdf`.

## Details

This function allows users to generate several standard diagnostic plots
from the MCMC literature, as implemented by Plummer et. al (2006).
Diagnostic plots implemented include trace plots, autocorrelation plots,
density plots, running means, and Gelman-Rubin convergence diagnostics
(Gelman & Rubin 1992).

## References

Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
(2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
Carlo." Working Paper. Available at
<http://imai.princeton.edu/research/files/redist.pdf>.

Gelman, Andrew and Donald Rubin. (1992) "Inference from iterative
simulations using multiple sequences (with discussion)." Statistical
Science.

Plummer, Martin, Nicky Best, Kate Cowles and Karen Vines. (2006) "CODA:
Convergence Diagnosis and Output Analysis for MCMC." R News.

## Examples

``` r
# \donttest{
data(fl25)
data(fl25_enum)
data(fl25_adj)

## Get an initial partition
init_plan <- fl25_enum$plans[, 5118]
fl25$init_plan <- init_plan

## 25 precinct, three districts - no pop constraint ##
fl_map <- redist_map(fl25, existing_plan = 'init_plan', adj = fl25_adj)
#> Projecting to CRS 3857
alg_253 <- redist_flip(fl_map, nsims = 10000)
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> Starting chain 1
#> ℹ Starting swMH().
#> ■                                  0% | ETA:10s
#> ■■■■■■■■■                         26% | ETA:  2s | MH Acceptance: 0.78
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.78
#> 

## Get Republican Dissimilarity Index from simulations
rep_dmi_253 <- redistmetrics::seg_dissim(alg_253, fl25, mccain, pop) |>
    redistmetrics::by_plan(ndists = 3)

## Generate diagnostic plots
redist.diagplot(rep_dmi_253, plot = "trace")

redist.diagplot(rep_dmi_253, plot = "autocorr")

redist.diagplot(rep_dmi_253, plot = "densplot")

redist.diagplot(rep_dmi_253, plot = "mean")


## Gelman Rubin needs two chains, so we run a second
alg_253_2 <- redist_flip(fl_map, nsims = 10000)
#> 
#> ── redist_flip() ───────────────────────────────────────────────────────────────
#> 
#> ── Automated Redistricting Simulation Using Markov Chain Monte Carlo ──
#> ℹ Preprocessing data.
#> Starting chain 1
#> ℹ Starting swMH().
#> ■                                  0% | ETA:11s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s | MH Acceptance: 0.78
#> 

rep_dmi_253_2 <- redistmetrics::seg_dissim(alg_253_2, fl25, mccain, pop) |>
    redistmetrics::by_plan(ndists = 3)

## Make a list out of the objects:
rep_dmi_253_list <- list(rep_dmi_253, rep_dmi_253_2)

## Generate Gelman Rubin diagnostic plot
redist.diagplot(sumstat = rep_dmi_253_list, plot = "gelmanrubin")


# }
```
