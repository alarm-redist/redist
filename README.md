
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **redist**: Simulation Methods for Legislative Redistricting

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/kosukeimai/redist.svg?branch=master)](https://travis-ci.org/kosukeimai/redist)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version-last-release/redist)](https://cran.r-project.org/package=redist)
![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/redist)

<!-- badges: end -->

<img src="man/figures/map_photo.jpg" style="width: 100%"/>

<img src="man/figures/logo.png" align="right" height="128" style="margin-left: 4px;"/>

This R package enables researchers to sample redistricting plans from a
pre-specified target distribution using Sequential Monte Carlo and
Markov Chain Monte Carlo algorithms. The package supports various
constraints in the redistricting process, such as geographic compactness
and population parity requirements. Tools for analysis, including
computation of various summary statistics and plotting functionality,
are also included.

Authors:

-   [Christopher T Kenny](https://www.christophertkenny.com),
    <christopherkenny@fas.harvard.edu> (Maintainer)
-   [Cory McCartan](https://corymccartan.github.io),
    <cmccartan@g.harvard.edu>
-   [Ben Fifield](https://benfifield.com), <benfifield@gmail.com>
-   [Kosuke Imai](https://imai.fas.harvard.edu), <imai@harvard.edu>

Contributors:

-   Jun Kawahara, <jkawahara@i.kyoto-u.ac.jp>
-   Alex Tarr, <atarr@princeton.edu>
-   [Michael Higgins](http://www-personal.k-state.edu/~mikehiggins/),
    <mjh5@princeton.edu>

Papers:

-   [Automated Redistricting Simulation Using Markov Chain Monte
    Carlo](https://doi.org/10.1080/10618600.2020.1739532) *Journal of
    Computational and Graphical Statistics*
-   [The Essential Role of Empirical Validation in Legislative
    Redistricting
    Simulation](https://doi.org/10.1080/2330443X.2020.1791773)
    *Statistics and Public Policy* Vol. 7, No. 1, pp. 52-68.
-   [Sequential Monte Carlo for Sampling Balanced and Compact
    Redistricting Plans](https://arxiv.org/pdf/2008.06131.pdf)

## Installation Instructions

`redist` is available on CRAN and can be installed using:

``` r
install.packages("redist")
```

You can also install the most recent development version of `redist`
using the \`remotes\`\` package.

``` r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("alarm-redist/redist", dependencies=TRUE)
```

## Getting started

A basic analysis has two steps. First, you define a redistricting plan
using `redist_map`. Then you simulate plans using one of the algorithm
functions: `redist_smc`, `redist_flip`, and `redist_mergesplit`.

``` r
library(redist)
library(dplyr)

data(iowa)

# set a 0.01% population constraint
iowa_map = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.0001, total_pop = pop)
# simulate 250 plans using the SMC algorithm
iowa_plans = redist_smc(iowa_map, nsims=250, verbose=FALSE)
#> SEQUENTIAL MONTE CARLO
#> Sampling 250 99-unit maps with 4 districts and population between 761513 and 761665.
#> Making split 1 of 3
#> Note: maximum hit; falling back to naive k estimator.
#> Resampling effective sample size: 246.1 (98.4% efficiency).
#> Making split 2 of 3
#> Resampling effective sample size: 246.2 (98.5% efficiency).
#> Making split 3 of 3
#> Resampling effective sample size: 244.3 (97.7% efficiency).
```

After generating plans, you can use `redist`’s plotting functions to
study the geographic and partisan characteristics of the simulated
ensemble.

``` r
library(ggplot2)
library(patchwork) # for plotting

redist.plot.plans(iowa_plans, draws=c("cd_2010", "1", "2", "3"),
                  geom=iowa_map)
```

![](man/figures/README-readme-plot-1.png)<!-- -->

``` r
dev_comp = iowa_plans %>%
    mutate(comp = distr_compactness(iowa_map)) %>%
    group_by(draw) %>%
    summarize(`Population deviation` = max(abs(total_pop/get_target(iowa_map) - 1)),
              Compactness = comp[1])

hist(dev_comp, `Population deviation`) + hist(dev_comp, Compactness) +
    plot_layout(guides="collect") +
    plot_annotation(title="Simulated plan characteristics")
```

![](man/figures/README-readme-plot-2.png)<!-- -->

``` r
redist.plot.scatter(dev_comp, `Population deviation`, Compactness) +
    labs(title="Population deviation and compactness by plan")
```

![](man/figures/README-readme-plot-3.png)<!-- -->

``` r
iowa_plans %>%
    mutate(`Democratic vote` = group_frac(iowa_map, dem_08, tot_08)) %>%
    plot(`Democratic vote`, size=0.5, color_thresh=0.5) +
    scale_color_manual(values=c("tomato2", "dodgerblue")) +
    labs(title="Democratic vote share by district")
```

![](man/figures/README-readme-plot-4.png)<!-- -->

A more detailed introduction to redistricting methods and the package
can be found in the [Get
Started](https://alarm-redist.github.io/redist/articles/redist.html)
page. The package
[vignettes](https://alarm-redist.github.io/redist/articles/) contain
more detailed information and guides to specific workflows.

## Recommended Workflows

In some cases, there are two functions with the same purpose that
operate on different inputs. For most users, we recommend working with
the versions that operate on `redist_map` and `redist_plan` objects.

| Recommended         | Alternate Options   | Notes                           |
|---------------------|---------------------|---------------------------------|
| `redist_smc`        | `redist.smc`        |                                 |
| `redist_mergesplit` | `redist.mergesplit` |                                 |
| `redist_flip`       | `redist.flip`       | (Formerly called `redist.mcmc`) |
