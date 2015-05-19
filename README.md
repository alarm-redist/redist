# redist 1.1
R package for simulating redistricting plans via Markov chain Monte
Carlo by Ben Fifield ([bfifield@princeton.edu](bfifield@princeton.edu)),
Alex Tarr ([atarr@princeton.edu](atarr@princeton.edu)), Michael
Higgins ([mjh5@princeton.edu](mjh5@princeton.edu)), and Kosuke Imai
([kimai@princeton.edu](kimai@princeton.edu)). Maintainer is Ben Fifield.

## Installation Instructions
The package is available on CRAN and can be installed using:

```
install.packages("redist")
```

Users can also install the most stable development release of the `redist` package using the `install_github()` function in the `devtools` package.

```
library(devtools)
install_devtools("redistricting/redist")
```

## src Folder
We hope the following guide will be of help to users who want to take a look at the original
`redist` source code:
- `sw_mh_alg.cpp`: Contains the `swMH()` function, which conducts
Markov chain Monte Carlo simulation of redistricting plans.
- `sw_mh_helper.cpp`: A series of functions to aid in simulating
  redistricting plans.
- `make_swaps_helper.cpp`: A series of functions to propose and make
swaps of geographic units in the primary redistricting algorithm.
- `constraint_calc_helper.cpp`: Functions to calculate the strength of
certain implemented constraints such as population and compactness
requirements.
- `rsg.cpp`: An implementation of the random seed-and-grow algorithm
described in detail in Chen and Rodden (2013).
- `check_contiguity.cpp`: A contiguity check for the implementation of
the Chen and Rodden (2013) algorithm in `rsg.cpp`.
- `redist_analysis.cpp`: Functions to aid in analysis of simulated
redistricting plans.
  

