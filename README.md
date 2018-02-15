redist: Markov Chain Monte Carlo Methods for Redistricting Simulation [![Build Status](https://travis-ci.org/kosukeimai/redist.svg?branch=master)](https://travis-ci.org/kosukeimai/redist) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/redist)](https://cran.r-project.org/package=redist)
===========================================================================================================================================================================================================================================================================================================

Authors:

-   [Ben Fifield](https://www.benfifield.com), <bfifield@princeton.edu> (Maintainer)
-   Alex Tarr, <atarr@princeton.edu>
-   [Michael Higgins](http://www-personal.k-state.edu/~mikehiggins/), <mjh5@princeton.edu>
-   [Kosuke Imai](https://imai.princeton.edu), <kimai@princeton.edu>

Paper: [A New Automated Redistricting Simulator Using Markov Chain Monte Carlo](http://imai.princeton.edu/research/files/redist.pdf)

Installation Instructions
-------------------------

`redist` is available on CRAN and can be installed using:

``` r
install.packages("redist")
```

You can also install the most recent development version of `redist` using the `devtools` package. First you have to install `devtools` using the following code. Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install `redist`:

``` r
library(devtools)
install_github("kosukeimai/redist",dependencies=TRUE)
```

Usage Examples - 25 Precincts into 3 Districts
==============================================

No Population Constraint
------------------------

``` r
## Load data
library(redist)
data(algdat.pfull)

## Run the simulations
mcmc.out <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                        popvec = algdat.pfull$precinct.data$pop,
                        nsims = 10000,
                        ndists = 3)
```

20% Population Constraint
-------------------------

``` r
## Load data
data(algdat.p20)

## -----------------------------------------------
## Run mcmc algorithm - hard population constraint
## (reject any sample where population parity > 20%)
## -----------------------------------------------
mcmc.out <- redist.mcmc(adjobj = algdat.p20$adjlist,
                        popvec = algdat.p20$precinct.data$pop,
                        nsims = 10000,
                        popcons = .2,
                        ndists = 3)

## ---------------------------------------------------------------------
## Run mcmc algorithm - draws from gibbs defined by distance from parity
## (Run with no tempering)
## ---------------------------------------------------------------------
mcmc.out.gb <- redist.mcmc(adjobj = algdat.p20$adjlist,
                           popvec = algdat.p20$precinct.data$pop,
                           ndists = 3,
                           nsims = 10000,
                           constraint = "population",
                           constraintweights = 5.4)
## Reweight draws back to the uniform distribution
mcmc.out.gb <- redist.ipw(mcmc.out.gb, targetpop = .2)

## ---------------------------------------------------------------------
## Run mcmc algorithm - draws from gibbs defined by distance from parity
## (Run with simulated tempering, betas power law sequence of length 10
## from 0 to 1)
## Also including optional beta weights, to upweight prob of transitions
## to colder temperatures
## ---------------------------------------------------------------------
betaweights <- rep(NA, 10); for(i in 1:10){betaweights[i] <- 2^i}
mcmc.out.st <- redist.mcmc(adjobj = algdat.p20$adjlist,
                           popvec = algdat.p20$precinct.data$pop,
                           ndists = 3,
                           nsims = 10000,
                           constraint = "population",
                           constraintweights = 5.4,
                           temper = TRUE,
                           betaweights = betaweights)
mcmc.out.st <- redist.ipw(mcmc.out.st, targetpop = .2)
```

20% Population Constraint and Compactness Constraint
----------------------------------------------------

``` r
## ----------------------------------------------------------
## Constrain on population and compactness with tempering,
## weight on population = 5.4 while weight on compactness = 3
## Also specifying argument for ssdmat, the distance matrix
## ----------------------------------------------------------
mcmc.out.st.multiple <- redist.mcmc(adjobj = algdat.p20$adjlist,
                                    popvec = algdat.p20$precinct.data$pop,
                                    ndists = 3,
                                    nsims = 10000,
                                    constraint = c("population", "compact"),
                                    constraintweights = c(5.4, 3),
                                    ssdmat = algdat.p20$distancemat,
                                    temper = TRUE,
                                    betaweights = betaweights)
mcmc.out.st <- redist.ipw(mcmc.out.st.multiple, targetpop = .2)
```

20% Population Constraint and Compactness Constraint, using MPI
---------------------------------------------------------------

``` r
## ----------------------------------------------------------
## Constrain on population and compactness with parallel tempering,
## weight on population = 5.4 while weight on compactness = 3
## Also specifying argument for ssdmat, the distance matrix.
## In addition, specifying a 20-beta tempering ladder.
## Save file as "redist_mpi.RData"
## ----------------------------------------------------------
redist.mcmc.mpi(
  adjobj = algdat.p20$adjlist,
  popvec = algdat.p20$precinct.data$pop,
  ndists = 3,
  nsims = 10000,
  constraint = c("population", "compact"),
  constraintweights = c(5.4, 3),
  ssdmat = algdat.p20$distancemat,
  betaseqlength = 20,
  savename = "redist_mpi",
  verbose = TRUE
)

## Note that reweighting using redist.ipw() currently
## has to happen in a separate analysis file after
## redist.mcmc.mpi() is run and the output file is saved. 
## We will change this in a future release to run the 
## inverse probability reweighting inside of
## redist.mcmc.mpi().
```

### A Sample slurm Script for Submitting MPI Jobs

``` bash
## ----------------------------------------------------------
## Sample slurm script for submitting an R script using 
## redist.mcmc.mpi(). Set ntasks equal to the length of the
## beta sequence ladder plus one extra core.
## ----------------------------------------------------------

#!/usr/bin/env bash
#SBATCH --ntasks=21
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -J [jobname here]
#SBATCH -o log.%j
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------------'
echo -e '\n\n'

mpiexec -np 1 R --no-save < [scriptname here].R
```
