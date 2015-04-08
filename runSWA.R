####################################
# Author: Ben Fifield
# Created: 2014/09/20
# Last Revision: 2014/11/09
# Institution: Princeton University
# Purpose: Sets parameters, calls script to run SWA,
#          calls script to combine data, and
#          calls script to generate diagnostics
####################################
rm(list = ls())

## Set parameters
a <- "Mac"
params <- expand.grid(state = "testset252",
                      eprob = 0.05, marginpct = 5,
                      lambda = 1, pnum = 1:10,
                      initbeta = 0,
                      initbetadiss = 0,
                      initbetapop = 0,
                      initbetaswitch = 0,
                      annealbeta = 0, annealbetadiss = 0,
                      annealbetapop = 1, annealbetaswitch = 0,
                      targbetapop = -10,
                      bybetapop = -1,
                      weightpow = 1,
                      nsims = 50, loop = 6,
                      thin = 10)

## Set working directory based on environmenten
if(a == "Tukey"){
    codedir <- "/scratch/network/bfifield/segregation/code/"
}
if(a == "Mac"){
    codedir <- "~/Dropbox/Graduate School/Projects/segregation/svn/code/"
}

## Run simulations
source(paste(codedir, "testMHmin.R", sep = ""))

## Combine data from loops
source(paste(codedir, "combineDataSplit.R", sep = ""))

## Generate diagnostics
source(paste(codedir, "genDiag.R", sep = ""))

