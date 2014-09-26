####################################
# Author: Ben Fifield
# Created: 2014/09/20
# Last Revision: 2014/09/22
# Institution: Princeton University
# Purpose: Sets parameters, calls script to run SWA,
#          calls script to combine data, and
#          calls script to generate diagnostics
####################################
rm(list = ls())

## Set parameters
a <- "Tukey"
params <- expand.grid(state = "OK",
                      eprob = 0.05, marginpct = 1,
                      lambda = c(15, 20, 25), pnum = 1:10,
                      initbeta = 0,
                      initbetadiss = 0,
                      initbetapop = 0,
                      initbetaswitch = 0,
                      annealbeta = 0, annealbetadiss = 0,
                      annealbetapop = 1, annealbetaswitch = 0,
                      targbetapop = -150,
                      bybetapop = c(-25, -30),
                      weightpow = 14,
                      nsims = 25000, loop = 4)

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

