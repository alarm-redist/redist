####################################
# Author: Ben Fifield
# Created: 2014/09/21
# Last revision: 2014/09/30
# Institution: Princeton University
# Purpose: Called by runSWA.R after testMHmin.R
#          to combine data
####################################

## Set up container objects
cds <- matrix(NA, nrow = nrow(geodat),
              ncol = (nsims * loop / thin))

ncc <- rep(NA, (nsims * loop / thin))
nacc <- rep(NA, (nsims * loop / thin))
rsplit <- rep(NA, (nsims * loop / thin))
rpar <- rep(NA, (nsims * loop / thin))
radj <- rep(NA, (nsims * loop / thin))
rlam <- rep(NA, (nsims * loop / thin))
beta <- rep(NA, (nsims * loop / thin))
betadiss <- rep(NA, (nsims * loop / thin))
betapop <- rep(NA, (nsims * loop / thin))
betaswitch <- rep(NA, (nsims * loop / thin))
ssdvec <- rep(NA, (nsims * loop / thin))
dissvec <- rep(NA, (nsims * loop / thin))
popvec <- rep(NA, (nsims * loop / thin))
switchvec <- rep(NA, (nsims * loop / thin))
decisionvec <- rep(NA, (nsims * loop / thin))

## Indices for thinning
indthin <- which((1:nsims) %% thin == 0)

#####################################
## Loop over simulation parameters ##
#####################################

for(i in 1:loop){

    ## Load data set
    load(paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
               "_", margin.pct, "_", lambda,
               "_b", params$initbeta[aid] * -1, "a", abeta,
               "_bDiss", params$initbetadiss[aid], "a", abetadiss,
               "_bPop", params$initbetapop[aid] * -1, "a", abetapop,
               "t", params$targbetapop[aid] * -1, "b", params$bybetapop[aid] * -1,
               "_bSwitch", params$initbetaswitch[aid] * -1, "a", abetaswitch,
               "_pow", wpow,
               "_par", pnum, "_loop", i, ".RData", sep = ""))

    ind <- ((i - 1) * (nsims / thin) + 1):(i * (nsims / thin))
    
    ## Store objects together
    cds[1:nrow(geodat), ind] <- ecuts[[1]][,indthin]

    ncc[ind] <- ecuts[[2]][indthin]
    nacc[ind] <- ecuts[[3]][indthin]
    rsplit[ind] <- ecuts[[4]][indthin]
    rpar[ind] <- ecuts[[5]][indthin]
    radj[ind] <- ecuts[[6]][indthin]
    rlam[ind] <- ecuts[[7]][indthin]
    beta[ind] <- ecuts[[8]][indthin]
    betadiss[ind] <- ecuts[[9]][indthin]
    betapop[ind] <- ecuts[[10]][indthin]
    betaswitch[ind] <- ecuts[[11]][indthin]
    ssdvec[ind] <- ecuts[[12]][indthin]
    dissvec[ind] <- ecuts[[13]][indthin]
    popvec[ind] <- ecuts[[14]][indthin]
    switchvec[ind] <- ecuts[[15]][indthin]
    decisionvec[ind] <- ecuts[[16]][indthin]

}

## Store data in edgecuts object
edgecuts <- vector(mode = "list", length = 16)
edgecuts[[1]] <- cds
edgecuts[[2]] <- ncc
edgecuts[[3]] <- nacc
edgecuts[[4]] <- rsplit
edgecuts[[5]] <- rpar
edgecuts[[6]] <- radj
edgecuts[[7]] <- rlam
edgecuts[[8]] <- beta
edgecuts[[9]] <- betadiss
edgecuts[[10]] <- betapop
edgecuts[[11]] <- betaswitch
edgecuts[[12]] <- ssdvec
edgecuts[[13]] <- dissvec
edgecuts[[14]] <- popvec
edgecuts[[15]] <- switchvec
edgecuts[[16]] <- decisionvec

## Save full object
save(edgecuts, file = paste(dwd, "edgecuts", state, "_", (1 - eprob) * 100,
               "_", margin.pct, "_", lambda,
               "_b", params$initbeta[aid] * -1, "a", abeta,
               "_bDiss", params$initbetadiss[aid], "a", abetadiss,
               "_bPop", params$initbetapop[aid] * -1, "a", abetapop,
               "t", params$targbetapop[aid] * -1, "b", params$bybetapop[aid] * -1,
               "_bSwitch", params$initbetaswitch[aid] * -1, "a", abetaswitch,
               "_pow", wpow,
               "_par", pnum, ".RData", sep = ""))

