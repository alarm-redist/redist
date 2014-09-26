####################################
# Author: Ben Fifield
# Created: 2014/09/21
# Last revision: 2014/09/21
# Institution: Princeton University
# Purpose: Called by runSWA.R after testMHmin.R
#          to combine data
####################################

## Set up container objects
cds <- matrix(NA, nrow = nrow(pcData),
              ncol = (nsims * loop))

ncc <- rep(NA, (nsims * loop))
nacc <- rep(NA, (nsims * loop))
rsplit <- rep(NA, (nsims * loop))
rpar <- rep(NA, (nsims * loop))
radj <- rep(NA, (nsims * loop))
rlam <- rep(NA, (nsims * loop))
beta <- rep(NA, (nsims * loop))
betadiss <- rep(NA, (nsims * loop))
betapop <- rep(NA, (nsims * loop))
betaswitch <- rep(NA, (nsims * loop))
ssdvec <- rep(NA, (nsims * loop))
dissvec <- rep(NA, (nsims * loop))
popvec <- rep(NA, (nsims * loop))
switchvec <- rep(NA, (nsims * loop))
decisionvec <- rep(NA, (nsims * loop))

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

    ind <- ((i - 1) * nsims + 1):(i * nsims)
    
    ## Store objects together
    cds[1:nrow(pcData), ind] <- ecuts[[1]]

    ncc[ind] <- ecuts[[2]]
    nacc[ind] <- ecuts[[3]]
    rsplit[ind] <- ecuts[[4]]
    rpar[ind] <- ecuts[[5]]
    radj[ind] <- ecuts[[6]]
    rlam[ind] <- ecuts[[7]]
    beta[ind] <- ecuts[[8]]
    betadiss[ind] <- ecuts[[9]]
    betapop[ind] <- ecuts[[10]]
    betaswitch[ind] <- ecuts[[11]]
    ssdvec[ind] <- ecuts[[12]]
    dissvec[ind] <- ecuts[[13]]
    popvec[ind] <- ecuts[[14]]
    switchvec[ind] <- ecuts[[15]]
    decisionvec[ind] <- ecuts[[16]]

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

