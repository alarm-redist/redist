#####################################
# Author: Ben Fifield
# Created: 2014/09/20
# Last Revision: 2014/09/21
# Institution: Princeton University
# Purpose: Called by runSWA.R to run simulations
#####################################

## Get parallel process
if(a == "Tukey"){
    aid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
    aid <- as.numeric(aid)
}
if(a == "Mac"){
    aid <- 1
}

## Get state
state <- params$state[aid]

## Set working directory | state
if(a == "Tukey"){
    setwd(paste("/scratch/network/bfifield/segregation/data/", state,
                "_final", sep = ""))
    dwd <- "/scratch/network/bfifield/segregation/data/simRuns/"
}
if(a == "Mac"){
    setwd(paste("~/Dropbox/Graduate School/Projects/segregation/svn/data/", state,
          "_final", sep = ""))
    dwd <- "~/Desktop/simRuns/"
}

## Load packages, data
library("redist"); library("BARD"); library("maptools"); library("methods")
load(paste(getwd(), "/preprocPrec.RData", sep = ""))
load(paste(getwd(), "/pwdPrecinct.RData", sep = ""))
        
######################
## Parameter values ##
######################
            
## District type
dists <- length(unique(pcData$cds))

## Edgecut prob - prob of turning edge off = 1 - eprob
eprob <- params$eprob[aid]

## Population constraint
margin.pct <- params$marginpct[aid]

## Lambda
lambda <- params$lambda[aid]

## Par num
pnum <- params$pnum[aid]

## Compactness beta
beta <- params$initbeta[aid]
betadiss <- params$initbetadiss[aid]
betapop <- params$initbetapop[aid]
betaswitch <- params$initbetaswitch[aid]

## Number of simulations
nsims <- params$nsims[aid]

## Repeat multiple times, save after
loop <- params$loop[aid]

## Annealing
abeta <- params$annealbeta[aid]
abetadiss <- params$annealbetadiss[aid]
abetapop <- params$annealbetapop[aid]
abetaswitch <- params$annealbetaswitch[aid]

## Weights power
wpow <- params$weightpow[aid]

parity <- sum(pcData$pop) / dists
margin <- round(parity * margin.pct)

## Beta vectors
betavec <- seq(0, -2000, by = -100)
betadissvec <- seq(0, 2000, by = 10)
betapopvec <- seq(0, params$targbetapop[aid], by = params$bybetapop[aid])
betaswitchvec <- seq(0, -1000, by = -1)

## Weights for beta
betaweights <- rep(NA, length(betapopvec))
for(i in 1:length(betaweights)){
    betaweights[i] <- wpow^i
}

#########################
## Run the simulations ##
#########################
for(i in 1:loop){

    ## Get starting values and rng state in place
    if(i > 1){
        
        ## Load previous simulation for starting values
        load(paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                   "_", margin.pct, "_", lambda,
                   "_b", params$initbeta[aid] * -1, "a", abeta,
                   "_bDiss", params$initbetadiss[aid], "a", abetadiss,
                   "_bPop", params$initbetapop[aid] * -1, "a", abetapop,
                   "t", params$targbetapop[aid] * -1,
                   "b", params$bybetapop[aid] * -1,
                   "_bSwitch", params$initbetaswitch[aid] * -1, "a", abetaswitch,
                   "_pow", wpow,
                   "_par", pnum, "_loop", (i - 1), ".RData", sep = ""))
        
        ## Get starting values
        cds <- ecuts[[1]][,nsims]
        if(abeta == 1){
            beta <- ecuts[[8]][nsims]
        }
        if(abetadiss == 1){
            betadiss <- ecuts[[9]][nsims]
        }
        if(abetapop == 1){
            betapop <- ecuts[[10]][nsims]
        }
        if(abetaswitch == 1){
            betaswitch <- ecuts[[11]][nsims]
        }
        ## Set seed, use rng state from previous sims
        set.seed(1)
        .Random.seed <- ecuts[[17]]
        
        
    } else{

        cds <- eval(parse(text = paste("pcData$cds", pnum, sep = "")))
        
    }
    
    #######################
    ## Run the algorithm ##
    #######################
    pcData$blackhisp <- pcData$BlackPop + pcData$HispPop

    set.seed(pnum)
    ecuts <- swMH(al.pc, cds, cds, nsims, eprob,
                  pcData$pop, pcData$blackhisp, parity, margin,
                  dists, lambda, pwdPrecinct,
                  beta, betadiss, betapop, betaswitch,
                  betavec, betadissvec, betapopvec, betaswitchvec, betaweights,
                  annealbeta = abeta, annealbetadiss = abetadiss,
                  annealbetapop = abetapop, annealbetaswitch = abetaswitch)

    ## Save seed and simulations
    r <- .Random.seed
        
    ecuts[[17]] <- r
    
    save(ecuts, file = paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                    "_", margin.pct, "_", lambda,
                    "_b", params$initbeta[aid] * -1, "a", abeta,
                    "_bDiss", params$initbetadiss[aid], "a", abetadiss,
                    "_bPop", params$initbetapop[aid] * -1, "a", abetapop,
                    "t", params$targbetapop[aid] * -1, "b",
                    params$bybetapop[aid] * -1,
                    "_bSwitch", params$initbetaswitch[aid] * -1, "a",
                    abetaswitch,
                    "_pow", wpow,
                    "_par", pnum, "_loop", i, ".RData", sep = ""))
    
    print(paste("loop", i, sep = " "))
    
}

