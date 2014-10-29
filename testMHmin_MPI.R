#####################################
# Authors: Ben Fifield & Alex Tarr
# Created: 2014/09/20
# Last Revision: 2014/10/30
# Institution: Princeton University
# Purpose: Called by runSWA.R to run simulations
#####################################

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

## Get state
state <- params$state

## Set working directory | state
if(a == "Tukey"){
  setwd(paste("/scratch/network/bfifield/segregation/data/", state,
              "_final", sep = ""))
  dwd <- "/scratch/network/bfifield/segregation/data/simRuns/"
}

## Load packages, data
library("redist"); library("BARD"); library("maptools"); library("methods");
if (!is.loaded("mpi_initialize")) { 
  library("Rmpi") 
}

load(paste(getwd(), "/preprocPrec.RData", sep = ""))
load(paste(getwd(), "/pwdPrecinct.RData", sep = ""))

# Generate swapping sequence (in future, allow for user to input how often swaps are made)
nits <- params[1,ncols(params)-1]*params[1,ncols(params)]
swaps <- matrix(0,2,nits)
for(i in 1:nits){
  swaps[,i] = sample(1:nrow(params),size=2)
}

## Functions

ecutsAppend <- function(ecuts,ndata){
  if(length(ecuts) == 0){
    ecuts <- ndata
  }
  else{
    ecuts[[1]] <- cbind(ecuts[[1]],ndata[[1]])
    ecuts[[2]] <- cbind(ecuts[[2]],ndata[[2]])
    ecuts[[3]] <- cbind(ecuts[[3]],ndata[[3]])
    ecuts[[4]] <- cbind(ecuts[[4]],ndata[[4]])
    ecuts[[5]] <- cbind(ecuts[[5]],ndata[[5]])
    ecuts[[6]] <- cbind(ecuts[[6]],ndata[[6]])
    ecuts[[7]] <- cbind(ecuts[[7]],ndata[[7]])
    ecuts[[8]] <- cbind(ecuts[[8]],ndata[[8]])
    ecuts[[9]] <- cbind(ecuts[[9]],ndata[[9]])
    ecuts[[10]] <- cbind(ecuts[[10]],ndata[[10]])
    ecuts[[11]] <- cbind(ecuts[[11]],ndata[[11]])
    ecuts[[12]] <- cbind(ecuts[[12]],ndata[[12]])
    ecuts[[13]] <- cbind(ecuts[[13]],ndata[[13]])
    ecuts[[14]] <- cbind(ecuts[[14]],ndata[[14]])
    ecuts[[15]] <- cbind(ecuts[[15]],ndata[[15]])
    ecuts[[16]] <- cbind(ecuts[[16]],ndata[[16]])
  }
  return(ecuts)
}

ecutsMPI <- function(){
  load(paste(getwd(), "/preprocPrec.RData", sep = ""))
  load(paste(getwd(), "/pwdPrecinct.RData", sep = ""))
  
  ######################
  ## Parameter values ##
  ######################
  
  ## District type
  dists <- length(unique(pcData$cds))
  
  ## Edgecut prob - prob of turning edge off = 1 - eprob
  eprob <- params$eprob
  
  ## Population constraint
  margin.pct <- params$marginpct
  
  ## Lambda
  lambda <- params$lambda
  
  ## Par num
  pnum <- params$pnum
  
  ## Compactness beta
  beta <- params$initbeta
  betadiss <- params$initbetadiss
  betapop <- params$initbetapop
  betaswitch <- params$initbetaswitch
  
  ## Number of simulations
  nsims <- params$nsims
  
  ## Repeat multiple times, save after
  loop <- params$loop
  
  ## Annealing
  abeta <- params$annealbeta
  abetadiss <- params$annealbetadiss
  abetapop <- params$annealbetapop
  abetaswitch <- params$annealbetaswitch
  
  ## Weights power
  wpow <- params$weightpow
  
  parity <- sum(pcData$pop) / dists
  margin <- round(parity * margin.pct)
  
  ## Beta vectors (NEEDS TO BE CHANGED FOR MPI IMPLEMENTATION)
  betavec <- seq(0, -2000, by = -100)
  betadissvec <- seq(0, 2000, by = 10)
  betapopvec <- seq(0, -450, by = -15)
  betaswitchvec <- seq(0, -1000, by = -1)
  
  betaweights <- rep(NA, length(betapopvec))
  for(i in 1:length(betaweights)){
    betaweights[i] <- wpow^i
  }
  
  #######################################
  ## Additional Pre-Processing for MPI ##  
  #######################################
  
  ## Find iterations for which a swap is proposed involving process procID
  swapIts <- which(swaps == procID, arr.ind = TRUE)[,2]
  
  ## Swap partners
  partner <- swaps[,swapIts][swaps[,swapIts] != procID]
  
  ## Initialize ecuts object
  ecuts <- list()
  
  
  #########################
  ## Run the simulations ##
  #########################
  for(i in 1:loop){
    
    ## Construct adjusted "nsims" vector
    temp <- swapIts[swapIts <= nsims*loop && swapIts > nsims*(loop-1)]
    nsimsAdj <- c(temp,nsims*loop) - c((loop-1)*nsims,temp)
    nsimsAdj[nsimsAdj > 0] # Corrects issue with swaps occurring on nsims*loop
    
    ## Get starting values and rng state in place
    if(i > 1){
      
      ## Load previous simulation for starting values
      load(paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                 "_", margin.pct, "_", lambda,
                 "_b", params$initbeta * -1, "a", abeta,
                 "_bDiss", params$initbetadiss, "a", abetadiss,
                 "_bPop", params$initbetapop * -1, "a", abetapop,
                 "t", params$targbetapop * -1,
                 "b", params$bybetapop * -1,
                 "_bSwitch", params$initbetaswitch * -1, "a", abetaswitch,
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
      .Random.seed <- ecuts[[19]]
      
      
    } else{
      
      cds <- eval(parse(text = paste("pcData$cds", pnum, sep = "")))
      
    }
    
    #######################
    ## Run the algorithm ##
    #######################
    pcData$blackhisp <- pcData$BlackPop + pcData$HispPop
    
    set.seed(pnum)
    
    for(j in 1:length(nsimsAdj)){
      
      temp <- swMH(al.pc, cds, cds, nsimsAdj[j], eprob,
                    pcData$pop, pcData$blackhisp, parity, margin,
                    dists, lambda, pwdPrecinct,
                    beta, betadiss, betapop, betaswitch,
                    betavec, betadissvec, betapopvec, betaswitchvec, betaweights,
                    annealbeta = abeta, annealbetadiss = abetadiss,
                    annealbetapop = abetapop, annealbetaswitch = abetaswitch)
      
      ## Get likelihood for each of the constraints
      likePop <- temp[[17]]
      likeDiss <- temp[[18]]
      
      ## Update ecuts object
      ecuts <- ecutsAppend(ecuts,temp)
      
      
      ## Use MPI to exchange swap information
      ##
      ## Tag guide: 1 -> likelihood sent
      ##            2 -> temperature sent
      ##            3 -> acceptance probability sent
      ##
      ## Note: Need to allow for multiple betas and likelihoods. As is, only likePop is
      ##       communicated and tested
      
      if(j != length(nsimsAdj) || length(nsimsAdj) == length(temp)){ # Swap proposed
        ## Send commands (blocking)
        mpi.send.Robj(likePop,dest=partner[j],tag=1)
        mpi.send.Robj(betaPop,dest=partner[j],tag=2)
        ## Receive commands (blocking)
        likePart <- mpi.recv.Robj(partner[j],tag=1)
        betaPart <- mpi.recv.Robj(partner[j],tag=2)
        
        ## Higher ranked process communicates random draw to lower ranked process
        if(partner[j] < procID){
          accept <- runif(1)
          mpi.send.Robj(accept,dest=partner[j],tag=3)
        }
        else{
          accept <- mpi.recv.Robj(partner[j],tag=3)
        }
        
        ## Compute acceptance probability (for now, population only)
        prob <- (likePop^betaPart*likePart^betaPop)/(likePop^betaPop*likePart^betaPart)
        if(prob > accept){
          ## Exchange temperature values
          betaPop <- betaPart
        }           
      }
      
      # Update inputs to swMH
      cds <- ecuts[[1]][,nsimAdj(j)]
      
    } 
    
    
    ## Save seed and simulations
    r <- .Random.seed
    
    ecuts[[19]] <- r
    
    save(ecuts, file = paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                             "_", margin.pct, "_", lambda,
                             "_b", params$initbeta * -1, "a", abeta,
                             "_bDiss", params$initbetadiss, "a", abetadiss,
                             "_bPop", params$initbetapop * -1, "a", abetapop,
                             "t", params$targbetapop * -1, "b",
                             params$bybetapop * -1,
                             "_bSwitch", params$initbetaswitch * -1, "a",
                             abetaswitch,
                             "_pow", wpow,
                             "_par", pnum, "_loop", i, ".RData", sep = ""))
    
    print(paste("loop", i, sep = " "))
    
  }
}

## Spawn slaves equal to number of rows of params (one for each temperature)
mpi.spawn.Rslaves(nslaves=nrow(params))

## Get processor ID
mpi.bcast.cmd(procID<-mpi.comm.rank())

## Send swaps variable to each slave
mpi.bcast.Robj2slave(swaps)
## Send parameters to each slave (different row for each slave)
mpi.scatter.Robj2slave(split(params,f=1:nrow(params)))

## Send ecutsMPI function to slaves
mpi.bcast.Robj2slave(ecutsMPI)

## Send ecutsAppend function to slaves
mpi.bcast.Robj2slave(ecutsAppend)

## Execute ecutsMPI program on each slave
mpi.remote.exec(ecutsMPI())

## Close slaves
mpi.close.Rslaves()

## Combine save files (NEED TO SET DIRECTORY)
source(paste(codedir, "combineDataSplit.R", sep = ""))

## Generate diagnostics
source(paste(codedir, "genDiag.R", sep = ""))

## Terminate MPI processes and close R
mpi.quit()
