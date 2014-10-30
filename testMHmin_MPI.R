#####################################
# Authors: Ben Fifield & Alex Tarr
# Created: 2014/09/20
# Last Revision: 2014/10/30
# Institution: Princeton University
# Purpose: Called by runSWA.R to run simulations
#####################################

params <- expand.grid(state = "OK_final",
                      eprob = 0.05, marginpct = 1,
                      lambda = 20, pnum = 1,
                      initbeta = 0,
                      initbetadiss = 0,
                      initbetapop = seq(0, -100, by = -10),
                      initbetaswitch = 0,
                      annealbeta = 0, annealbetadiss = 0,
                      annealbetapop = 0, annealbetaswitch = 0,
                      targbetapop = 0,
                      bybetapop = 0,
                      weightpow = 0,
                      nsims = 50, loop = 2,
                      wd = "/scratch/network/bfifield/segregation/data/",
                      dwd = "/scratch/network/bfifield/segregation/data/simRuns/",
                      codedir = "/scratch/network/bfifield/segregation/code/")

## Get state
state <- params$state[1]

## Modify pnum
params$pnum <- 1:nrow(params)

## Set working directory | state
setwd(paste(params$wd[1], state, sep = ""))

## Load packages, data
if (!is.loaded("mpi_initialize")) { 
    library("Rmpi") 
}

load(paste(getwd(), "/preprocPrec.RData", sep = ""))
load(paste(getwd(), "/pwdPrecinct.RData", sep = ""))

## Generate swapping sequence (in future, allow for user to input how often swaps are made)
nits <- params$nsims[1] * params$loop[1]
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

    library("redist"); library("BARD"); library("maptools"); library("methods")

    load(paste(getwd(), "/preprocPrec.RData", sep = ""))
    load(paste(getwd(), "/pwdPrecinct.RData", sep = ""))
    
    ######################
    ## Parameter values ##
    ######################

    ## Set state
    state <- params$state

    ## Set directories
    dwd <- params$dwd
    codedir <- params$codedir
    
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
    
    ## Weights power
    wpow <- params$weightpow
    
    parity <- sum(pcData$pop) / dists
    margin <- round(parity * margin.pct)
    
    ## Empty beta vector and weights for simulated tempering
    bvec <- rep(1, 5)  
    betaweights <- rep(1, length(bvec))
    
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
        nsimsAdj <- nsimsAdj[nsimsAdj > 0] # Corrects issue with swaps occurring on nsims*loop
        
        ## Get starting values and rng state in place
        if(i > 1){
            
            ## Load previous simulation for starting values
            load(paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                       "_", margin.pct, "_", lambda,
                       "_b", params$initbeta * -1, 
                       "_bDiss", params$initbetadiss, 
                       "_bPop", params$initbetapop * -1, 
                       "_bSwitch", params$initbetaswitch * -1,
                       "_pow", wpow,
                       "_par", pnum, "_loop", (i - 1), ".RData", sep = ""))
            
            ## Get starting values
            cds <- ecuts[[1]][,nsims]
            ## Set seed, use rng state from previous sims
            set.seed(1)
            .Random.seed <- ecuts[[20]]
            
            
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
                         beta = beta, betadiss = betadiss, betapop = betapop,
                         betaswitch = betaswitch,
                         betavec = bvec, betadissvec = bvec,
                         betapopvec = bvec, betaswitchvec = bvec,
                         betaweights = betaweights)
            
            ## Get likelihood for each of the constraints
            likeComp <- temp[[17]]
            likePop <- temp[[18]]
            likeDiss <- temp[[19]]
            
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
                mpi.send.Robj(betapop,dest=partner[j],tag=2)
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
                prob <- (likePop^betaPart*likePart^betapop)/(likePop^betapop*likePart^betaPart)
                if(prob > accept){
                    ## Exchange temperature values
                    betapop <- betaPart
                }           
            }
            
                                        # Update inputs to swMH
            cds <- ecuts[[1]][,nsimsAdj[j]]
            
        } 
        
        
        ## Save seed and simulations
        r <- .Random.seed
        
        ecuts[[20]] <- r
t        
        save(ecuts, file = paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                        "_", margin.pct, "_", lambda,
                        "_b", params$initbeta * -1, 
                        "_bDiss", params$initbetadiss, 
                        "_bPop", params$initbetapop * -1, 
                        "_bSwitch", params$initbetaswitch * -1,
                        "_pow", wpow,
                        "_par", pnum, "_loop", i, ".RData", sep = ""))
        
        print(paste("loop", i, sep = " "))
        
    }
}

## Spawn slaves equal to number of rows of params (one for each temperature)
mpi.spawn.Rslaves(nslaves = nrow(params))

## Get processor ID
mpi.bcast.cmd(procID<-mpi.comm.rank())

## Send swaps variable to each slave
mpi.bcast.Robj2slave(swaps)

## Send parameters to each slave (different row for each slave)
# mpi.bcast.cmd(params<-mpi.scatter.Robj2slave())
params <- split(params, f=1:nrow(params))
mpi.scatter.Robj2slave(params)

## Send ecutsMPI function to slaves
mpi.bcast.Robj2slave(ecutsMPI)

## Send ecutsAppend function to slaves
mpi.bcast.Robj2slave(ecutsAppend)

## Execute ecutsMPI program on each slave
mpi.remote.exec(swaps)
mpi.bcast.cmd(ecutsMPI())

## ## Combine save files (NEED TO SET DIRECTORY)
## mpi.bcast.cmd(source(paste(codedir, "combineDataSplit.R", sep = "")))

## ## Generate diagnostics
## mpi.bcast.cmd(source(paste(codedir, "genDiag.R", sep = "")))

## Close slaves
mpi.close.Rslaves()

## Terminate MPI processes and close R
mpi.quit()
