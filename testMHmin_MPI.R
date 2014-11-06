#####################################
# Authors: Ben Fifield & Alex Tarr
# Created: 2014/09/20
# Last Revision: 2014/10/30
# Institution: Princeton University
# Purpose: Called by runSWA.R to run simulations
#####################################

params <- expand.grid(state = "testset_25_2",
                      testset = TRUE,
                      eprob = 0.05, marginpct = 1,
                      lambda = 1, pnum = 1,
                      initbeta = 0,
                      initbetadiss = 0,
                      initbetapop = seq(0, -100, by = -10),
                      initbetaswitch = 0,
                      annealbeta = 0, annealbetadiss = 0,
                      annealbetapop = 0, annealbetaswitch = 0,
                      targbetapop = 0,
                      bybetapop = 0,
                      weightpow = 0,
                      nsims = 50000, loop = 1, thin = 1,
                      wd = "/scratch/network/bfifield/segregation/data/",
                      dwd = "/scratch/network/bfifield/segregation/data/simRuns/",
                      codedir = "/scratch/network/bfifield/segregation/code/redist-pkg/")

## Get state
state <- params$state[1]
testset <- params$testset[1]

## Modify pnum
params$pnum <- 1:nrow(params)

## Set working directory | state
if(testset == FALSE){
    setwd(paste(params$wd[1], state, sep = ""))
} else{
    setwd(paste(params$wd[1], "testsets/", sep = ""))
}
## Load packages, data
if (!is.loaded("mpi_initialize")) { 
    library("Rmpi") 
}

if(testset == FALSE){
    load(paste(getwd(), "/preprocPrec.RData", sep = ""))
    load(paste(getwd(), "/pwdPrecinct.RData", sep = ""))
} else{
    load(paste(getwd(), "/", params$state[1], ".RData", sep = ""))
}
## Generate swapping sequence (in future, allow for user to input how often swaps are made)
freq <- 10
nits <- params$nsims[1] * params$loop[1]
swaps <- matrix(NA,2,nits)
for(i in 1:nits){
    if(i %% freq == 0){
        swaps[,i] = sample(1:nrow(params),size=2)
    }
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

    if(params$testset == FALSE){
        load(paste(getwd(), "/preprocPrec.RData", sep = ""))
        load(paste(getwd(), "/pwdPrecinct.RData", sep = ""))
    } else{
        load(paste(getwd(), "/", params$state, ".RData", sep = ""))
    }
    ######################
    ## Parameter values ##
    ######################

    ## Set state
    state <- params$state

    ## Set directories
    dwd <- params$dwd
    codedir <- params$codedir
    
    ## District type
    if(params$testset == FALSE){
        dists <- length(unique(pcData$cds))
    } else{
        dists <- length(unique(cdmat[,1]))
    }
    
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

    ## Thinning amount
    thin <-  params$thin
    
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
    
    ## Test whether zero indexing causes error
    for(i in 1:length(al.pc)){
        al.pc[[i]] <- al.pc[[i]] - 1
    }
    
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
            ## set.seed(1)
            ## .Random.seed <- ecuts[[20]]
            
            
        } else{

            if(params$testset == FALSE){
                cds <- eval(parse(text = paste("pcData$cds", pnum, sep = "")))
            } else{
                cds <- cdmat[,sample(1:ncol(cdmat), 1)]
            }
        }

        ## Initialize ecuts object
        ecuts <- list()
        
        #######################
        ## Run the algorithm ##
        #######################
        pcData$blackhisp <- pcData$BlackPop + pcData$HispPop
        
        ## set.seed(pnum)
        
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

        nbetapop <- betapop
        save(ecuts, nbetapop,
             file = paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                        "_", margin.pct, "_", lambda,
                        "_b", params$initbeta * -1, 
                        "_bDiss", params$initbetadiss, 
                        "_bPop", params$initbetapop * -1, 
                        "_bSwitch", params$initbetaswitch * -1,
                        "_pow", wpow,
                        "_par", pnum, "_loop", i, ".RData", sep = ""))
        
        print(paste("loop", i, sep = " "))
        
    }

    ######################
    ## Combine the data ##
    ######################
    
    ## Set up container objects
    cds <- matrix(NA, nrow = nrow(pcData),
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

    ## Loop over simulation parameters ##
    for(i in 1:loop){

        ## Load data set
        load(paste(dwd, "ecuts", state, "_", (1 - eprob) * 100,
                   "_", margin.pct, "_", lambda,
                   "_b", params$initbeta * -1,
                   "_bDiss", params$initbetadiss,
                   "_bPop", params$initbetapop * -1,
                   "_bSwitch", params$initbetaswitch * -1,
                   "_pow", wpow,
                   "_par", pnum, "_loop", i, ".RData", sep = ""))

        ind <- ((i - 1) * (nsims / thin) + 1):(i * (nsims / thin))
        
        ## Store objects together
        cds[1:nrow(pcData), ind] <- ecuts[[1]][,indthin]

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
    save(edgecuts,
         file = paste(dwd, "edgecuts", state, "_", (1 - eprob) * 100,
                       "_", margin.pct, "_", lambda,
                       "_b", params$initbeta * -1,
                       "_bDiss", params$initbetadiss,
                       "_bPop", nbetapop, 
                       "_bSwitch", params$initbetaswitch * -1, 
                       "_pow", wpow, ".RData", sep = ""))

    ###############################
    ## Generate diagnostic plots ##
    ###############################

    ## Number of acceptances
    accept <- 0
    for(z in 2:length(edgecuts[[10]])){
        if(edgecuts[[10]][z-1] != edgecuts[[10]][z]){
            accept <- accept + 1
        }
    }

    ## Distance from population parity
    paritydist <- distParity(edgecuts[[1]], pcData$pop, parity)

    ## For subsetting accept/reject plot
    ind <- seq(1, length(edgecuts[[16]]), by = 25)

    ## Diagnostic plot ##
    pdf(file = paste(dwd, "diagPlot", state, "_", (1 - eprob) * 100,
        "_", margin.pct, "_", lambda,
        "_b", params$initbeta * -1, 
        "_bDiss", params$initbetadiss,
        "_bPop", params$initbetapop * -1,
        "_bSwitch", params$initbetaswitch * -1,
        "_pow", wpow,
        "_par", pnum, ".pdf", sep = ""),
        height = 6, width = 12)
    par(mfrow = c(1,2))
    ## Trace of beta
    plot(edgecuts[[10]], type = "l",
         xlab = "Iterations", ylab = "",
         main = paste("Trace of Beta \n Lambda =",
             params$lambda, "pow =",
             params$weightpow, sep = " "))
    ## Distribution of beta values
    barplot(table(edgecuts[[10]]),
            main = paste("Distribution of Beta_pop \n Acceptance Rate =",
                accept / length(edgecuts[[10]]), sep = " "))

    ## Accept/Reject
    plot(edgecuts[[16]][ind], pch = 16,
         cex = .5,
         yaxt = "n",
         col = ifelse(edgecuts[[16]][ind] == 0, "red", "black"),
         xlab = "Iterations", ylab = "",
         main = paste("Rejections from SWA \n Acceptance Rate =",
             sum(edgecuts[[16]]) / length(edgecuts[[16]]), sep = " "))
    axis(2, c(0, 1), c("Reject", "Accept"))
    ## Trace of district population
    plot(paritydist, type = "l",
         xlab = "Iterations", ylab = "",
         main = "Trace of Distance from Parity")
    dev.off() 
    
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
mpi.bcast.cmd(ecutsMPI())

## Combine save files (NEED TO SET DIRECTORY)
## mpi.bcast.cmd(source(paste(params$codedir[1], "combineDataSplit.R", sep = "")))

## Generate diagnostics
## mpi.bcast.cmd(source(paste(params$codedir[1], "genDiag.R", sep = "")))

## Close slaves
mpi.close.Rslaves()

## Terminate MPI processes and close R
mpi.quit()
