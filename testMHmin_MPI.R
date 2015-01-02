#####################################
# Authors: Ben Fifield & Alex Tarr
# Created: 2014/09/20
# Last Revision: 2014/10/30
# Institution: Princeton University
# Purpose: Called by runSWA.R to run simulations
#####################################
aid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

betaseq <- rep(NA, 11)
for(i in 1:11){
    betaseq[i] <- -(0.1^((i-1) / (length(betaseq) - 1)))
}
betaseq <- 400 * (betaseq + .1)

params <- expand.grid(state = "ms", adjSwap = TRUE, freq = 1000, 
                      eprob = 0.05, marginpct = 1,
                      lambda = 18, pnum = 1,
                      initbeta = 0,
                      initbetadiss = 0,
                      initbetapop = betaseq,
                      initbetaswitch = 0,
                      annealbeta = 0, annealbetadiss = 0,
                      annealbetapop = 1, annealbetaswitch = 0,
                      targbetapop = 0,
                      bybetapop = 0,
                      weightpow = 0,
                      nsims = 10000, loop = 10, thin = 1,
                      wd = "/scratch/network/bfifield/segregation/data/",
                      logdir = "/scratch/network/bfifield/segregation/code/slurm/",
                      dwd = "/scratch/network/bfifield/segregation/data/simRuns/",
                      codedir = "/scratch/network/bfifield/segregation/code/redist-pkg/")

## Get state
state <- params$state[1]

## Adjacent swap flag
adjSwap <- params$adjSwap[1]

## Modify pnum
params$pnum <- 1:nrow(params)

## Set working directory | state
setwd(paste(params$wd[1], state, sep = ""))

## Load packages, data
if (!is.loaded("mpi_initialize")) { 
    library("Rmpi") 
}

load(paste(getwd(), "/algdat.RData", sep = ""))

## Generate swapping sequence (in future, allow for user to input how often swaps are made)
freq <- params$freq[1]
nits <- params$nsims[1] * params$loop[1]

if(adjSwap){
  swaps <- matrix(NA,1,nits)
  partner <- matrix(NA,1,nits)
  for(i in 1:nits){
    if(i %% freq == 0){
      swaps[1,i] = sample(1:nrow(params),size=1)
      swaps[2,i] = sample(c(-2,0),size=1)
    }
  }
}
else{
  swaps <- matrix(NA,2,nits)
  for(i in 1:nits){
      if(i %% freq == 0){
          swaps[,i] = sample(1:nrow(params),size=2)
      }
  }
}

## Initial temperature adjacency list
tempAdj <- 1:nrow(params)

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

    aid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

    library("redist"); library("BARD"); library("maptools"); library("methods")

    load(paste(getwd(), "/algdat.RData", sep = ""))

    ######################
    ## Parameter values ##
    ######################

    ## Set state
    state <- params$state
    
    ## Set adjacent swap flag
    adjSwap <- params$adjSwap
    
    ## Swapping frequency
    freq <- params$freq

    ## Set directories
    dwd <- params$dwd
    codedir <- params$codedir
    logdir <- params$logdir
    
    ## District type
    if(substr(state, 1, 7) == "testset"){
        dists <- length(unique(cdmat[,1]))
    } else{
        dists <- length(unique(geodat$cds))
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
    
    parity <- sum(geodat$pop) / dists
    margin <- round(parity * margin.pct)
    
    ## Empty beta vector and weights for simulated tempering
    bvec <- rep(1, 5)  
    betaweights <- rep(1, length(bvec))
    
    #######################################
    ## Additional Pre-Processing for MPI ##  
    #######################################
    
    ## Find iterations for which a swap is proposed involving process procID
    swapIts <- which(swaps == procID, arr.ind = TRUE)[,2]
      
    #########################
    ## Run the simulations ##
    #########################

    ## Filename for storing progress
    fname <- paste(logdir,"progress_log",procID,"_par",aid,sep="")
    
    for(i in 1:loop){
        
      if(adjSwap){
        nsimsAdj <- rep(freq,nsims/freq)
      }
      else{
        ## Construct adjusted "nsims" vector
        tempIts <- swapIts[swapIts <= nsims*i & swapIts > nsims*(i-1)]
        ## Swap partners
        partner <- swaps[,tempIts][swaps[,tempIts] != procID]
        nsimsAdj <- c(tempIts,nsims*i) - c((i-1)*nsims,tempIts)
        nsimsAdj <- nsimsAdj[nsimsAdj > 0] # Corrects issue with swaps occurring on nsims*loop
      }
        ## Get starting values and rng state in place
        if(i > 1){
            
            ## Load previous simulation for starting values
            load(paste(dwd, "ecutsMPI", state, "_", (1 - eprob) * 100,
                       "_", margin.pct, "_", lambda,
                       "_b", params$initbeta * -1, 
                       "_bDiss", params$initbetadiss, 
                       "_bPop", params$initbetapop * -1, 
                       "_bSwitch", params$initbetaswitch * -1,
                       "_pow", wpow,
                       "_par", aid, "_loop", (i - 1), ".RData", sep = ""))
            
            ## Get starting values
            cds <- ecuts[[1]][,nsims]

            ## Set seed, use rng state from previous sims
            ## set.seed(1)
            ## .Random.seed <- ecuts[[20]]
            
            
        } else{

            if(substr(state, 1, 7) == "testset"){
                cds <- cdmat[,sample(1:ncol(cdmat), 1)]
            } else{
                cds <- eval(parse(text = paste("geodat$cds", aid, sep = "")))
            }
        }

        ## Initialize ecuts object
        ecuts <- list()
        
        #######################
        ## Run the algorithm ##
        #######################
        geodat$blackhisp <- geodat$BlackPop + geodat$HispPop
        
        set.seed(pnum)
        sink(fname)
        cat("Length nsimsAdj = ",length(nsimsAdj), "\n", append = TRUE)
        for(j in 1:length(nsimsAdj)){

            cat("Iteration ", j, "\n", append = TRUE)

            cat("Start swMH.\n", append = TRUE)
            temp <- swMH(al.pc, cds, cds, nsimsAdj[j], eprob,
                         geodat$pop, geodat$blackhisp, parity, margin,
                         dists, lambda, ssdmat,
                         beta = beta, betadiss = betadiss, betapop = betapop,
                         betaswitch = betaswitch,
                         betavec = bvec, betadissvec = bvec,
                         betapopvec = bvec, betaswitchvec = bvec,
                         betaweights = betaweights)
            cat("End swMH.\n", append = TRUE)
            
            ## Get likelihood for each of the constraints
            likeComp <- temp[[17]]
            likePop <- temp[[18]]
            likeDiss <- temp[[19]]
            
            ## Update ecuts object
            cat("Start Append.\n", append = TRUE)
            ecuts <- ecutsAppend(ecuts,temp)
            cat("End Append.\n", append = TRUE)
            
            ## Use MPI to exchange swap information
            ##
            ## Tag guide: 1 -> likelihood sent
            ##            2 -> temperature sent
            ##            3 -> acceptance probability sent
            ##
            ## Note: Need to allow for multiple betas and likelihoods. As is, only likePop is
            ##       communicated and tested
            
            if(adjSwap){
              ## Determine swapping partner
              tempPart <- swaps[2,(i-1)*nsims + j*freq]
              
              ## Account for edge cases if necessary
              if(tempAdj[1] == swaps[1,(i-1)*nsims + j*freq]){
                tempPart <- tempAdj[2] #Edge case
              }
              else if(tempAdj[length(tempAdj)] == swaps[1,(i-1)*nsims + j*freq]){
                tempPart <- tempAdj[length(tempAdj)-1] #Edge case
              }
              else{
                tempInd <- which(tempAdj == swaps[1,(i-1)*nsims + j*freq],arr.ind=TRUE)
                tempPart <- tempAdj[tempInd + tempPart + 1]
              }
              
              ## Set partner variable
              if(swaps[1,(i-1)*nsims + j*freq] == procID){
                partner <- tempPart
              }
              else if(tempPart == procID){
                partner <- swaps[1,(i-1)*nsims + j*freq]
              }
              
              ## Communication step
              if(swaps[1,(i-1)*nsims + j*freq] == procID || tempPart == procID){
                cat("Start send to ", partner, ".\n", append = TRUE)
                ## Send commands (blocking)
                mpi.send.Robj(likePop,dest=partner,tag=1)
                mpi.send.Robj(betapop,dest=partner,tag=2)
                cat("End send to ", partner, ".\n", append = TRUE)
                ## Receive commands (blocking)
                cat("Start receive to ", partner, ".\n", append = TRUE)
                likePart <- mpi.recv.Robj(partner,tag=1)
                betaPart <- mpi.recv.Robj(partner,tag=2)
                cat("End receive to ", partner, ".\n", append = TRUE)
                
                ## Higher ranked process communicates random draw to lower ranked process
                cat("Start mh step.\n", append = TRUE)
                if(partner < procID){
                  accept <- runif(1)
                  mpi.send.Robj(accept,dest=partner,tag=3)
                }
                else{
                  accept <- mpi.recv.Robj(partner,tag=3)
                }
                
                ## Compute acceptance probability (for now, population only)
                prob <- (likePop^betaPart*likePart^betapop)/(likePop^betapop*likePart^betaPart)
                if(prob > accept){
                  ## Exchange temperature values
                  betapop <- betaPart
                  
                  ## Adjust temperature adjacency list and distribute among processes
                  tempAdj[tempAdj == tempPart] <- swaps[1,(i-1)*nsims + j*freq]
                  tempAdj[tempAdj == swaps[1,(i-1)*nsims + j*freq]] <- tempPart
                }
                cat("End mh step.\n", append = TRUE)
                
                ## Send temperature adjacency list
                if(swaps[1,(i-1)*nsims + j*freq] == procID){
                  oProcs <- tempAdj[tempAdj != tempPart & tempAdj != swaps[1,(i-1)*nsims + j*freq]]
                  for(k in 1:length(oProcs)){
                    mpi.send.Robj(tempAdj,dest=oProcs[k],tag=4)
                  }
                }
              }
              else{
                tempAdj <- mpi.recv.Robj(swaps[1,(i-1)*nsims + j*freq],tag=4)
              }
            }
            
            else{
              if(j != length(nsimsAdj) || length(nsimsAdj) == length(tempIts)){ # Swap proposed
                  cat("Start send to ", partner[j], ".\n", append = TRUE)
                  ## Send commands (blocking)
                  mpi.send.Robj(likePop,dest=partner[j],tag=1)
                  mpi.send.Robj(betapop,dest=partner[j],tag=2)
                  cat("End send to ", partner[j], ".\n", append = TRUE)
                  ## Receive commands (blocking)
                  cat("Start receive to ", partner[j], ".\n", append = TRUE)
                  likePart <- mpi.recv.Robj(partner[j],tag=1)
                  betaPart <- mpi.recv.Robj(partner[j],tag=2)
                  cat("End receive to ", partner[j], ".\n", append = TRUE)
                
                  ## Higher ranked process communicates random draw to lower ranked process
                  cat("Start mh step.\n", append = TRUE)
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
                  cat("End mh step.\n", append = TRUE)
              }
            }
            
            ## Update inputs to swMH
            cds <- ecuts[[1]][,nsimsAdj[j]]

            ## Write progress to file
            prog <- ( (i-1)*nsims + sum(nsimsAdj[1:j]) )/(loop*nsims/100) 
            cat(prog,"% of task on processor",procID,"has completed.\n",
                append = TRUE)
            
        } 
        sink() 
        
        ## Save seed and simulations
        r <- .Random.seed
        
        ecuts[[20]] <- r

        nbetapop <- betapop
        save(ecuts, nbetapop,
             file = paste(dwd, "ecutsMPI", state, "_", (1 - eprob) * 100,
                        "_", margin.pct, "_", lambda,
                        "_b", params$initbeta * -1, 
                        "_bDiss", params$initbetadiss, 
                        "_bPop", params$initbetapop * -1, 
                        "_bSwitch", params$initbetaswitch * -1,
                        "_pow", wpow,
                        "_par", aid, "_loop", i, ".RData", sep = ""))
        
        print(paste("loop", i, sep = " "))
        
    }

    ######################
    ## Combine the data ##
    ######################
    
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

    ## Loop over simulation parameters ##
    for(i in 1:loop){

        ## Load data set
        load(paste(dwd, "ecutsMPI", state, "_", (1 - eprob) * 100,
                   "_", margin.pct, "_", lambda,
                   "_b", params$initbeta * -1,
                   "_bDiss", params$initbetadiss,
                   "_bPop", params$initbetapop * -1,
                   "_bSwitch", params$initbetaswitch * -1,
                   "_pow", wpow,
                   "_par", aid, "_loop", i, ".RData", sep = ""))

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
    save(edgecuts,
         file = paste(dwd, "edgecutsMPI", state, "_", (1 - eprob) * 100,
             "_", margin.pct, "_", lambda,
             "_b", params$initbeta * -1,
             "_bDiss", params$initbetadiss,
             "_bPop", nbetapop, 
             "_bSwitch", params$initbetaswitch * -1, 
             "_pow", wpow, "_par", aid,
             ".RData", sep = ""))

    ###############################
    ## Generate diagnostic plots ##
    ###############################

    if(substr(state, 1, 7) == "testset"){

        ##     ## Number of acceptances
        ##     accept <- 0
        ##     for(z in 2:length(edgecuts[[10]])){
        ##         if(edgecuts[[10]][z-1] != edgecuts[[10]][z]){
        ##             accept <- accept + 1
        ##         }
        ##     }

        ##     ## Distance from population parity
        ##     paritydist <- distParity(edgecuts[[1]], geodat$pop, parity)

        ##     ## For subsetting accept/reject plot
        ##     ind <- seq(1, length(edgecuts[[16]]), by = 25)

        simdist <- sumstat(edgecuts[[1]],
                           geodat$HispPop,
                           geodat$BlackPop,
                           geodat$obama,
                           geodat$mccain,
                           geodat$pop)
        
        ## Diagnostic plot ##
        pdf(file = paste(dwd, "diagPlotMPI", state, "_", (1 - eprob) * 100,
                "_", margin.pct, "_", lambda,
                "_b", params$initbeta * -1, 
                "_bDiss", params$initbetadiss,
                "_bPop", params$initbetapop * -1,
                "_bSwitch", params$initbetaswitch * -1,
                "_pow", wpow, "_par", aid, ".pdf", sep = ""),
            height = 6, width = 12)
        ## par(mfrow = c(1,2))
        ## ## Trace of beta
        ## plot(edgecuts[[10]], type = "l",
        ##      xlab = "Iterations", ylab = "",
        ##      main = paste("Trace of Beta \n Lambda =",
        ##          params$lambda, "pow =",
        ##          params$weightpow, sep = " "))
        ## ## Distribution of beta values
        ## barplot(table(edgecuts[[10]]),
        ##         main = paste("Distribution of Beta_pop \n Acceptance Rate =",
        ##             accept / length(edgecuts[[10]]), sep = " "))

        ## ## Accept/Reject
        ## plot(edgecuts[[16]][ind], pch = 16,
        ##      cex = .5,
        ##      yaxt = "n",
        ##      col = ifelse(edgecuts[[16]][ind] == 0, "red", "black"),
        ##      xlab = "Iterations", ylab = "",
        ##      main = paste("Rejections from SWA \n Acceptance Rate =",
        ##          sum(edgecuts[[16]]) / length(edgecuts[[16]]), sep = " "))
        ## axis(2, c(0, 1), c("Reject", "Accept"))
        ## ## Trace of district population
        ## plot(paritydist, type = "l",
        ##      xlab = "Iterations", ylab = "",
        ##      main = "Trace of Distance from Parity")

        ## compare to null distribution if a test set
        par(mfrow = c(2,2))
        hist(nulldist[,1], freq = FALSE,
             ylim = c(0,15),
             xlab = "",
             ylab = "Hispanic Dissimilarity",
             breaks = 20)
        lines(density(simdist[,1], from = 0))
        hist(nulldist[,2], freq = FALSE,
             ylim = c(0,15),
             xlab = "",
             ylab = "Afam Dissimilarity",
             breaks = 20)
        lines(density(simdist[,2], from = 0))
        hist(nulldist[,3], freq = FALSE,
             ylim = c(0,15),
             xlab = "",
             ylab = "Democratic Dissimilarity",
             breaks = 20)
        lines(density(simdist[,3], from = 0))
        hist(nulldist[,4], freq = FALSE,
             ylim = c(0,15),
             xlab = "",
             ylab = "Republican Dissimilarity",
             breaks = 20)
        lines(density(simdist[,4], from = 0))
        dev.off()
    }

}

## Spawn slaves equal to number of rows of params (one for each temperature)
mpi.spawn.Rslaves(nslaves = nrow(params))

## Get processor ID
mpi.bcast.cmd(procID<-mpi.comm.rank())

## Send swaps variable to each slave
mpi.bcast.Robj2slave(swaps)

## Send temperature adjacency if doing adjacent swaps
if(params$adjSwap[1]){
  mpi.bcastRobj2slave(tempAdj)
}

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
