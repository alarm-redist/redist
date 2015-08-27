##############################################
# Authors: Ben Fifield & Alex Tarr           #
# Created: 2015/06/01                        #
# Last Revision: N/A                         #
# Institution: Princeton University          #
# Purpose: R wrapper to run swMH() code w/   #
# parallel tempering (mpi)                   #
##############################################

ecutsMPI <- function(){
    ## Load redist library
    library(redist)

    fname <- paste("log", procID, sep = "")
    
    if(params$verbose){
        sink(fname)
        ## Initialize ##
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")
        
        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        cat("redist.mcmc.mpi(): Automated Redistricting Simulation Using
        Markov Chain Monte Carlo w/ Parallel Tempering \n\n", append = TRUE)
    }
    
    ## Extract variables
    if(is.na(grouppopvec)){
        grouppopvec <- NULL
    }
    if(is.na(ssdmat)){
        ssdmat <- NULL
    }
    if(is.na(params$adjswaps)){
        adjswaps <- NULL
    }else{
        adjswaps <- params$adjswaps
    }
    if(is.na(params$freq)){
        freq <- NULL
    }else{
        freq <- params$freq
    }
    if(is.na(params$constraint)){
        constraint <- NULL
    }else{
        constraint <- params$constraint
    }
    if(is.na(params$nsims)){
        nsims <- NULL
    }else{
        nsims <- params$nsims
    }
    if(is.na(params$nloop)){
        nloop <- NULL
    }else{
        nloop <- params$nloop
    }
    if(is.na(params$eprob)){
        eprob <- NULL
    }else{
        eprob <- params$eprob
    }
    if(is.na(params$popcons)){
        popcons <- NULL
    }else{
        popcons <- params$popcons
    }
    if(is.na(params$lambda)){
        lambda <- NULL
    }else{
        lambda <- params$lambda
    }
    if(is.na(params$maxiterrsg)){
        maxiterrsg <- NULL
    }else{
        maxiterrsg <- params$maxiterrsg
    }
    if(is.na(params$contiguitymap)){
        contiguitymap <- NULL
    }else{
        contiguitymap <- params$contiguitymap
    }
    if(is.na(params$loopscompleted)){
        loopscompleted <- NULL
    }else{
        loopscompleted <- params$loopscompleted
    }
    if(is.na(params$rngseed)){
        rngseed <- NULL
    }else{
        rngseed <- params$rngseed
    }
    if(is.na(params$ndists)){
        ndists <- NULL
    }else{
        ndists <- params$ndists
    }
    if(is.na(params$savename)){
        savename <- NULL
    }else{
        savename <- params$savename
    }
    ## Run redist preprocessing function
    preprocout <- redist.preproc(adjobj = adjobj, popvec = popvec,
                                 initcds = initcds, ndists = ndists,
                                 popcons = popcons,
                                 grouppopvec = grouppopvec, ssdmat = ssdmat,
                                 beta = params$beta, temper = "parallel",
                                 constraint = constraint,
                                 betaseq = NULL, betaweights = NULL,
                                 adjswaps = adjswaps,
                                 maxiterrsg = maxiterrsg,
                                 contiguitymap = contiguitymap)

    ## Set betas - if tempering, modified later
    betapop <- preprocout$params$betapop
    betacompact <- preprocout$params$betacompact
    betaseg <- preprocout$params$betaseg
    betasimilar <- preprocout$params$betasimilar
    temper <- "parallel"
    
    ## Find procID involved in swaps (non-adjacent only)
    if(!adjswaps){
        swapIts <- which(swaps == procID, arr.ind = TRUE)[,2]
    }
    
    ## Set seed before first iteration of algorithm if provided by user
    if(!is.null(rngseed) & is.numeric(rngseed)){
        set.seed(rngseed)
    }
    
    ## Get starting loop value
    loopstart <- loopscompleted + 1

    for(i in loopstart:nloop){
        
        if(adjswaps){
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

        ## Get initial partition
        if(i > loopstart){
            
            cds <- algout$partitions[,nsims]
            
            if(temper != "none" & constraint == "compact"){
                betacompact <- algout$beta_sequence[nsims]
            }
            if(temper != "none" & constraint == "segregation"){
                betaseg <- algout$beta_sequence[nsims]
            }
            if(temper != "none" & constraint == "population"){
                betapop <- algout$beta_sequence[nsims]
            }
            if(temper != "none" & constraint == "similarity"){
                betasimilar <- algout$beta_sequence[nsims]
            }
            if(!is.null(rngseed) & is.numeric(rngseed)){
                .Random.seed <- algout$randseed
            }
            
            rm(list = "algout")
            
        } else{ 
            
            ## Reload the data if restarting
            if(loopstart > 1){
                
                ## Load the data (CHANGE)
                load(paste(savename, "_loop", i - 1, "_chain", procID, ".RData", sep = ""))
                
                ## NOTE tempadj NEEDS TO BE PART OF SAVE DATA
                
                ## Stop if number of simulations per loop is different
                if(nsims != ncol(algout[[1]])){
                    stop("Please specify the same number of simulations per
               loop across all loops")
                }
                
                cds <- algout$partitions[,nsims]
                
                if(temper != "none" & constraint == "compact"){
                    betacompact <- algout$beta_sequence[nsims]
                }
                if(temper != "none" & constraint == "segregation"){
                    betaseg <- algout$beta_sequence[nsims]
                }
                if(temper != "none" & constraint == "population"){
                    betapop <- algout$beta_sequence[nsims]
                }
                if(temper != "none" & constraint == "similarity"){
                    betasimilar <- algout$beta_sequence[nsims]
                }
                if(!is.null(rngseed) & is.numeric(rngseed)){
                    .Random.seed <- algout$randseed
                }
                
                rm(list = "algout")
                
            }else{
                cds <- preprocout$data$initcds
                ## Initialize algout object (for use in ecutsAppend)
                algout <- list()
            }
            
        } 
        
        #######################
        ## Run the algorithm ##
        #######################
        for(j in 1:length(nsimsAdj)){

            cat("Iteration ", j, "\n", append = TRUE)
            ## Run algorithm
            temp <- swMH(aList = preprocout$data$adjlist,
                         cdvec = cds,
                         cdorigvec = preprocout$data$initcds,
                         popvec = preprocout$data$popvec,
                         grouppopvec = preprocout$data$grouppopvec,
                         nsims = nsimsAdj[j],
                         eprob = eprob,
                         pct_dist_parity = preprocout$params$pctdistparity,
                         beta_sequence = preprocout$params$betaseq,
                         beta_weights = preprocout$params$betaweights,
                         ssdmat = preprocout$data$ssdmat,
                         lambda = lambda,
                         beta_population = betapop,
                         beta_compact = betacompact,
                         beta_segregation = betaseg,
                         beta_similar = betasimilar,
                         anneal_beta_population = preprocout$params$temperbetapop,
                         anneal_beta_compact = preprocout$params$temperbetacompact,
                         anneal_beta_segregation = preprocout$params$temperbetaseg,
                         anneal_beta_similar = preprocout$params$temperbetasimilar,
                         adjswap = preprocout$params$adjswaps)
            
            ## Combine data
            algout <- ecutsAppend(algout,temp)
            
            ## Get temperature
            beta <- temp$beta_sequence[nsimsAdj[j]]
            
            ## Check average MH ratio for target beta chain
            if(beta == target.beta & sum(nsimsAdj[1:j]) >= 10e3 & sum(nsimAdj[1:j]) < 10e3+freq & i == 1){
              if(mean(algout$mhprob) <= 0.2){
                stop("Target beta is too small. Please increase beta to a larger value")
              }
              else if(mean(algout$mhprob) >= 0.4){
                stop("Target beta is too large. Please decrease beta to a smaller value")
              }
            }
            
            ## Get likelihood
            if(constraint == "compact"){
                like <- exp(temp$constraint_compact[nsimsAdj[j]]) 
            }
            else if(constraint == "population"){
                like <- exp(temp$constraint_pop[nsimsAdj[j]]) 
            }
            else if(constraint == "segregation"){
                like <- exp(temp$constraint_segregation[nsimsAdj[j]]) 
            }
            else{
                like <- exp(temp$constraint_similar[nsimsAdj[j]]) 
            }
            ## Use MPI to exchange swap information
            ##
            ## Tag guide: 1 -> likelihood sent
            ##            2 -> temperature sent
            ##            3 -> acceptance probability sent
            ##
            
            if(adjswaps){

                ## Determine which nodes are swapping
                tempseg <- swaps[(i-1)*nsims + j*freq]
                ## Get node indices
                temps <- tempadj[tempseg:(tempseg+1)]
                ## Communication step        
                if(procID %in% temps){
                    ## Determine partner
                    partner <- temps[procID != temps]
                    ## Send commands (blocking)
                    Rmpi::mpi.send.Robj(like,dest=partner,tag=1)
                    Rmpi::mpi.send.Robj(beta,dest=partner,tag=2)
                    ## Receive commands (blocking)
                    likePart <- Rmpi::mpi.recv.Robj(partner,tag=1)
                    betaPart <- Rmpi::mpi.recv.Robj(partner,tag=2)
                    
                    ## Higher ranked process communicates random
                    ## draw to lower ranked process
                    if(partner < procID){
                        accept <- runif(1)
                        Rmpi::mpi.send.Robj(accept,dest=partner,tag=3)
                    }else{
                        accept <- Rmpi::mpi.recv.Robj(partner,tag=3)
                    }
                    
                    ## Compute acceptance probability (for now, population only)
                    prob <- (like^betaPart*likePart^beta)/(like^beta*likePart^betaPart)
                    if(prob > accept){
                        ## Exchange temperature values
                        beta <- betaPart
                        
                        ## Adjust temperature adjacency list
                        tempadj[tempseg:(tempseg+1)] <- tempadj[(tempseg+1):tempseg]
                        ## Send temperature adjacency list
                        if(procID == tempadj[tempseg+1]){
                            oProcs <- tempadj[!(tempadj %in% temps)]
                            for(k in 1:length(oProcs)){
                                Rmpi::mpi.send.Robj(tempadj,dest=oProcs[k],tag=4)
                            }
                        }
                    }else{
                        if(procID == tempadj[tempseg]){
                            oProcs <- tempadj[!(tempadj %in% temps)]
                            for(k in 1:length(oProcs)){
                                Rmpi::mpi.send.Robj(tempadj,dest=oProcs[k],tag=4)
                            }
                        }
                    } 
                }else{
                    tempadj <- Rmpi::mpi.recv.Robj(tempadj[tempseg],tag=4)
                }
            }else{
                if(j != length(nsimsAdj) || length(nsimsAdj) == length(tempIts)){
                    ## Swap proposed
                    ## Send commands (blocking)
                    Rmpi::mpi.send.Robj(like,dest=partner[j],tag=1)
                    Rmpi::mpi.send.Robj(beta,dest=partner[j],tag=2)
                    ## Receive commands (blocking)
                    likePart <- Rmpi::mpi.recv.Robj(partner[j],tag=1)
                    betaPart <- Rmpi::mpi.recv.Robj(partner[j],tag=2)
                    
                    ## Higher ranked process communicates random
                    ## draw to lower ranked process
                    if(partner[j] < procID){
                        accept <- runif(1)
                        Rmpi::mpi.send.Robj(accept,dest=partner[j],tag=3)
                    }
                    else{
                        accept <- Rmpi::mpi.recv.Robj(partner[j],tag=3)
                    }
                    
                    ## Compute acceptance probability (for now, population only)
                    prob <- (like^betaPart*likePart^beta)/(like^beta*likePart^betaPart)
                    if(prob > accept){
                        ## Exchange temperature values
                        beta <- betaPart
                    }
                }
            }
            ## Update inputs to swMH
            cds <- temp$partitions[,nsimsAdj[j]]
            ## End loop over j
        }
        
        class(algout) <- "redist"
        
        ## Save random number state if setting the seed
        if(!is.null(rngseed)){
            algout$randseed <- .Random.seed
        }

        ## Save output
        if(nloop > 1){
            save(algout, file = paste(savename, "_loop", i,"_temp",
                             procID, ".RData", sep = ""))
        }else if(!is.null(savename)){
            save(algout, file = paste(savename, "_temp",
                             algout$beta_sequence[nsims],
                             ".RData", sep = ""))
        }
        ## End loop over i
    }
    if(params$verbose){
        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        
        cat("redist.mcmc.mpi() simulations finished.\n", append = TRUE)
        sink()
    }
    ## End function
}

ecutsAppend <- function(algout,ndata){
    if(length(algout) == 0){
        algout <- ndata
    }else{
        algout$partitions <- cbind(algout$partitions,ndata$partitions)
        algout$distance_parity <- c(algout$distance_parity,ndata$distance_parity)
        algout$distance_original <- c(algout$distance_original, ndata$distance_original)
        algout$mhdecisions <- c(algout$mhdecisions,ndata$mhdecisions)
        algout$mhprob <- c(algout$mhprob,ndata$mhprob)
        algout$pparam <- c(algout$pparam,ndata$pparam)
        algout$constraint_pop <- c(algout$constraint_pop,ndata$constraint_pop)
        algout$constraint_compact <- c(algout$constraint_compact,ndata$constraint_compact)
        algout$constraint_segregation <- c(algout$constraint_segregation,ndata$constraint_segregation)
        algout$constraint_similar <- c(algout$constraint_similar,ndata$constraint_similar)
        algout$beta_sequence<- c(algout$beta_sequence,ndata$beta_sequence)
    }
    return(algout)
}

redist.mcmc.mpi <- function(adjobj, popvec, nsims, ndists = NA, initcds = NULL,
                            loopscompleted = 0, nloop = 1, nthin = 1,
                            eprob = 0.05,
                            lambda = 0, popcons = NA, grouppopvec = NA,
                            ssdmat = NA,rngseed = NA,
                            beta = -10, constraint = "population",  
                            betaseqlength = 10, adjswaps = TRUE,
                            freq = 100, savename = NA, maxiterrsg = 5000,
                            contiguitymap = "rooks", verbose = FALSE
){
    
    #########################
    ## Inputs to function: ##
    #########################
    ## adjobj - adjacency object of geographic units. Accepts adjlist or adjmat
    ## popvec - population of each of the units
    ## nsims - number of iterations to run the algorithm (per loop)
    ## ndists - Target number of congressional districts. Default is NULL
    ## initcds - initial congressional units. Must be contiguous partitions. Default is Null
    ## loopscompleted - number of loops completed. Default is 0
    ## nloop - Total number of loops desired across all simulations.
    ## nthin - How much to thin the chain. Defaulted to 1 (keep everything)
    ## eprob - edgecut probability. Defaulted to 0.05
    ## lambda - expected number of swaps attempted per algorithm iteration.
    ##          Default set to 0 (single swap).
    ## popcons - strength of hard population constraint. Defaulted to no
    ##           constraint. popcons = 0.01 implies a 1% population constraint.
    ## grouppopvec - vector of populations for a minority group. To be used
    ##               in conjunction with the segregation M-H constraint
    ## ssdmat - matrix of squared distances between population units.
    ##          To be used when applying the compactness constraint.
    ## beta - Constraint strength
    ## constraint - Population, compactness, segregation, or dissimilarity
    ## betaseqlength - Number of temperatures used in tempering (i.e. number of chains). Default is
    ##                 ten
    ## freq - Frequency of between-chain swaps. Default to once every 100 iterations
    ## savename - Where to save the simulations
    ## maxiterrsg - maximum number of iterations for random seed and grow starts
    ## verbose - whether to print initialization script

    ## Check if Rmpi library is installed
    if (!requireNamespace("Rmpi", quietly = TRUE)) {
        stop("You must install package 'Rmpi' to use this function. Please install it if you wish to continue."
            ,call. = FALSE)
    }
    
    ## Load Rmpi library
    if (!is.loaded("mpi_initialize")) { 
        library("Rmpi") 
    }

    ##########################
    ## Is anything missing? ##
    ##########################
    if(missing(adjobj)){
        stop("Please supply adjacency matrix or list")
    }
    if(missing(popvec)){
        stop("Please supply vector of geographic unit populations")
    }
    if(missing(nsims)){
        stop("Please supply number of simulations to run algorithm")
    }
    if(is.null(ndists) & is.null(initcds)){
        stop("Please provide either the desired number of congressional districts
         or an initial set of congressional district assignments")
    }
    if(nloop > 1 & missing(savename)){
        stop("Please supply save directory if saving simulations at checkpoints")
    }
    
    ###################
    ## Preprocessing ##
    ###################
  
    ## Augment initcds if necessary
    nrow.init <- ifelse(is.null(initcds), 0, nrow(initcds))
    ncol.init <- ifelse(is.null(initcds), ndists, ncol(initcds))
    if(nrow.init < betaseqlength){
        initcds <- rbind(initcds,matrix(NA,betaseqlength-nrow.init,ncol.init))
    }
    
    ## Generate temperature sequence (power law)
    temp <- rep(NA, betaseqlength)
    for(i in 1:betaseqlength){
        temp[i] <- 0.1^((i-1) / (betaseqlength - 1)) - .1
    }
    beta <- temp*beta/0.9
    target.beta <- beta[1]
  
    ## Generate swapping sequence
    if(adjswaps){
        swaps <- matrix(NA,1,nsims*(nloop-loopscompleted))
        ## partner <- matrix(NA,1,nits)
        for(i in 1:length(swaps)){
            if(i %% freq == 0){
                swaps[i] = sample(1:(betaseqlength-1),size=1)
            }
            ## Initial temperature adjacency
            tempadj <- 1:betaseqlength
        }
    }else{
        swaps <- matrix(NA,2,nsims*(nloop-loopscompleted))
        for(i in 1:ncol(swaps)){
            if(i %% freq == 0){
                swaps[,i] = sample(1:betaseqlength,size=2)
            }
        }
    }
    
    ## Create parameters list to distribute across nodes
    params <- expand.grid(nsims = nsims,nloop = nloop,eprob = eprob,
                          ndists = ndists,lambda = lambda,popcons = popcons,
                          beta = beta,target.beta = target.beta,constraint = constraint,
                          betaseqlength = betaseqlength,adjswaps = adjswaps,
                          nthin = nthin,freq = freq,maxiterrsg = maxiterrsg,
                          contiguitymap = contiguitymap,verbose = verbose,
                          loopscompleted = loopscompleted,rngseed = rngseed,
                          savename = savename)

    ##################
    ## Spawn Slaves ##
    ##################
    ## Note this will not work on Windows platform
    Rmpi::mpi.spawn.Rslaves(nslaves = betaseqlength)
    
    ## Get processor ID for each slave
    Rmpi::mpi.bcast.cmd(procID <- mpi.comm.rank())
    
    
    #########################
    ## Send Data to Slaves ##
    #########################
    ## Swapping Schedule
    Rmpi::mpi.bcast.Robj2slave(swaps)
    
    ## Temperature adjacency
    if(adjswaps){
        Rmpi::mpi.bcast.Robj2slave(tempadj)
    }
    
    ## Adjacency Object
    Rmpi::mpi.bcast.Robj2slave(adjobj)
    
    ## Population Vector
    Rmpi::mpi.bcast.Robj2slave(popvec)
    
    ## Initial Plans
    initcds <- split(initcds, f=1:nrow(initcds))
    Rmpi::mpi.scatter.Robj2slave(initcds)
    
    ## Group population vector
    Rmpi::mpi.bcast.Robj2slave(grouppopvec)
    
    ## Squared-distance matrix
    Rmpi::mpi.bcast.Robj2slave(ssdmat)
    
    ## Parameters List
    params <- split(params, f=1:nrow(params))
    Rmpi::mpi.scatter.Robj2slave(params)
    
    ## Send ecutsMPI function to slaves
    Rmpi::mpi.bcast.Robj2slave(ecutsMPI)
    
    ## Send ecutsAppend function to slaves
    Rmpi::mpi.bcast.Robj2slave(ecutsAppend)
    
    ## Execute ecutsMPI program on each slave
    Rmpi::mpi.bcast.cmd(ecutsMPI())
    
    ## Close slaves
    Rmpi::mpi.close.Rslaves()
    
    ## Terminate MPI processes and close R
    Rmpi::mpi.quit()
  
}
