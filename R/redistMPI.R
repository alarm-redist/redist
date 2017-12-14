##############################################
# Authors: Ben Fifield & Alex Tarr           #
# Created: 2015/06/01                        #
# Last Revision: N/A                         #
# Institution: Princeton University          #
# Purpose: R wrapper to run swMH() code w/   #
# parallel tempering (mpi)                   #
##############################################

ecutsMPI <- function(procID = procID, params = params, adjobj = adjobj, popvec = popvec, initcds = initcds, swaps = swaps){
    ## Load redist library
    library(redist)

    if(is.na(params$savename)){
        fname <- paste0("log", procID)
    }else{
        fname <- paste0("log", procID, "_", params$savename)
    }
    
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
    if(sum(is.na(initcds)) == length(initcds)){
        initcds <- NULL
    }

    nthin <- params$nthin
    
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
            algout <- list()
                        
        } else{ 
            
            ## Reload the data if restarting
            if(loopstart > 1){
                
                ## Load the data
                load(paste(savename,"_proc", procID, "_loop", i - 1, ".RData", sep = ""))
                
                ## Load the temperature adjacency matrix (need to specify WD)
                load(paste0(savename,"_tempadjMat.RData"))
                tempadj <- tempadjMat[nrow(tempadjMat),]
                
                ## Load the swapping schedule (need to specify WD)
                load(paste0(savename,"_swaps.RData"))
                
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
                algout <- list()
                                
            }else{
                cds <- preprocout$data$initcds
                ## Initialize algout object (for use in ecutsAppend)
                algout <- list()
                ## Temperature Adjacency Matrix
                if(adjswaps){
                  tempadjMat <- tempadj
                }
            }
            
        } 
        
        #######################
        ## Run the algorithm ##
        #######################
        for(j in 1:length(nsimsAdj)){

            cat("Swap ", j, "\n", append = TRUE)
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
            cat("Likelihood is", like, "\n", append = TRUE)
            ## Use MPI to exchange swap information
            ##
            ## Tag guide: 1 -> likelihood sent
            ##            2 -> temperature sent
            ##            3 -> acceptance probability sent
            ##

            if(adjswaps){

                cat("Start adjswaps\n", append = TRUE)
                ## Determine which nodes are swapping
                tempseg <- swaps[(i-1)*nsims + j*freq]
                ## Get node indices
                temps <- tempadj[tempseg:(tempseg+1)]
                ## Communication step
                cat("Is procID in temps? If true, enter communication",
                    procID %in% temps, "\n", append = TRUE)
                if(procID %in% temps){
                    cat("Enter communication step\n", append = TRUE)
                    ## Determine partner
                    cat("Start determine partner\n", append = TRUE)
                    partner <- temps[procID != temps]
                    ## Send commands (blocking)
                    cat("Start send commands\n", append = TRUE)
                    Rmpi::mpi.send.Robj(like,dest=partner,tag=1)
                    Rmpi::mpi.send.Robj(beta,dest=partner,tag=2)
                    ## Receive commands (blocking)
                    cat("Start receive commands\n", append = TRUE)
                    likePart <- Rmpi::mpi.recv.Robj(partner,tag=1)
                    betaPart <- Rmpi::mpi.recv.Robj(partner,tag=2)
                    
                    ## Higher ranked process communicates random
                    ## draw to lower ranked process
                    cat("Start communicate draw\n", append = TRUE)
                    if(partner < procID){
                        accept <- runif(1)
                        Rmpi::mpi.send.Robj(accept,dest=partner,tag=3)
                    }else{
                        accept <- Rmpi::mpi.recv.Robj(partner,tag=3)
                    }
                    
                    ## Compute acceptance probability (for now, population only)
                    cat("Start compute acceptance prob\n", append = TRUE)
                    
                    cat("Components of likelihood:\n", append = TRUE)
                    cat("like =", like, "\n", append = TRUE)
                    cat("beta =", beta, "\n", append = TRUE)
                    cat("betaPart =", betaPart, "\n", append = TRUE)
                    cat("likePart =", likePart, "\n", append = TRUE)
                    cat("Constraint in likelihood =", log(like), "\n", append = TRUE)
                    cat("Constraint in likelihoodPart =", log(likePart), "\n", append = TRUE)

                    prob <- (like^betaPart*likePart^beta)/(like^beta*likePart^betaPart)

                    cat("Prob =", prob, "and accept =", accept, "\n", append = TRUE)
                    if(prob > accept){
                        cat("Prob > accept\n", append = TRUE)
                        ## Exchange temperature values
                        cat("Start exchange temps\n", append = TRUE)
                        beta <- betaPart
                        
                        ## Adjust temperature adjacency list
                        cat("Start adjust tempAdj list\n", append = TRUE)
                        tempadj[tempseg:(tempseg+1)] <- tempadj[(tempseg+1):tempseg]
                        ## Send temperature adjacency list
                        cat("Start send tempAdj\n", append = TRUE)
                        if(procID == tempadj[tempseg+1]){
                            oProcs <- tempadj[!(tempadj %in% temps)]
                            for(k in 1:length(oProcs)){
                                Rmpi::mpi.send.Robj(tempadj,dest=oProcs[k],tag=4)
                            }
                        }
                    }else{
                        cat("Prob < accept\n", append = TRUE)
                        if(procID == tempadj[tempseg]){
                            oProcs <- tempadj[!(tempadj %in% temps)]
                            for(k in 1:length(oProcs)){
                                Rmpi::mpi.send.Robj(tempadj,dest=oProcs[k],tag=4)
                            }
                        }
                    }
                    cat("Exit communication step\n", append = TRUE)
                }else{
                    cat("Not in temps, receive tempadj from other nodes\n",
                        append = TRUE)
                    tempadj <- Rmpi::mpi.recv.Robj(tempadj[tempseg],tag=4)
                    cat("End receive tempadj from other nodes\n", append = TRUE)
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
            
            ## Update tempadjMat
            if(adjswaps){
              ## Update temperature adjacency matrix
              tempadjMat <- rbind(tempadjMat,tempadj)
            }
            ## End loop over j
        }
        
        class(algout) <- "redist"
        
        ## Save random number state if setting the seed
        if(!is.null(rngseed)){
            algout$randseed <- .Random.seed
        }

        ## Save output
        if(nloop > 1){
            save(algout, file = paste(savename, "_proc", procID,"_loop", 
                                      i, ".RData", sep = ""))
          ## Save temperature adjacency matrix
          if(adjswaps){
              save(tempadjMat, file = paste0(savename, "_tempadjMat.RData"))
          }
          
          ## Save swaps
          save(swaps,file = paste0(savename, "_swaps.RData"))
          
        }else if(!is.null(savename)){
            save(algout, file = paste(savename, "_chain",
                             algout$beta_sequence[nsims],
                             ".RData", sep = ""))
        }
        ## End loop over i
    }
    
    ###############################
    ## Combine and save the data ##
    ###############################
    if(nloop > 1){
      if(procID == tempadj[1]){
        redist.combine.mpi(savename = savename, nsims = nsims, nloop = nloop,
                           nthin = nthin, nunits = length(preprocout$data$adjlist),
                           tempadj = tempadj)
      }
    }else if(!is.null(savename)){
      save(algout, file = paste(savename, ".RData", sep = ""))
  }
    
    if(params$verbose){
        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        
        cat("redist.mcmc.mpi() simulations finished.\n", append = TRUE)
        sink()
    }
    ## End function
}

#' Combine successive runs of \code{redist.mcmc.mpi}
#'
#' \code{redist.combine.mpi} is used to combine successive runs of
#' \code{redist.mcmc.mpi} into a single data object
#'
#' @usage redist.combine.mpi(savename, nsims, nloop, nthin, nunits, tempadj)
#'
#' @param savename The name (without the loop or \code{.RData} suffix)
#' of the saved simulations.
#' @param nsims The number of simulations in each loop.
#' @param nloop The number of loops being combined.
#' @param nthin How much to thin the simulations being combined.
#' @param nunits The number of geographic units from the simulations.
#' @param tempadj The temperature adjacency object saved by
#' \code{redist.mcmc.mpi}.
#'
#' @details This function allows users to combine multiple successive runs of
#' \code{redist.mcmc.mpi} into a single \code{redist} object for analysis.
#'
#' @return \code{redist.combine.mpi} returns an object of class "redist".
#' The object \code{redist} is a list that contains the folowing components (the
#' inclusion of some components is dependent on whether tempering
#' techniques are used):
#' \item{partitions}{Matrix of congressional district assignments generated by the
#' algorithm. Each row corresponds to a geographic unit, and each column
#' corresponds to a simulation.}
#' \item{distance_parity}{Vector containing the maximum distance from parity for
#' a particular simulated redistricting plan.}
#' \item{mhdecisions}{A vector specifying whether a proposed redistricting plan
#' was accepted (1) or rejected (0) in a given iteration.}
#' \item{mhprob}{A vector containing the Metropolis-Hastings acceptance
#' probability for each iteration of the algorithm.}
#' \item{pparam}{A vector containing the draw of the \code{p} parameter for each
#' simulation, which dictates the number of swaps attempted.}
#' \item{constraint_pop}{A vector containing the value of the population
#' constraint for each accepted redistricting plan.}
#' \item{constraint_compact}{A vector containing the value of the compactness
#' constraint for each accepted redistricting plan.}
#' \item{constraint_segregation}{A vector containing the value of the
#' segregation constraint for each accepted redistricting plan.}
#' \item{constraint_similar}{A vector containing the value of the similarity
#' constraint for each accepted redistricting plan.}
#' \item{beta_sequence}{A vector containing the value of beta for each iteration
#' of the algorithm. Returned when tempering is being used.}
#' \item{mhdecisions_beta}{A vector specifying whether a proposed beta value was
#' accepted (1) or rejected (0) in a given iteration of the algorithm. Returned
#' when tempering is being used.}
#' \item{mhprob_beta}{A vector containing the Metropolis-Hastings acceptance
#' probability for each iteration of the algorithm. Returned when tempering
#' is being used.}
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' @examples
#' \dontrun{
#' data(algdat.pfull)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins, Imai and
#' ## Tarr (2015)
#'
#' ## Get an initial partition
#' set.seed(1)
#' initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]
#'
#' ## Run the algorithm
#' redist.mcmc.mpi(adjobj = algdat.pfull$adjlist,
#' popvec = algdat.pfull$precinct.data$pop,
#' initcds = initcds,
#' nsims = 10000, nloops = 2, savename = "test")
#' out <- redist.combine.mpi(savename = "test", nsims = 10000, nloop = 2,
#' nthin = 10, nunits = length(algdat.pfull$adjlist), tempadj = tempAdjMat)
#' }
#' @export
redist.combine.mpi <- function(savename, nsims, nloop, nthin, nunits, tempadj){
    
    ##############################
    ## Set up container objects ##
    ##############################
    partitions <- matrix(NA, nrow = nunits,
                         ncol = (nsims * nloop / nthin))
    
    distance_parity <- rep(NA, (nsims * nloop / nthin))
    distance_original <- rep(NA, (nsims * nloop / nthin))
    mhdecisions <- rep(NA, (nsims * nloop / nthin))
    mhprob <- rep(NA, (nsims * nloop / nthin))
    pparam <- rep(NA, (nsims * nloop / nthin))
    constraint_pop <- rep(NA, (nsims * nloop / nthin))
    constraint_compact <- rep(NA, (nsims * nloop / nthin))
    constraint_segregation <- rep(NA, (nsims * nloop / nthin))
    constraint_similar <- rep(NA, (nsims * nloop / nthin))
    
    beta_sequence <- rep(NA, (nsims * nloop / nthin))
    
    ## Indices for thinning
    indthin <- which((1:nsims) %% nthin == 0)
    
    ####################################
    ## Combine data in multiple loops ##
    ####################################

    for(i in 1:nloop){

        ## Load data
        load(paste(savename, "_proc", tempadj[1], "_loop", i, ".RData", sep = ""))

        ind <- ((i - 1) * (nsims / nthin) + 1):(i * (nsims / nthin))
        
        ## Store objects together
        partitions[1:nunits, ind] <- algout$partitions[,indthin]
        
        distance_parity[ind] <- algout$distance_parity[indthin]
        distance_original[ind] <- algout$distance_original[indthin]
        mhdecisions[ind] <- algout$mhdecisions[indthin]
        mhprob[ind] <- algout$mhprob[indthin]
        pparam[ind] <- algout$pparam[indthin]
        constraint_pop[ind] <- algout$constraint_pop[indthin]
        constraint_compact[ind] <- algout$constraint_compact[indthin]
        constraint_segregation[ind] <- algout$constraint_segregation[indthin]
        constraint_similar[ind] <- algout$constraint_similar[indthin]
        
        beta_sequence[ind] <- algout$beta_sequence[indthin]
        
    }
    
    #################################
    ## Store data in algout object ##
    #################################

    algout <- vector(mode = "list")
    
    algout$partitions <- partitions
    algout$distance_parity <- distance_parity
    algout$distance_original <- distance_original
    algout$mhdecisions <- mhdecisions
    algout$mhprob <- mhprob
    algout$pparam <- pparam
    algout$constraint_pop <- constraint_pop
    algout$constraint_compact <- constraint_compact
    algout$constraint_segregation <- constraint_segregation
    algout$constraint_similar <- constraint_similar

    algout$beta_sequence <- beta_sequence

    
    #########################
    ## Set class of object ##
    #########################
    class(algout) <- "redist"
    
    #################
    ## Save object ##
    #################
    save(algout, file = paste(savename, ".RData", sep = ""))
    return(algout)
    
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

#' MCMC Redistricting Simulator using MPI
#'
#' \code{redist.mcmc.mpi} is used to simulate Congressional redistricting
#' plans using Markov Chain Monte Carlo methods.
#'
#' @usage redist.mcmc.mpi(adjobj, popvec, nsims, ndists = NA, initcds = NULL,
#' loopscompleted = 0, nloop = 1, nthin = 1,
#' eprob = 0.05,
#' lambda = 0, popcons = NA, grouppopvec = NA,
#' ssdmat = NA, rngseed = NA,
#' beta = -10, constraint = "population",
#' betaseqlength = 10, adjswaps = TRUE,
#' freq = 100, savename = NA, maxiterrsg = 5000,
#' contiguitymap = "rooks", verbose = FALSE)
#'
#' @param adjobj An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param popvec A vector containing the populations of each geographic
#' unit.
#' @param nsims The number of simulations run before a save point.
#' @param ndists The number of congressional districts. The default is
#' \code{NULL}.
#' @param initcds A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided, random
#' and contiguous congressional district assignments will be generated using
#' \code{redist.rsg}.
#' @param loopscompleted Number of save points reached by the
#' algorithm. The default is \code{0}.
#' @param nloop The total number of save points for the algorithm. The
#' default is \code{1}. Note that the total number of simulations run
#' will be \code{nsims} * \code{nloop}.
#' @param nthin The amount by which to thin the Markov Chain. The default
#' is \code{1}.
#' @param eprob The probability of keeping an edge connected. The default
#' is \code{0.05}.
#' @param lambda The parameter determining the number of swaps to attempt
#' each iteration of the algorithm. The number of swaps each iteration is
#' equal to Pois(\code{lambda}) + 1. The default is \code{0}.
#' @param popcons The strength of the hard population
#' constraint. \code{popcons} = 0.05 means that any proposed swap that
#' brings a district more than 5\% away from population parity will be
#' rejected. The default is \code{NULL}.
#' @param grouppopvec A vector of populations for some sub-group of
#' interest. The default is \code{NULL}.
#' @param ssdmat A matrix of squared distances between geographic
#' units. The default is \code{NULL}.
#' @param rngseed Allows the user to set the seed for the
#' simulations. Default is \code{NULL}.
#' @param beta The strength of the target strength in the MH ratio. The default
#' is 0.
#' @param constraint Which constraint to apply. Accepts \code{compact},
#' \code{segregation}, \code{population}, \code{similarity}, or \code{none}
#' (no constraint applied). The default is \code{none}.
#' @param betaseqlength Length of beta sequence desired for
#' tempering. The default is \code{10}.
#' @param adjswaps Flag to restrict swaps of beta so that only
#' values adjacent to current constraint are proposed. The default is
#' \code{TRUE}.
#' @param freq Frequency of between-chain swaps. Default to once every 100
#' iterations
#' @param savename Filename to save simulations. Default is \code{NULL}.
#' @param maxiterrsg Maximum number of iterations for random seed-and-grow
#' algorithm to generate starting values. Default is 5000.
#' @param contiguitymap Use queens or rooks distance criteria for generating an
#' adjacency list from a "SpatialPolygonsDataFrame" data type. Default is "rooks". 
#' @param verbose Whether to print initialization statement. Default is
#' \code{TRUE}.
#'
#' @details This function allows users to simulate redistricting plans
#' using Markov Chain Monte Carlo methods. Several constraints
#' correspoding to substantive requirements in the redistricting process
#' are implemented, including population parity and geographic
#' compactness. In addition, the function includes multiple-swap and
#' parallel tempering functionality in MPI to improve the mixing of the Markov
#' Chain.
#'
#' @return \code{redist.mcmc.mpi} returns an object of class "redist". The object
#' \code{redist} is a list that contains the folowing components (the
#' inclusion of some components is dependent on whether tempering
#' techniques are used):
#' \item{partitions}{Matrix of congressional district assignments generated by the
#' algorithm. Each row corresponds to a geographic unit, and each column
#' corresponds to a simulation.}
#' \item{distance_parity}{Vector containing the maximum distance from parity for
#' a particular simulated redistricting plan.}
#' \item{mhdecisions}{A vector specifying whether a proposed redistricting plan
#' was accepted (1) or rejected (0) in a given iteration.}
#' \item{mhprob}{A vector containing the Metropolis-Hastings acceptance
#' probability for each iteration of the algorithm.}
#' \item{pparam}{A vector containing the draw of the \code{p} parameter for each
#' simulation, which dictates the number of swaps attempted.}
#' \item{constraint_pop}{A vector containing the value of the population
#' constraint for each accepted redistricting plan.}
#' \item{constraint_compact}{A vector containing the value of the compactness
#' constraint for each accepted redistricting plan.}
#' \item{constraint_segregation}{A vector containing the value of the
#' segregation constraint for each accepted redistricting plan.}
#' \item{constraint_similar}{A vector containing the value of the similarity
#' constraint for each accepted redistricting plan.}
#' \item{beta_sequence}{A vector containing the value of beta for each iteration
#' of the algorithm. Returned when tempering is being used.}
#' \item{mhdecisions_beta}{A vector specifying whether a proposed beta value was
#' accepted (1) or rejected (0) in a given iteration of the algorithm. Returned
#' when tempering is being used.}
#' \item{mhprob_beta}{A vector containing the Metropolis-Hastings acceptance
#' probability for each iteration of the algorithm. Returned when tempering
#' is being used.}
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' @examples
#' \dontrun{
#' data(algdat.pfull)
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins,
#' ## Imai and Tarr (2015)
#'
#' ## Get an initial partition
#' set.seed(1)
#' initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]
#'
#' ## Run the algorithm
#' redist.mcmc.mpi(adjobj = algdat.pfull$adjlist,
#' popvec = algdat.pfull$precinct.data$pop,
#' initcds = initcds,
#' nsims = 10000, savename = "test")
#' }
#' @export
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

    ## Check if Rmpi library is installed
    if (!requireNamespace("Rmpi", quietly = TRUE)) {
        stop("You must install package 'Rmpi' to use this function. Please install it if you wish to continue."
            ,call. = FALSE)
    }
    
    ## ## Load Rmpi library
    ## if (!is.loaded("mpi_initialize")) { 
    ##     library("Rmpi") 
    ## }

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
    Rmpi::mpi.bcast.cmd(procID <- Rmpi::mpi.comm.rank())
    
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
    Rmpi::mpi.bcast.cmd(ecutsMPI(procID, params, adjobj, popvec, initcds, swaps))
    
    ## Close slaves
    Rmpi::mpi.close.Rslaves()
    
    ## Terminate MPI processes and close R
    Rmpi::mpi.quit()
  
}
