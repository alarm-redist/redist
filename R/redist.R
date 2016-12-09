###########################################
## Author: Ben Fifield
## Institution: Princeton University
## Date Created: 2015/02/04
## Date Modified: 2015/03/09
## Purpose: R wrapper to run swMH() code (non-mpi)
###########################################

redist.preproc <- function(adjobj, popvec, initcds = NULL, ndists = NULL,
                           popcons = NULL, grouppopvec = NULL, ssdmat = NULL,
                           beta = 0, temper = NULL, constraint = NULL,
                           betaseq = NULL, betaseqlength = NULL,
                           betaweights = NULL, adjswaps = TRUE, maxiterrsg = NULL,
                           contiguitymap = NULL
                           ){

    #########################
    ## Inputs to function: ##
    #########################
    ## adjobj - adjacency object of geographic units. Accepts adjlist or adjmat
    ## popvec - population of each of the units
    ## initcds - initial congressional units. Must be contiguous partitions. Default is NULL
    ## ndists - number of desired congressional units. Default is NULL
    ## popcons - strength of hard population constraint. Defaulted to no
    ##           constraint. popcons = 0.01 implies a 1% population constraint.
    ## grouppopvec - vector of populations for a minority group. To be used
    ##               in conjunction with the segregation M-H constraint
    ## ssdmat - matrix of squared distances between population units.
    ##          To be used when applying the compactness constraint.
    ## beta - target strength of constraint in MH ratio. Defaults to 0.
    ## temper - whether to use tempering (parallel or simulated) algorithms.
    ##          Defaults to `none` (no tempering)
    ## constraint - which constraint to apply. Defaults to `none` (no tempering)
    ## betaseq - Spacing for beta sequence if tempering. Default is power law
    ##           spacing, but can also be provided by user
    ## betaseqlength - Number of temperatures in the beta sequence. Default is
    ##                 ten
    ## betaweights - Vector of weights for beta sequence. Provided by user
    ## adjswaps - Flag for adjacent swaps for geyer-thompson tempering or MPI
    ##            parallel tempering. Default to TRUE
    ## maxiterrsg - Maximum number of iterations for RSG algorithm
    ## contiguitymap - Distance criteria for adjacency list from input map
    
    #######################
    ## Check missingness ##
    #######################
    if(missing(adjobj)){
        stop("Please supply adjacency matrix or list")
    }
    if(missing(popvec)){
        stop("Please supply vector of geographic unit populations")
    }
    if((beta == 0 & temper == "none" & constraint != "none") |
       (beta == 0 & temper == "simulated") |
       (temper != "none" & constraint == "none")){
        stop("If applying constraints or using simulated tempering, please set non-zero constraint and specify the constraint.")
    }
    if(!(temper %in% c("none", "simulated", "parallel"))){
        stop("Please specify either `none`, `simulated` or `parallel` for tempering argument")
    }
    if(!(constraint %in% c("none", "compact", "segregation",
                           "population", "similarity"))){
        stop("Please specify `none`, `compact`, `segregation`, `population`, or `similarity` for constraint")
    }
    
    ############################################
    ## If not a list, convert adjlist to list ##
    ############################################
    if(!is.list(adjobj)){

        ## If a matrix, check to see if adjacency matrix
        if(is.matrix(adjobj)){

            ## Is it square?
            squaremat <- (nrow(adjobj) == ncol(adjobj))
            ## All binary entries?
            binary <- ((length(unique(c(adjobj))) == 2) &
                           (sum(unique(c(adjobj)) %in% c(0, 1)) == 2))
            ## Diagonal elements all 1?
            diag <- (sum(diag(adjobj)) == nrow(adjobj))
            ## Symmetric?
            symmetric <- isSymmetric(adjobj)

            ## If all are true, change to adjlist and automatically zero-index
            if(squaremat & binary & diag & symmetric){
                
                ## Initialize object
                adjlist <- vector("list", nrow(adjobj))

                ## Loop through rows in matrix
                for(i in 1:nrow(adjobj)){

                    ## Extract row
                    adjvec <- adjobj[,i]
                    ## Find elements it is adjacent to
                    inds <- which(adjobj == 1)
                    ## Remove self-adjacency
                    inds <- inds[inds != i,]
                    ## Zero-index
                    inds <- inds - 1
                    ## Put in adjlist
                    adjlist[[i]] <- inds
                    
                }
                
            }else { ## If not valid adjacency matrix, throw error
                stop("Please input valid adjacency matrix")
            }
        }else if(class(adjobj) == "SpatialPolygonsDataFrame"){ ## shp object

            ## Distance criterion
            queens <- ifelse(contiguitymap == "rooks", FALSE, TRUE)
            
            ## Convert shp object to adjacency list
            adjlist <- poly2nb(adjobj, queen = queens)
            
            ## Zero-index list
            for(i in 1:length(adjlist)){
                adjlist[[i]] <- adjlist[[i]] - 1
            }
            
            ## Change class to list
            class(adjlist) <- "list"
            
        }else{ ## If neither list, matrix, or shp, throw error
            stop("Please input an adjacency list, adjacency matrix, or Spatial
                 Polygons shp file")
        }

    }else{

        ## Rename adjacency object as list
        adjlist <- adjobj

        ## Is list zero-indexed?
        minlist <- min(unlist(adjlist))
        maxlist <- max(unlist(adjlist))
        oneind <- (sum(minlist == 1, maxlist == length(adjlist)) == 2)
        zeroind <- (sum(minlist == 0, maxlist == (length(adjlist) - 1)) == 2)
        
        if(oneind){
            ## Zero-index list
            for(i in 1:length(adjlist)){
                adjlist[[i]] <- adjlist[[i]] - 1
            }
        }else if(!(oneind | zeroind)){
            ## if neither oneind or zeroind, then stop
            stop("Adjacency list must be one-indexed or zero-indexed")
        }
        
    }

    ###################################################################
    ## Check whether initial partitions (if provided) are contiguous ##
    ###################################################################
    if(!is.null(initcds)){
        if(!is.na(initcds)[1]){
            ndists <- length(unique(initcds))
            divlist <- genAlConn(adjlist, initcds)
            ncontig <- countpartitions(divlist)
            
            if(ncontig != ndists){
                stop(paste("Your initial congressional districts have ", ndists,
                           " unique districts but ",
                           ncontig, " contigous connected components. Please provide a starting map with contigous districts.", sep = ""))
            }
        }
    }

    ##############################################################################
    ## If no initial congressional districts provided, use Random Seed and Grow ##
    ## (Chen and Rodden 2013) algorithm                                         ##
    ##############################################################################
    if(is.null(initcds) | is.na(initcds)[1]){

        ## Set up target pop, strength of constraint (10%)
        if(is.null(popcons)){
            popcons_rsg <- .1
        }else{
            popcons_rsg <- popcons
        }

        ## Print start
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")
        
        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        cat("redist.rsg(): Automated Redistricting Starts\n\n", append= TRUE)
        
        ## Run the algorithm
        initout <- redist.rsg(adj.list = adjlist,
                              population = popvec,
                              ndists = ndists,
                              thresh = popcons_rsg,
                              verbose = FALSE,
                              maxiter = maxiterrsg)
        ## Get initial cds
        initcds <- initout$district_membership
        
    }
    ###########################################################
    ## Check other inputs to make sure they are right length ##
    ###########################################################
    if((length(popvec) != length(adjlist)) | (sum(is.na(popvec)) > 0)){
        stop("Each entry in adjacency list must have a corresponding entry
              in vector of populations")
    }
    if((length(initcds) != length(adjlist)) | (sum(is.na(initcds)) > 0)){
        stop("Each entry in adjacency list must have an initial congressional
             district assignment")
    }
    if(constraint == "segregation" & is.null(grouppopvec)){
        stop("If applying the segregation constraint, please provide a vector
             of subgroup populations")
    }
    if(constraint == "segregation" & !(is.null(grouppopvec))){
        if((length(grouppopvec) != length(adjlist)) |
           (sum(is.na(grouppopvec)) > 0)){
            stop("If applying the segregation constraint, each entry in adjacency
              list must have corresponding entry in vector of group populations")
        }
    }

    ####################
    ## Zero-index cds ##
    ####################
    if(min(initcds) != 0){
        initcds <- initcds - min(initcds)
    }
    if(length(unique(initcds)) != (max(initcds) + 1)){
        stop("Need congressional assignment ids to be sequence increasing by 1")
    }

    ####################################################
    ## Calculate parity and population margin allowed ##
    ####################################################
    dists <- length(unique(initcds))
    if(is.null(popcons)){
        popcons <- 100
    }

    #####################################
    ## Set grouppopvec if not provided ##
    #####################################
    if(is.null(grouppopvec)){
        grouppopvec <- popvec
    }

    ################################
    ## Set ssdmat if not provided ##
    ################################
    if(is.null(ssdmat) & constraint == "compact"){
        if(class(adjobj) == "SpatialPolygonsDataFrame"){
            centroids <- coordinates(adjobj)
            ssdmat <- calcPWDh(centroids)
        }else{
            stop("Provide squared distances matrix if constraining compactness")
        }
    }else if(is.null(ssdmat)){
        ssdmat <- matrix(1, length(adjlist), length(adjlist))
    }

    ########################
    ## Set up constraints ##
    ########################
    betapop <- 0; betacompact <- 0; betaseg <- 0; betasimilar <- 0
    temperbetapop <- 0; temperbetacompact <- 0
    temperbetaseg <- 0; temperbetasimilar <- 0
    if(constraint == "compact"){
        betacompact <- beta
        if(temper == "simulated"){
            temperbetacompact <- 1
        }
    }
    if(constraint == "segregation"){
        betaseg <- beta
        if(temper == "simulated"){
            temperbetaseg <- 1
        }
    }
    if(constraint == "population"){
        betapop <- beta
        if(temper == "simulated"){
            temperbetapop <- 1
        }
    }
    if(constraint == "similarity"){
        betasimilar <- beta
        if(temper == "simulated"){
            temperbetasimilar <- 1
        }
    }
    
    ###################################
    ## Check if betaspacing provided ##
    ###################################
    if(temperbetacompact == 1 | temperbetapop == 1 |
       temperbetaseg == 1 | temperbetasimilar == 1){
        if(betaseq[1] == "powerlaw"){

            ## Stop if no target provided
            if(betaseq[1] == "powerlaw" & betapop == 0 & temperbetapop == 1){
                stop("Provide target beta value for population constraint")
            }
            if(betaseq[1] == "powerlaw" & betacompact == 0 & temperbetacompact == 1){
                stop("Provide target beta value for compactness constraint")
            }
            if(betaseq[1] == "powerlaw" & betaseg == 0 & temperbetaseg == 1){
                stop("Provide target beta value for segregation constraint")
            }
            if(betaseq[1] == "powerlaw" & betasimilar == 0 & temperbetasimilar == 1){
                stop("Provide target beta value for similarity constraint")
            }

            ## Generate power law sequence
            betaseq <- rep(NA, betaseqlength)
            for(i in 1:length(betaseq)){
                betaseq[i] <- 0.1^((i-1) / (length(betaseq) - 1)) - .1
            }

            ## Get multiplicative constant to get desired sequence
            if(temperbetacompact == 1){
                multip <- betacompact / .9
                betacompact <- 0
            }
            if(temperbetapop == 1){
                multip <- betapop / .9
                betapop <- 0
            }
            if(temperbetaseg == 1){
                multip <- betaseg / .9
                betaseg <- 0
            }
            if(temperbetasimilar == 1){
                multip <- betasimilar / .9
                betasimilar <- 0
            }

            ## Multiply the sequence by constant
            betaseq <- betaseq * multip
            
        }else if(is.vector(betaseq)){
            betaseq <- betaseq
        }else if(!is.vector(betaseq) & betaseq[1] != "powerlaw"){
            stop("Please provide valid sequence of betas")
        }
        if(is.null(betaweights)){
            betaweights <- rep(1, length(betaseq))
        }
    }else{
        betaseq <- c(1, 1, 1, 1)
        betaweights <- c(1, 1, 1, 1)
    }

    ## Reverse beta sequence
    betaseq <- rev(betaseq)

    ########################################
    ## Convert adjacent swaps flag to 0/1 ##
    ########################################
    adjswaps <- adjswaps * 1

    #################
    ## Return list ##
    #################
    preprocout <- list(data = list(adjlist = adjlist,
                           popvec = popvec,
                           initcds = initcds,
                           grouppopvec = grouppopvec,
                           ssdmat = ssdmat
                                   ),
                       params = list(pctdistparity = popcons,
                           dists = dists,
                           betacompact = betacompact,
                           betapop = betapop,
                           betaseg = betaseg,
                           betasimilar = betasimilar,
                           temperbetacompact = temperbetacompact,
                           temperbetapop = temperbetapop,
                           temperbetaseg = temperbetaseg,
                           temperbetasimilar = temperbetasimilar,
                           betaseq = betaseq,
                           betaweights = betaweights,
                           adjswaps = adjswaps
                                     )
                       )

    class(preprocout) <- "redist"
    
    return(preprocout)
    
}

redist.combine <- function(savename, nsims, nloop, nthin, nunits, temper = 0
                           ){

    ########################
    ## Inputs to function ##
    ########################
    ## savename - filename
    ## nsims - number of simulations in each loop
    ## nloop - number of loops to run simulations
    ## nthin - how much to thin the simulations
    ## nunits - number of geographic units
    ## temper - whether tempering is used in the simulations

    ##############################
    ## Set up container objects ##
    ##############################
    partitions <- matrix(NA, nrow = nunits,
                         ncol = (nsims * nloop / nthin))

    distance_parity <- rep(NA, (nsims * nloop / nthin))
    mhdecisions <- rep(NA, (nsims * nloop / nthin))
    mhprob <- rep(NA, (nsims * nloop / nthin))
    pparam <- rep(NA, (nsims * nloop / nthin))
    constraint_pop <- rep(NA, (nsims * nloop / nthin))
    constraint_compact <- rep(NA, (nsims * nloop / nthin))
    constraint_segregation <- rep(NA, (nsims * nloop / nthin))
    constraint_similar <- rep(NA, (nsims * nloop / nthin))

    if(temper == 1){
        beta_sequence <- rep(NA, (nsims * nloop / nthin))
        mhdecisions_beta <- rep(NA, (nsims * nloop / nthin))
        mhprob_beta <- rep(NA, (nsims * nloop / nthin))
    }
    
    ## Indices for thinning
    indthin <- which((1:nsims) %% nthin == 0)

    ####################################
    ## Combine data in multiple loops ##
    ####################################
    for(i in 1:nloop){

        ## Load data
        load(paste(savename, "_loop", i, ".RData", sep = ""))

        ind <- ((i - 1) * (nsims / nthin) + 1):(i * (nsims / nthin))
    
        ## Store objects together
        partitions[1:nunits, ind] <- algout$partitions[,indthin]

        distance_parity[ind] <- algout$distance_parity[indthin]
        mhdecisions[ind] <- algout$mhdecisions[indthin]
        mhprob[ind] <- algout$mhprob[indthin]
        pparam[ind] <- algout$pparam[indthin]
        constraint_pop[ind] <- algout$constraint_pop[indthin]
        constraint_compact[ind] <- algout$constraint_compact[indthin]
        constraint_segregation[ind] <- algout$constraint_segregation[indthin]
        constraint_similar[ind] <- algout$constraint_similar[indthin]

        if(temper == 1){
            beta_sequence[ind] <- algout$beta_sequence[indthin]
            mhdecisions_beta[ind] <- algout$mhdecisions_beta[indthin]
            mhprob_beta[ind] <- algout$mhprob_beta[indthin]
        }
        
    }

    #################################
    ## Store data in algout object ##
    #################################
    if(temper == 1){
        algout <- vector(mode = "list")
    }else{
        algout <- vector(mode = "list")
    }
    algout$partitions <- partitions
    algout$distance_parity <- distance_parity
    algout$mhdecisions <- mhdecisions
    algout$mhprob <- mhprob
    algout$pparam <- pparam
    algout$constraint_pop <- constraint_pop
    algout$constraint_compact <- constraint_compact
    algout$constraint_segregation <- constraint_segregation
    algout$constraint_similar <- constraint_similar
    if(temper == 1){
        algout$beta_sequence <- beta_sequence
        algout$mhdecisions_beta <- mhdecisions_beta
        algout$mhprob_beta <- mhprob_beta
    }
    
    #########################
    ## Set class of object ##
    #########################
    class(algout) <- "redist"

    #################
    ## Save object ##
    #################
    save(algout, file = paste(savename, ".RData", sep = ""))
    
}

redist.mcmc <- function(adjobj, popvec, nsims, ndists = NULL, initcds = NULL,
                        loopscompleted = 0, nloop = 1, nthin = 1, eprob = 0.05,
                        lambda = 0, popcons = NULL, grouppopvec = NULL,
                        ssdmat = NULL, beta = 0, temper = "none",
                        constraint = "none",
                        betaseq = "powerlaw", betaseqlength = 10,
                        betaweights = NULL, 
                        adjswaps = TRUE, rngseed = NULL, maxiterrsg = 5000,
                        contiguitymap = "rooks", exact_mh = 0, savename = NULL, verbose = TRUE
                        ){

    #########################
    ## Inputs to function: ##
    #########################
    ## adjobj - adjacency object of geographic units. Accepts adjlist or adjmat
    ## popvec - population of each of the units
    ## nsims - number of iterations to run the algorithm
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
    ## beta - target strength of constraint in MH ratio. Defaults to 0.
    ## temper - whether to use tempering (parallel or simulated) algorithms.
    ##          Defaults to `none` (no tempering)
    ## constraint - which constraint to apply. Defaults to `none` (no tempering)
    ## betaseq - Spacing for beta sequence if tempering. Default is power law
    ## betaseqlength - Number of temperatures in the beta sequence. Default is
    ##                 ten
    ## betaweights - Vector of weights for beta sequence. Provided by user
    ## adjswaps - Flag for adjacent swaps for geyer-thompson tempering or MPI
    ##            parallel tempering. Default to TRUE
    ## rngseed - Random number generator seed number. Provided by user
    ## maxiterrsg - Maximum number of iterations for random seed-and-grow algorithm
    ## contiguitymap - Use queens distance or rooks distance for creating an adjlist
    ##                 from map. Defaults to rooks.
    ## exact_mh - Use exact metropolis hastings with boundary correction. Default to 0
    ## savename - Where to save the simulations
    ## verbose - whether to print initialization script

    if(verbose){
        ## Initialize ##
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")
        
        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        cat("redist.mcmc(): Automated Redistricting Simulation Using
         Markov Chain Monte Carlo\n\n", append = TRUE)
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
    if(!(contiguitymap %in% c("queens", "rooks"))){
        stop("Please supply `queens` or `rooks` for a distance criteria")
    }

    ## Set seed before first iteration of algorithm if provided by user
    if(!is.null(rngseed) & is.numeric(rngseed)){
        set.seed(rngseed)
    }
    
    #####################
    ## Preprocess data ##
    #####################
    preprocout <- redist.preproc(adjobj = adjobj, popvec = popvec,
                                 initcds = initcds, ndists = ndists,
                                 popcons = popcons,
                                 grouppopvec = grouppopvec, ssdmat = ssdmat,
                                 beta = beta, temper = temper,
                                 constraint = constraint,
                                 betaseq = betaseq, betaseqlength = betaseqlength,
                                 betaweights = betaweights,
                                 adjswaps = adjswaps, maxiterrsg = maxiterrsg,
                                 contiguitymap = contiguitymap)

    ## Set betas - if tempering, modified later
    betapop <- preprocout$params$betapop
    betacompact <- preprocout$params$betacompact
    betaseg <- preprocout$params$betaseg
    betasimilar <- preprocout$params$betasimilar

    ## Get starting loop value
    loopstart <- loopscompleted + 1
    
    #######################
    ## Run the algorithm ##
    #######################
    for(i in loopstart:nloop){

        ## Get congressional districts, tempered beta values
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

            ## Reload the data if re-startomg
            if(loopstart > 1){
                
                ## Load the data
                load(paste(savename, "_loop", i - 1, ".RData", sep = ""))
                
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
            }
            
        }        

        ## Run algorithm
        algout <- swMH(aList = preprocout$data$adjlist,
                       cdvec = cds,
                       cdorigvec = preprocout$data$initcds,
                       popvec = preprocout$data$popvec,
                       grouppopvec = preprocout$data$grouppopvec,
                       nsims = nsims,
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
                       adjswap = preprocout$params$adjswaps,
                       exact_mh = exact_mh)

        class(algout) <- "redist"

        ## Save random number state if setting the seed
        if(!is.null(rngseed)){
            algout$randseed <- .Random.seed
        }
        
        ## Save output
        if(nloop > 1){
            save(algout, file = paste(savename, "_loop", i, ".RData", sep = ""))
        }
        
    }

    ####################
    ## Annealing flag ##
    ####################
    temperflag <- ifelse(preprocout$params$temperbetapop == 1 |
                             preprocout$params$temperbetacompact == 1 |
                                 preprocout$params$temperbetaseg == 1 |
                                     preprocout$params$temperbetasimilar == 1,
                         1, 0)

    ###############################
    ## Combine and save the data ##
    ###############################
    if(nloop > 1){
        redist.combine(savename = savename, nsims = nsims, nloop = nloop,
                       nthin = nthin, nunits = length(preprocout$data$adjlist),
                       temper = temperflag)
    }else if(!is.null(savename)){
        save(algout, file = paste(savename, ".RData", sep = ""))
    }

    ## Examine the data
    if(nloop == 1){
        return(algout)
    }

}

redist.diagplot <- function(sumstat,
                            plot = c("trace", "autocorr", "densplot",
                                "mean", "gelmanrubin"),
                            logit = FALSE, savename = NULL
                            ){

    ########################
    ## Inputs to function ##
    ########################
    ## rediststat - a vector of a summary stat for each redistricting plan
    ## plot - type of plot to create. Inputs are "trace," "autocorr", "densplot",
    ##        "mean", or "gelmanrubin"
    ## logit - logit transformation on dissimilarity index. Default = FALSE
    ## savename - Name to save under. Default = NULL

    ##############
    ## Warnings ##
    ##############
    if(missing(sumstat)){
        stop("Please provide a vector or list of summary statistics to the function")
    }
    if(!(class(sumstat) %in% c("numeric", "list", "mcmc", "mcmc.list"))){
        stop("Please provide either a numeric vector, list, or mcmc object")
    }
    if(!(plot %in% c("trace", "autocorr", "densplot",
                     "mean", "gelmanrubin"))){
        stop("Sorry. We don't currently support that MCMC diagnostic.")
    }
    if(plot == "gelmanrubin" & !(class(sumstat) %in% c("list", "mcmc.list"))){
        stop("If generating a Gelman-Rubin plot, please provide an object of class list or mcmc.list")
    }
    
    ########################
    ## Create mcmc object ##
    ########################
    if(class(sumstat) == "numeric"){
        segout <- mcmc(sumstat)
    }else if(class(sumstat) == "list"){
        for(i in 1:length(sumstat)){
            sumstat[[i]] <- mcmc(sumstat[[i]])
        }       
        segout <- mcmc.list(sumstat)
    }else if(class(sumstat) %in% c("mcmc", "mcmc.list")){
        segout <- sumstat
    }
    
    ## Logit transform
    if(logit){
        if(class(segout) == "mcmc"){
            segout <- log(segout / (1 - segout))
        }else if(class(segout) == "mcmc.list"){
            for(i in 1:length(segout)){
                segout[[i]] <- log(segout[[i]] / (1 - segout[[i]]))
            }
        }
    }

    ##################
    ## Create plots ##
    ##################
    if(plot == "trace"){
        if(!is.null(savename)){
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        traceplot(segout)
        if(!is.null(savename)){
            dev.off()
        }
    }
    if(plot == "autocorr"){
        if(!is.null(savename)){
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        autocorr.plot(segout, lag.max = 50)
        if(!is.null(savename)){
            dev.off()
        }
    }
    if(plot == "densplot"){
        if(!is.null(savename)){
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        densplot(segout)
        if(!is.null(savename)){
            dev.off()
        }
    }
    if(plot == "mean"){
        if(!is.null(savename)){
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        cumuplot(segout, probs = .5, type = "l", lty = 1)
        if(!is.null(savename)){
            dev.off()
        }
    }
    if(plot == "gelmanrubin" & class(segout) == "mcmc.list"){
        if(!is.null(savename)){
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        gelman.plot(segout, transform = FALSE)
        if(!is.null(savename)){
            dev.off()
        }
    }

}

redist.ipw <- function(algout,
                       resampleconstraint = c("pop", "compact",
                           "segregation", "similar"),
                       targetbeta,
                       targetpop = NULL,
                       temper = 0)
{

    #########################
    ## Inputs to function: ##
    #########################
    ## algout - output of redist.mcmc. Must be of class redist
    ## resampleconstraint - index we constraint on to form target distribution.
    ## targetbeta - target beta value
    ## targetpop - specified target population. Default is null
    ## temper - whether we are tempering on the resampleconstraint

    ## Warnings:
    if(missing(algout) | class(algout) != "redist"){
        stop("Please provide a proper redist object")
    }
    if(!(resampleconstraint %in% c("pop", "compact", "segregation", "similar"))){
        stop("We do not provide support for that constraint at this time")
    }
    if(missing(targetbeta)){
        stop("Please specify the target beta value")
    }

    ## Get indices drawn under target beta if tempering
    if(temper == 1){
        indbeta <- which(algout$beta_sequence == targetbeta)
    }else{
        indbeta <- 1:ncol(algout$partitions)
    }

    ## Get indices of draws that meet target population
    if(!is.null(targetpop)){
        indpop <- which(algout$distance_parity <= targetpop)
    }else{
        indpop <- 1:ncol(algout$partitions)
    }
    
    ## Get intersection of indices
    inds <- intersect(indpop, indbeta)

    ## Construct weights
    psi <- algout[[paste("constraint_", resampleconstraint, sep = "")]][inds]
    weights <- 1 / exp(targetbeta * psi)

    ## Resample indices
    inds <- sample(1:length(inds), length(inds), replace = TRUE, prob = weights)

    ## Subset the entire list
    algout_new <- vector(mode = "list", length = length(algout))
    for(i in 1:length(algout_new)){

        ## Subset the matrix first, then the vectors
        if(i == 1){
            algout_new[[i]] <- algout[[i]][,inds]
        }else{
            algout_new[[i]] <- algout[[i]][inds]
        }
        
    }
    names(algout_new) <- names(algout)

    ## Change class
    class(algout_new) <- "redist"

    return(algout_new)
    
}

redist.segcalc <- function(algout,
                           grouppop,
                           fullpop)
{

    ########################
    ## Inputs to function ##
    ########################
    ## algout: a matrix of congressional district assignments or a redist object
    ## grouppop: a vector of populations for some subgroup
    ## fullpop: a vector of populations for a geographic district

    ## Warnings
    if(missing(algout) | !(class(algout) %in% c("data.frame", "matrix", "redist"))){
        stop("Please provide either a redist object or a proper matrix of congessional districts")
    }
    if(missing(grouppop)){
        stop("Please provide a vector of sub-group populations to calculate
the segregation index")
    }
    if(missing(fullpop)){
        stop("Please provide a vector of populations for each geographic unit")
    }

    ## If redist object, get the partitions entry
    if(class(algout) == "redist"){
        algout <- algout$partitions
    }
    
    if(!((nrow(algout) == length(grouppop)) &
             (length(grouppop) == length(fullpop)) &
                 (length(fullpop) == nrow(algout)))){
        stop("Please make sure there is a population entry for each geographic unit")
    }

    ## Calculate dissimilarity index
    seg.out <- segregationcalc(algout,
                               grouppop,
                               fullpop)

    ## Return
    return(seg.out)
    
}
    
