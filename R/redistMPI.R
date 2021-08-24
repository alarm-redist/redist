##############################################
# Authors: Ben Fifield & Alex Tarr           #
# Created: 2015/06/01                        #
# Last Revision: N/A                         #
# Institution: Princeton University          #
# Purpose: R wrapper to run swMH() code w/   #
# parallel tempering (mpi)                   #
##############################################

ecutsMPI <- function(procID = procID, params = params, adj = adj, total_pop = total_pop, init_plan = init_plan, swaps = swaps){
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
    if(is.na(group_pop)){
        group_pop <- NULL
    }
    if(is.na(counties)){
        counties <- NULL
    }
    if(is.na(areasvec)){
        areasvec <- NULL
    }
    if(is.na(ssdmat)){
        ssdmat <- matrix(1, 2, 2)
    }
    if(is.na(borderlength_mat)){
        borderlength_mat <- NULL
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
        constraint <- strsplit(params$constraint, ",")[[1]]
    }
    if(is.na(params$constraintweights)){
        constraintweights <- NULL
    }else{
        constraintweights <- as.numeric(strsplit(params$constraintweights, ",")[[1]])
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
    if(is.na(params$pop_tol)){
        pop_tol <- NULL
    }else{
        pop_tol <- params$pop_tol
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
    if(sum(is.na(init_plan)) == length(init_plan)){
        init_plan <- NULL
    }
    if(is.na(params$compactness_metric)){
        compactness_metric <- NULL
    }else{
        compactness_metric <- params$compactness_metric
    }

    nthin <- params$nthin

    ## Run redist preprocessing function
    preprocout <- redist.preproc(adj = adj, total_pop = total_pop,
                                 init_plan = init_plan, ndists = ndists,
                                 pop_tol = pop_tol,
                                 group_pop = group_pop,
                                 areasvec = areasvec, #ctk-cran-note
                                 borderlength_mat = borderlength_mat,#ctk-cran-note
                                 counties = counties,#ctk-cran-note
                                 ssdmat = ssdmat,
                                 compactness_metric = compactness_metric,#ctk-cran-note
                                 temper = FALSE,
                                 constraint = constraint,
                                 constraintweights = constraintweights,#ctk-cran-note
                                 betaseq = NULL, betaweights = NULL,
                                 adjswaps = adjswaps,
                                 maxiterrsg = maxiterrsg,
                                 contiguitymap = contiguitymap)

    ## Set betas - if tempering, modified later
    beta <- params$beta

    weightpop <- preprocout$params$weightpop
    weightcompact <- preprocout$params$weightcompact
    weightvra <- preprocout$params$weightvra
    weightsimilar <- preprocout$params$weightsimilar
    weightcountysplit <- preprocout$params$weightcountysplit
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

            if(temper == "parallel"){
                beta <- algout$beta_sequence[nsims]
            }

            if(!is.null(rngseed) & is.numeric(rngseed)){
                set.seed(algout$randseed)
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

                if(temper == "parallel"){
                    beta <- algout$beta_sequence[nsims]
                }

                if(!is.null(rngseed) & is.numeric(rngseed)){
                    set.seed(algout$randseed)
                }

                rm(list = "algout")
                algout <- list()

            }else{
                cds <- preprocout$data$init_plan
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

            cat("Swap ", j, " out of ", length(nsimsAdj), " swaps.\n",
                append = TRUE)

            ## Run algorithm
            temp <- swMH(aList = preprocout$data$adjlist,
                         cdvec = cds,
                         cdorigvec = preprocout$data$init_plan,
                         popvec = preprocout$data$total_pop,
                         grouppopvec = preprocout$data$group_pop,
                         areas_vec = preprocout$data$areasvec,
                         county_membership = preprocout$data$counties,
                         borderlength_mat = preprocout$data$borderlength_mat,
                         nsims = nsimsAdj[j],
                         eprob = eprob,
                         pct_dist_parity = preprocout$params$pctdistparity,
                         beta_sequence = preprocout$params$betaseq,
                         beta_weights = preprocout$params$betaweights,
                         ssdmat = preprocout$data$ssdmat,
                         lambda = lambda,
                         beta = beta,
                         weight_population = weightpop,
                         weight_compact = weightcompact,
                         weight_vra = weightvra,
                         weight_similar = weightsimilar,
                         weight_countysplit = weightcountysplit,
                         adapt_beta = "none",
                         adjswap = preprocout$params$adjswaps,
                         exact_mh = 0,
                         adapt_eprob = 0,
                         adapt_lambda = 0,
                         compactness_measure = compactness_metric)

            ## Combine data
            algout <- ecutsAppend(algout,temp)

            ## Get temperature
            beta <- temp$beta_sequence[nsimsAdj[j]]

            ## Get likelihood
            like <- temp$energy_psi[nsimsAdj[j]]
            ## Use MPI to exchange swap information
            ##
            ## Tag guide: 1 -> likelihood sent
            ##            2 -> temperature sent
            ##            3 -> acceptance probability sent
            ##

            if(adjswaps){

                ## Determine which nodes are swapping
                tempvra <- swaps[(i-1)*nsims + j*freq]
                ## Get node indices
                temps <- tempadj[tempvra:(tempvra+1)]
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
                    cat("Components of likelihood:\n", append = TRUE)
                    cat("like =", like, "\n", append = TRUE)
                    cat("beta =", beta, "\n", append = TRUE)
                    cat("betaPart =", betaPart, "\n", append = TRUE)
                    cat("likePart =", likePart, "\n", append = TRUE)

                    a_like <- -1 * betaPart * like
                    b_like <- -1 * beta * likePart
                    x_like <- -1 * beta * like
                    y_like <- -1 * betaPart * likePart
                    prob <- exp(a_like + b_like - x_like - y_like)

                    if(prob > accept){
                        ## Exchange temperature values
                        beta <- betaPart

                        ## Adjust temperature adjacency list
                        tempadj[tempvra:(tempvra+1)] <- tempadj[(tempvra+1):tempvra]
                        ## Send temperature adjacency list
                        if(procID == tempadj[tempvra+1]){
                            oProcs <- tempadj[!(tempadj %in% temps)]
                            for(k in 1:length(oProcs)){
                                Rmpi::mpi.send.Robj(tempadj,dest=oProcs[k],tag=4)
                            }
                        }
                    }else{
                        if(procID == tempadj[tempvra]){
                            oProcs <- tempadj[!(tempadj %in% temps)]
                            for(k in 1:length(oProcs)){
                                Rmpi::mpi.send.Robj(tempadj,dest=oProcs[k],tag=4)
                            }
                        }
                    }
                }else{
                    tempadj <- Rmpi::mpi.recv.Robj(tempadj[tempvra],tag=4)
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

                    ## Compute acceptance probability
                    a_like <- -1 * betaPart * like
                    b_like <- -1 * beta * likePart
                    x_like <- -1 * beta * like
                    y_like <- -1 * betaPart * likePart
                    prob <- exp(a_like + b_like - x_like - y_like)
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
            algout$randseed <- .Random.seed[3]
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
        redist.combine.mpi(savename = savename, nloop = nloop,
                           nthin = nthin, tempadj = tempadj)
      }
    }else if(!is.null(savename)){
        if(procID == tempadj[1]){
            save(algout, file = paste(savename, ".RData", sep = ""))
        }
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
#' @usage redist.combine.mpi(savename, nloop, nthin, tempadj)
#'
#' @param savename The name (without the loop or \code{.RData} suffix)
#' of the saved simulations.
#' @param nloop The number of loops being combined.
#' @param nthin How much to thin the simulations being combined.
#' @param tempadj The temperature adjacency object saved by
#' \code{redist.mcmc.mpi}.
#'
#' @details This function allows users to combine multiple successive runs of
#' \code{redist.mcmc.mpi} into a single \code{redist} object for analysis.
#'
#' @return \code{redist.combine.mpi} returns an object of class "redist".
#' The object \code{redist} is a list that contains the following components (the
#' inclusion of some components is dependent on whether tempering
#' techniques are used):
#' \item{plans}{Matrix of congressional district assignments generated by the
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
#' \item{constraint_vra}{A vector containing the value of the
#' vra constraint for each accepted redistricting plan.}
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
#' # Cannot run on machines without Rmpi
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins, Imai and
#' ## Tarr (2015)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## Run the algorithm
#' redist.mcmc.mpi(adj = fl25_adj, total_pop = fl25$pop,
#'                 init_plan = init_plan, nsims = 10000, nloops = 2, savename = "test")
#' out <- redist.combine.mpi(savename = "test", nloop = 2,
#'                           nthin = 10, tempadj = tempAdjMat)
#' }
#' @concept post
#' @export
redist.combine.mpi <- function(savename, nloop, nthin, tempadj){

    ##############################
    ## Set up container objects ##
    ##############################
    load(paste(savename, "_proc", tempadj[1], "_loop1.RData", sep = ""))
    names_obj <- names(algout)

    ## Create containers
    nr <- nrow(algout$partitions)
    nc <- ncol(algout$partitions)
    partitions <- matrix(NA, nrow = nr, ncol = (nc * nloop / nthin))

    veclist <- vector(mode = "list", length = length(algout)-1)
    for(i in 1:length(veclist)){
        veclist[[i]] <- rep(NA, (nc * nloop / nthin))
    }

    ## Indices for thinning
    indthin <- which((1:nc) %% nthin == 0)

    ####################################
    ## Combine data in multiple loops ##
    ####################################

    for(i in 1:nloop){

        ## Load data
        load(paste(savename, "_proc", tempadj[1], "_loop", i, ".RData", sep = ""))

        ind <- ((i - 1) * (nc / nthin) + 1):(i * (nc / nthin))

        ## Store objects together
        for(j in 1:length(algout)){
            if(j == 1){
                partitions[1:nr, ind] <- algout$partitions[,indthin]
            }else{
                veclist[[j-1]][ind] <- algout[[j]][indthin]
            }
        }

    }

    #################################
    ## Store data in algout object ##
    #################################
    algout <- vector(mode = "list", length = length(algout))
    for(i in 1:length(algout)){
        if(i == 1){
            algout[[i]] <- partitions
        }else{
            algout[[i]] <- veclist[[i-1]]
        }
    }
    names(algout) <- names_obj


    #########################
    ## Set class of object ##
    #########################
    class(algout) <- "redist"

    #################
    ## Save object ##
    #################
    save(algout, file = paste(savename, ".RData", sep = ""))

}

ecutsAppend <- function(algout,ndata){
    if(length(algout) == 0){
        algout <- ndata
    }else{
        names_obj <- names(algout)
        for(i in 1:length(algout)){
            if(i == 1){
                algout[[i]] <- cbind(algout[[i]], ndata[[i]])
            }else{
                algout[[i]] <- c(algout[[i]], ndata[[i]])
            }
        }
        names(algout) <- names_obj
    }
    return(algout)
}

#' MCMC Redistricting Simulator using MPI
#'
#' \code{redist.mcmc.mpi} is used to simulate Congressional redistricting
#' plans using Markov Chain Monte Carlo methods.
#'
#'
#' @param adj An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param total_pop A vector containing the populations of each geographic
#' unit.
#' @param nsims The number of simulations run before a save point.
#' @param ndists The number of congressional districts. The default is
#' \code{NULL}.
#' @param init_plan A vector containing the congressional district labels
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
#' @param pop_tol The strength of the hard population
#' constraint. \code{pop_tol} = 0.05 means that any proposed swap that
#' brings a district more than 5\% away from population parity will be
#' rejected. The default is \code{NULL}.
#' @param group_pop A vector of populations for some sub-group of
#' interest. The default is \code{NULL}.
#' @param areasvec A vector of precinct areas for discrete Polsby-Popper.
#' The default is \code{NULL}.
#' @param counties A vector of county membership assignments. The default is \code{NULL}.
#' @param borderlength_mat A matrix of border length distances, where
#' the first two columns are the indices of precincts sharing a border and
#' the third column is its distance. Default is \code{NULL}.
#' @param ssdmat A matrix of squared distances between geographic
#' units. The default is \code{NULL}.
#' @param compactness_metric The compactness metric to use when constraining on
#' compactness. Default is \code{fryer-holden}, the other implemented option
#' is \code{polsby-popper}.
#' @param rngseed Allows the user to set the seed for the
#' simulations. Default is \code{NULL}.
#' @param constraint Which constraint to apply. Accepts any combination of \code{compact},
#' \code{vra}, \code{population}, \code{similarity}, or \code{none}
#' (no constraint applied). The default is NULL.
#' @param constraintweights The weights to apply to each constraint. Should be a vector
#' the same length as constraint. Default is NULL.
#' @param betaseq Sequence of beta values for tempering. The default is
#' \code{powerlaw} (see Fifield et. al (2015) for details).
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
#' @param verbose Whether to print initialization statement. Default is
#' \code{TRUE}.
#'
#'
#' @details This function allows users to simulate redistricting plans
#' using Markov Chain Monte Carlo methods. Several constraints
#' corresponding to substantive requirements in the redistricting process
#' are implemented, including population parity and geographic
#' compactness. In addition, the function includes multiple-swap and
#' parallel tempering functionality in MPI to improve the mixing of the Markov
#' Chain.
#'
#' @return \code{redist.mcmc.mpi} returns an object of class "redist". The object
#' \code{redist} is a list that contains the following components (the
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
#' \item{constraint_vra}{A vector containing the value of the
#' vra constraint for each accepted redistricting plan.}
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
#' @concept simulate
#' @examples
#' \dontrun{
#' # Cannot run on machines without Rmpi
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins, Imai and
#' ## Tarr (2015)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## Run the algorithm
#' redist.mcmc.mpi(adj = fl25_adj, total_pop = fl25$pop,
#'                 init_plan = init_plan, nsims = 10000, savename = "test")
#' }
#' @export
redist.mcmc.mpi <- function(adj, total_pop, nsims, ndists = NA,
                            init_plan = NULL,
                            loopscompleted = 0, nloop = 1, nthin = 1,
                            eprob = 0.05, lambda = 0,
                            pop_tol = NA, group_pop = NA, areasvec = NA,
                            counties = NA, borderlength_mat = NA,
                            ssdmat = NA, compactness_metric = "fryer-holden", rngseed = NA,
                            constraint = NA, constraintweights = NA,
                            betaseq = "powerlaw", betaseqlength = 10,
                            adjswaps = TRUE, freq = 100, savename = NA,
                            maxiterrsg = 5000, verbose = FALSE){
  contiguitymap <- 'rooks'

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
    if(missing(adj)){
        stop("Please supply adjacency matrix or list")
    }
    if(missing(total_pop)){
        stop("Please supply vector of geographic unit populations")
    }
    if(missing(nsims)){
        stop("Please supply number of simulations to run algorithm")
    }
    if(is.na(ndists) & is.null(init_plan)){
        stop("Please provide either the desired number of congressional districts
         or an initial set of congressional district assignments")
    }
    if(nloop > 1 & missing(savename)){
        stop("Please supply save directory if saving simulations at checkpoints")
    }

    ###################
    ## Preprocessing ##
    ###################

    ## Augment init_plan if necessary
    nrow.init <- ifelse(is.null(init_plan), 0, nrow(init_plan))
    ncol.init <- ifelse(is.null(init_plan), ndists, ncol(init_plan))
    if(nrow.init < betaseqlength){
        init_plan <- rbind(init_plan,matrix(NA,betaseqlength-nrow.init,ncol.init))
    }

    ## Generate temperature sequence (power law)
    temp <- rep(NA, betaseqlength)
    for(i in 1:betaseqlength){
        temp[i] <- 0.1^((i-1) / (betaseqlength - 1)) - .1
    }
    beta <- temp/0.9
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
                          ndists = ndists,lambda = lambda,pop_tol = pop_tol,
                          beta = beta,target.beta = target.beta,
                          constraint = paste(constraint, collapse = ","),
                          constraintweights = paste(constraintweights, collapse = ","),
                          compactness_metric = compactness_metric,
                          betaseqlength = betaseqlength,adjswaps = adjswaps,
                          nthin = nthin,freq = freq,maxiterrsg = maxiterrsg,
                          contiguitymap = contiguitymap,verbose = verbose,
                          loopscompleted = loopscompleted,rngseed = rngseed,
                          savename = savename, stringsAsFactors = FALSE)

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
    Rmpi::mpi.bcast.Robj2slave(adj)

    ## Population Vector
    Rmpi::mpi.bcast.Robj2slave(total_pop)

    ## Initial Plans
    init_plan <- split(init_plan, f=1:nrow(init_plan))
    Rmpi::mpi.scatter.Robj2slave(init_plan)

    ## Group population vector
    Rmpi::mpi.bcast.Robj2slave(group_pop)

    ## Areas vector
    Rmpi::mpi.bcast.Robj2slave(areasvec)

    ## County memberhsip vector
    Rmpi::mpi.bcast.Robj2slave(counties)

    ## Border distance matrix
    Rmpi::mpi.bcast.Robj2slave(borderlength_mat)

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
    Rmpi::mpi.bcast.cmd(ecutsMPI(procID, params, adj, total_pop, init_plan, swaps))

    ## Close slaves
    Rmpi::mpi.close.Rslaves()

    ## Terminate MPI processes and close R
    Rmpi::mpi.quit()

}
