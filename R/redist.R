###########################################
## Author: Ben Fifield
## Institution: Princeton University
## Date Created: 2015/02/04
## Date Modified: 2015/03/09
## Purpose: R wrapper to run swMH() code (non-mpi)
###########################################

combine.par.anneal <- function(a, b){

  ## Names of object
  name_out <- names(a)

  ## Create output object
  output_obj <- vector(mode = "list", length = length(a))

  ## Combine partitions
  for(i in 1:length(a)){
    if(i == i){
      output_obj[[i]] <- cbind(a[[i]], b[[i]])
    }else{
      output_obj[[i]] <- c(a[[i]], b[[i]])
    }
  }

  names(output_obj) <- name_out
  return(output_obj)

}




redist.preproc <- function(adjobj, popvec, initcds = NULL, ndists = NULL,
                           popcons = NULL,
                           countymembership = NULL,
                           grouppopvec = NULL,
                           areasvec = NULL,
                           borderlength_mat = NULL, ssdmat = NULL,
                           compactness_metric = NULL,
                           temper = NULL, constraint = NULL,
                           constraintweights = constraintweights,
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
  if(!is.null(constraintweights)){
    if((any(constraintweights == 0) & !is.null(constraint))){
      stop("If applying constraints or using simulated tempering, please set non-zero constraint by specifying the 'constraintweight' argument, and specify the names of the constraints in 'constraint'.")
    }
    if(any(!(constraint %in% c("compact", "segregation", "population", "similarity", "countysplit")))){
      stop("Please specify any combination of `compact`, `segregation`, `population`, `countysplit`, or `similarity` for constraint")
    }
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
      if(sum(is.na(initcds)) > 0){
        stop("You have NA's in your congressional districts. Please check the provided initcds vector for NA entries.")
      }

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
  if(is.null(initcds)){
    ## Set up target pop, strength of constraint (5%)
    if(is.null(popcons)){
      popcons_rsg <- .05
    }else{
      popcons_rsg <- popcons
    }

    ## Print start
    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    cat("\n", append = TRUE)
    cat(divider, append = TRUE)
    cat("Using redist.rsg() to generate starting values.\n\n", append= TRUE)

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
  if("segregation" %in% constraint & is.null(grouppopvec)){
    stop("If applying the segregation constraint, please provide a vector
             of subgroup populations")
  }
  if("countysplit" %in% constraint & is.null(countymembership)){
    stop("If applying the county split constraint, please provide a numeric vector indicating county membership.")
  }
  if("segregation" %in% constraint & !(is.null(grouppopvec))){
    if((length(grouppopvec) != length(adjlist)) |
       (sum(is.na(grouppopvec)) > 0)){
      stop("If applying the segregation constraint, each entry in adjacency
              list must have corresponding entry in vector of group populations")
    }
  }
  if("countysplit" %in% constraint & !is.null(countymembership)){
    if(length(countymembership) != length(adjlist) | sum(is.na(countymembership)) > 0){
      stop("You do not have a county membership assigned for every unit.")
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

  ## -------------------------------------
  ## Set county membership if not provided
  ## -------------------------------------
  if(is.null(countymembership)){
    countymembership <- c(0, 0, 0, 0)
  }else{
    if(is.factor(countymembership)){
      countymembership <- as.numeric(countymembership)
    }
    countymembership <- countymembership - min(countymembership)
  }

  ################################
  ## Set ssdmat if not provided ##
  ################################
  if(is.null(ssdmat) & "compact" %in% constraint & compactness_metric == "fryer-holden"){
    if(class(adjobj) == "SpatialPolygonsDataFrame"){
      centroids <- coordinates(adjobj)
      ssdmat <- calcPWDh(centroids)
    }else{
      stop("Provide squared distances matrix if constraining compactness using the Fryer-Holden metric.")
    }
  }else if(is.null(ssdmat)){
    ssdmat <- matrix(1, 2, 2)
  }

  ## ------------------------------------
  ## Set Polsby-Popper compactness inputs
  ## ------------------------------------
  if("compact" %in% constraint & compactness_metric == "polsby-popper"){
    if(is.null(areasvec) | is.null(borderlength_mat)){
      stop("If constraining on Polsby-Popper compactness, please provide both a vector of the areas of each geographic unit and a list with the border lengths of each pair of points.")
    }
    if(length(areasvec) != length(adjlist)){
      stop("The lengths of the areas vector and the adjacency list do not add up.")
    }
  }else{
    areasvec <- c(0, 0, 0, 0)
    borderlength_mat <- matrix(0, 2, 2)

  }

  ########################
  ## Set up constraints ##
  ########################
  beta <- ifelse(is.null(constraint) | temper, 0, 1)
  temperbeta <- ifelse(temper, "tempering", "none")

  if("population" %in% constraint){
    weightpop <- constraintweights[which(constraint == "population")]
  }else{
    weightpop <- 0
  }
  if("compact" %in% constraint){
    weightcompact <- constraintweights[which(constraint == "compact")]
  }else{
    weightcompact <- 0
  }
  if("segregation" %in% constraint){
    weightseg <- constraintweights[which(constraint == "segregation")]
  }else{
    weightseg <- 0
  }
  if("similarity" %in% constraint){
    weightsimilar <- constraintweights[which(constraint == "similarity")]
  }else{
    weightsimilar <- 0
  }
  if("countysplit" %in% constraint){
    weightcountysplit <- constraintweights[which(constraint == "countysplit")]
  }else{
    weightcountysplit <- 0
  }

  ###################################
  ## Check if betaspacing provided ##
  ###################################
  if(temperbeta == "tempering"){
    if(betaseq[1] == "powerlaw"){

      ## Generate power law sequence
      betaseq <- rep(NA, betaseqlength)
      for(i in 1:length(betaseq)){
        betaseq[i] <- (0.1^((i-1) / (length(betaseq) - 1)) - .1) / .9
      }

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
  preprocout <- list(
    data = list(
      adjlist = adjlist,
      popvec = popvec,
      initcds = initcds,
      grouppopvec = grouppopvec,
      areasvec = areasvec,
      borderlength_mat = borderlength_mat,
      ssdmat = ssdmat,
      countymembership = countymembership
    ),
    params = list(
      pctdistparity = popcons,
      dists = dists,
      beta = beta,
      temperbeta = temperbeta,
      betaseq = betaseq,
      betaweights = betaweights,
      adjswaps = adjswaps,
      weightpop = weightpop,
      weightcompact = weightcompact,
      weightseg = weightseg,
      weightsimilar = weightsimilar,
      weightcountysplit = weightcountysplit
    )
  )

  class(preprocout) <- "redist"

  return(preprocout)
}


#' Combine successive runs of \code{redist.mcmc}
#'
#' \code{redist.combine} is used to combine successive runs of \code{redist.mcmc}
#' into a single data object
#'
#' @usage redist.combine(savename, nloop, nthin, temper)
#'
#' @param savename The name (without the loop or \code{.RData} suffix)
#' of the saved simulations.
#' @param nloop The number of loops being combined.
#' @param nthin How much to thin the simulations being combined.
#' @param temper Wheterh simulated tempering was used (1) or not (0)
#' in the simulations. Default is 0.
#'
#' @details This function allows users to combine multiple successive runs of
#' \code{redist.mcmc} into a single \code{redist} object for analysis.
#'
#' @return \code{redist.combine} returns an object of class "redist". The object
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
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' @examples
#' \dontrun{
#' data(algdat.pfull)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins,Imai and
#' Tarr (2015)
#'
#' ## Get an initial partition
#' set.seed(1)
#' initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]
#'
#' ## Run the algorithm
#' alg_253 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
#' popvec = algdat.pfull$precinct.data$pop,
#' initcds = initcds,
#' nsims = 10000, nloops = 2, savename = "test")
#' out <- redist.combine(savename = "test", nloop = 2,
#' nthin = 10)
#' }
#' @export
redist.combine <- function(savename, nloop, nthin, temper = 0){

  ##############################
  ## Set up container objects ##
  ##############################
  load(paste(savename, "_loop1.RData", sep = ""))
  names_obj <- names(algout)

  ## Create containers
  nr <- nrow(algout$partitions)
  nc <- ncol(algout$partitions)
  partitions <- matrix(NA, nrow = nr,
                       ncol = (nc * nloop / nthin))

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
    load(paste(savename, "_loop", i, ".RData", sep = ""))

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


#' redist.combine.anneal
#'
#' Combine files generated by redist.mcmc(algorithm='anneal')
#'
#' @usage redist.combine.anneal(file_name)
#'
#' @param file_name The file name to search for in current working directory.
#'
#' @export
redist.combine.anneal <- function(file_name){

  ## List files
  fn <- list.files()[grep(file_name, list.files())]
  if(length(fn) == 0){
    stop("Can't find any files in current working directory with that name.")
  }
  load(fn[1])
  names_obj <- names(algout)

  # Create containers
  nr <- length(algout$partitions)
  nc <- length(fn)
  partitions <- matrix(NA, nrow = nr, ncol = nc)

  veclist <- vector(mode = "list", length = length(algout)-1)
  for(i in 1:length(veclist)){
    veclist[[i]] <- rep(NA, nc)
  }

  ## ------------
  ## Combine data
  ## ------------
  for(i in 1:length(fn)){
    ## Load data
    load(fn[i])

    ## Store objects together
    for(j in 1:length(algout)){
      if(j == 1){
        partitions[1:nr, i] <- algout$partitions
      }else{
        veclist[[j-1]][i] <- algout[[j]]
      }
    }
  }

  ## ---------------------------
  ## Store data in algout object
  ## ---------------------------
  algout <- vector(mode = "list", length = length(algout))
  for(i in 1:length(algout)){
    if(i == 1){
      algout[[i]] <- partitions
    }else{
      algout[[i]] <- veclist[[i-1]]
    }
  }
  names(algout) <- names_obj

  ## -------------
  ## Output object
  ## -------------
  class(algout) <- "redist"
  return(algout)

}

# Helpers for MPI
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

} # ends function redist.combine.mpi

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
} # ends function ecutsMPI




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
} # ends function ecutsAppend


#' MCMC Redistricting Simulator
#'
#' \code{redist.mcmc} is used to simulate Congressional redistricting
#' plans using Markov Chain Monte Carlo methods.
#'
#' @usage redist.mcmc(algorithm, adjobj, popvec, nsims, ndists = NA,
#' initcds = NULL, loopscompleted = 0, nloop = 1, nthin = 1,
#' eprob = 0.05, lambda = 0, popcons = NA, grouppopvec = NA,
#' areasvec, countymembership, borderlength_mat, ssdmat = NA,
#' temper, constraint = "population", constraintweights, compactness_metric,
#' betaseq, betaseqlength = 10, betaweights, adjswaps = TRUE,
#' rngseed = NA, beta = -10, maxiterrsg = 5000, freq = 100,
#' adapt_lambda, adapt_eprob, contiguitymap = "rooks", exact_mh, savename = NULL,
#' verbose = TRUE, num_hot_steps, num_annealing_steps, num_cold_steps, ncores)
#'
#' @param algorithm The algorithm to run. Default is \code{mcmc}, the other
#' implemented are \code{anneal} and \code{mpi}.
#' @param adjobj An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param popvec A vector containing the populations of each geographic
#' unit
#' @param nsims The number of simulations run before a save point.
#' @param ndists The number of congressional districts. The default is
#' \code{NULL}.
#' @param initcds A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided,
#' random and contiguous congressional district assignments will be generated
#' using \code{redist.rsg}.
#' @param loopscompleted Number of save points reached by the
#' algorithm. The default is \code{0}.
#' @param nloop The total number of save points for the algorithm. The
#' default is \code{1}. Note that the total number of simulations run
#' will be \code{nsims} * \code{nloop}.
#' @param nthin The amount by which to thin the Markov Chain. The
#' default is \code{1}.
#' @param eprob The probability of keeping an edge connected. The
#' default is \code{0.05}.
#' @param lambda The parameter determining the number of swaps to attempt
#' each iteration of the algorithm. The number of swaps each iteration is
#' equal to Pois(\code{lambda}) + 1. The default is \code{0}.
#' @param popcons The strength of the hard population
#' constraint. \code{popcons} = 0.05 means that any proposed swap that
#' brings a district more than 5\% away from population parity will be
#' rejected. The default is \code{NULL}.
#' @param grouppopvec A vector of populations for some sub-group of
#' interest. The default is \code{NULL}.
#' @param areasvec A vector of precinct areas for discrete Polsby-Popper.
#' The default is \code{NULL}.
#' @param countymembership A vector of county membership assignments. The default is \code{NULL}.
#' @param borderlength_mat A matrix of border length distances, where
#' the first two columns are the indices of precincts sharing a border and
#' the third column is its distance. Default is \code{NULL}.
#' @param ssdmat A matrix of squared distances between geographic
#' units. The default is \code{NULL}.
#' @param temper Whether to use simulated tempering algorithm. Default is FALSE.
#' @param constraint Which constraint to apply. Accepts any combination of \code{compact},
#' \code{segregation}, \code{population}, \code{similarity}, or \code{none}
#' (no constraint applied). The default is NULL.
#' @param constraintweights The weights to apply to each constraint. Should be a vector
#' the same length as constraint. Default is NULL.
#' @param compactness_metric The compactness metric to use when constraining on
#' compactness. Default is \code{fryer-holden}, the other implemented option
#' is \code{polsby-popper}.
#' @param betaseq Sequence of beta values for tempering. The default is
#' \code{powerlaw} (see Fifield et. al (2015) for details).
#' @param betaseqlength Length of beta sequence desired for
#' tempering. The default is \code{10}.
#' @param betaweights Sequence of weights for different values of
#' beta. Allows the user to upweight certain values of beta over
#' others. The default is \code{NULL} (equal weighting).
#' @param adjswaps Flag to restrict swaps of beta so that only
#' values adjacent to current constraint are proposed. The default is
#' \code{TRUE}.
#' @param rngseed Allows the user to set the seed for the
#' simulations. Default is \code{NULL}.
#' @param beta The strength of the target strength in the MH ratio. The default
#' is 0. This is only used for algorithm mpi.
#' @param maxiterrsg Maximum number of iterations for random seed-and-grow
#' algorithm to generate starting values. Default is 5000.
#' @param freq Frequency of between-chain swaps. Default to once every 100
#' iterations. This is only used for algorithm mpi.
#' @param adapt_lambda Whether to adaptively tune the lambda parameter so that the Metropolis-Hastings
#' acceptance probability falls between 20\% and 40\%. Default is FALSE.
#' @param adapt_eprob Whether to adaptively tune the edgecut probability parameter so that the
#' Metropolis-Hastings acceptance probability falls between 20\% and 40\%. Default is
#' FALSE.
#' @param contiguitymap Use queens or rooks distance criteria for generating an
#' adjacency list from a "SpatialPolygonsDataFrame" data type.
#' Default is "rooks".
#' @param exact_mh Whether to use the approximate (0) or exact (1)
#' Metropolis-Hastings ratio calculation for accept-reject rule. Default is FALSE.
#' @param savename Filename to save simulations. Default is \code{NULL}.
#' @param verbose Whether to print initialization statement.
#' Default is \code{TRUE}.
#' @param num_hot_steps The number of steps to run the simulator at beta = 0.
#' Default is 40000. This is only used for algorithm anneal.
#' @param num_annealing_steps The number of steps to run the simulator with
#' linearly changing beta schedule. Default is 60000. This is only used for
#' algorithm anneal.
#' @param num_cold_steps The number of steps to run the simulator at beta = 1.
#' Default is 20000. This is only used for algorithm anneal.
#' @param ncores The number of cores available to parallelize over. Default is 1.
#' This is only used for algorithm anneal.
#'
#' @details This function allows users to simulate redistricting plans
#' using Markov Chain Monte Carlo methods. Several constraints
#' correspoding to substantive requirements in the redistricting process
#' are implemented, including population parity and geographic
#' compactness. In addition, the function includes multiple-swap and
#' simulated tempering functionality to improve the mixing of the Markov
#' Chain.
#'
#' @return \code{redist.mcmc} returns an object of class "redist". The object
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
#' alg_253 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
#' popvec = algdat.pfull$precinct.data$pop,
#' initcds = initcds,
#' nsims = 10000)
#' }
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
#' redist.mcmc(algorithm = "mpi", adjobj = algdat.pfull$adjlist,
#' popvec = algdat.pfull$precinct.data$pop,
#' initcds = initcds,
#' nsims = 10000, savename = "test")
#' }
#' @export
redist.mcmc <- function(adjobj, popvec, nsims, algorithm = c("mcmc", "anneal", "mpi"),
                        ndists = NULL, initcds = NULL,
                        loopscompleted = 0, nloop = 1, nthin = 1, eprob = 0.05,
                        lambda = 0, popcons = NULL, grouppopvec = NULL,
                        areasvec = NULL, countymembership = NULL,
                        borderlength_mat = NULL, ssdmat = NULL, temper = FALSE,
                        constraint = NULL, constraintweights = NULL,
                        compactness_metric = "fryer-holden",
                        betaseq = "powerlaw", betaseqlength = 10,
                        betaweights = NULL, beta = -10, freq = 100,
                        adjswaps = TRUE, rngseed = NULL, maxiterrsg = 5000,
                        adapt_lambda = FALSE, adapt_eprob = FALSE,
                        contiguitymap = "rooks", exact_mh = FALSE, savename = NULL,
                        verbose = TRUE,
                        num_hot_steps = 40000, num_annealing_steps = 60000,
                        num_cold_steps = 20000, ncores = 1
){

  if(verbose){
    ## Initialize ##
    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    cat("\n", append = TRUE)
    cat(divider, append = TRUE)
    cat("redist.mcmc(): Automated Redistricting Simulation Using
         Markov Chain Monte Carlo\n\n", append = TRUE)
  }

  if(algorithm == "mpi"){
  ## Check if Rmpi library is installed
  if (!requireNamespace("Rmpi", quietly = TRUE)) {
      stop("You must install package 'Rmpi' to use this function. Please install it if you wish to continue."
          ,call. = FALSE)
  }

  ## ## Load Rmpi library
  ## if (!is.loaded("mpi_initialize")) {
  ##     library("Rmpi")
  ## }
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
  algorithm <- match.arg(algorithm)
  if(missing(nsims) & (algorithm %in% c('mcmc','mpi'))){
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
  if(!is.null(constraint) & is.null(constraintweights)){
    stop("Please provide a weight value in 'constraintweights' for each constraint specified in 'constraint'.")
  }
  if(!(compactness_metric %in% c("fryer-holden", "polsby-popper"))){
    stop("We only support either 'fryer-holden' or 'polsby-popper' as compactness metrics.")
  }

  ## Set seed before first iteration of algorithm if provided by user
  if(!is.null(rngseed) & is.numeric(rngseed)){
    set.seed(rngseed)
  }

  if(adapt_lambda){
    adapt_lambda <- 1
  }else{
    adapt_lambda <- 0
  }
  if(adapt_eprob){
    adapt_eprob <- 1
  }else{
    adapt_eprob <- 0
  }
  if(exact_mh){
    exact_mh <- 1
  }else{
    exact_mh <- 0
  }

  #####################
  ## Preprocess data ##
  #####################
  if(algorithm %in% c('mcmc','anneal')){
  cat("Preprocessing data.\n\n")
  preprocout <- redist.preproc(adjobj = adjobj, popvec = popvec,
                               initcds = initcds, ndists = ndists,
                               popcons = popcons,
                               countymembership = countymembership,
                               grouppopvec = grouppopvec,
                               areasvec = areasvec,
                               borderlength_mat = borderlength_mat,
                               ssdmat = ssdmat,
                               compactness_metric = compactness_metric,
                               temper = temper, # Default anneal is False
                               constraint = constraint, constraintweights = constraintweights,
                               betaseq = betaseq, #Default anneal is 'powerlaw'
                               betaseqlength = betaseqlength, #Default anneal is 10
                               betaweights = betaweights, #Default anneal is NULL
                               adjswaps = adjswaps, #Default anneal is TRUE
                               maxiterrsg = maxiterrsg,
                               contiguitymap = contiguitymap)

  ## Set betas - if tempering, modified later
  weightpop <- preprocout$params$weightpop
  weightcompact <- preprocout$params$weightcompact
  weightseg <- preprocout$params$weightseg
  weightsimilar <- preprocout$params$weightsimilar
  weightcountysplit <- preprocout$params$weightcountysplit
  } else {
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
  }

  if(algorithm == 'mpi'){
      ## Create parameters list to distribute across nodes
      params <- expand.grid(nsims = nsims,nloop = nloop,eprob = eprob,
                            ndists = ndists,lambda = lambda,popcons = popcons,
                            beta = beta,target.beta = target.beta,constraint = constraint,
                            betaseqlength = betaseqlength,adjswaps = adjswaps,
                            nthin = nthin,freq = freq,maxiterrsg = maxiterrsg,
                            contiguitymap = contiguitymap,verbose = verbose,
                            loopscompleted = loopscompleted,rngseed = rngseed,
                            savename = savename)
  }

  if(algorithm == 'mcmc'){
  ## Get starting loop value
  loopstart <- loopscompleted + 1

  #######################
  ## Run the algorithm ##
  #######################
  for(i in loopstart:nloop){

    ## Get congressional districts, tempered beta values
    if(i > loopstart){

      cds <- algout$partitions[,nsims]

      if(temper){
        beta <- algout$beta_sequence[nsims]
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

        if(temper){
          beta <- algout$beta_sequence[nsims]
        }

        if(!is.null(rngseed) & is.numeric(rngseed)){
          .Random.seed <- algout$randseed
        }

        rm(list = "algout")

      }else{
        cds <- preprocout$data$initcds
      }
    } #ends else
  }
  } #ends if mcmc

  if(algorithm %in% c('mcmc','anneal')){
    ## Run algorithm
    cat("Starting swMH().\n")
    algout <- swMH(aList = preprocout$data$adjlist,
                   cdvec = preprocout$data$initcds,
                   cdorigvec = preprocout$data$initcds,
                   popvec = preprocout$data$popvec,
                   grouppopvec = preprocout$data$grouppopvec,
                   areas_vec = preprocout$data$areasvec,
                   county_membership = preprocout$data$countymembership,
                   borderlength_mat = preprocout$data$borderlength_mat,
                   nsims = 100,
                   eprob = eprob,
                   pct_dist_parity = preprocout$params$pctdistparity,
                   beta_sequence = preprocout$params$betaseq,
                   beta_weights = preprocout$params$betaweights,
                   ssdmat = preprocout$data$ssdmat,
                   lambda = lambda,
                   beta = 0,
                   weight_population = weightpop,
                   weight_compact = weightcompact,
                   weight_segregation = weightseg,
                   weight_similar = weightsimilar,
                   weight_countysplit = weightcountysplit,
                   adapt_beta = algorithm,
                   adjswap = preprocout$params$adjswaps,
                   exact_mh = exact_mh,
                   adapt_lambda = adapt_lambda,
                   adapt_eprob = adapt_eprob,
                   compactness_measure = compactness_metric,
                   num_hot_steps = num_hot_steps,
                   num_annealing_steps = num_annealing_steps,
                   num_cold_steps = num_cold_steps)

    class(algout) <- "redist"

    ## Save random number state if setting the seed
    if(!is.null(rngseed)&algorithm=='mcmc'){
      algout$randseed <- .Random.seed
    }

    ## Save output
    if(nloop > 1){
      save(algout, file = paste(savename, "_loop", i, ".RData", sep = ""))
    }

  ####################
  ## Annealing flag ##
  ####################
  temperflag <- ifelse(preprocout$params$temperbeta == "tempering", 1, 0)

  ###############################
  ## Combine and save the data ##
  ###############################
  if(nloop > 1){
    redist.combine(savename = savename, nloop = nloop,
                   nthin = nthin,
                   temper = temperflag)
  }else if(!is.null(savename)){
    save(algout, file = paste(savename, ".RData", sep = ""))
  }

  ## Examine the data
  if(nloop == 1 | algorithm == 'anneal'){
    return(algout)
  }
    } else{
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
}

#' Inverse probability reweighting for MCMC Redistricting
#'
#' \code{redist.ipw} properly weights and resamples simulated redistricting plans
#' so that the set of simulated plans resemble a random sample from the
#' underlying distribution. \code{redist.ipw} is used to correct the sample when
#' population parity, geographic compactness, or other constraints are
#' implemented.
#'
#' @usage redist.ipw(algout, targetpop = NULL)
#'
#' @param algout An object of class "redist".
#' @param targetpop The desired level of population parity. \code{targetpop} =
#' 0.01 means that the desired distance from population parity is 1\%. The
#' default is \code{NULL}.
#'
#' @details This function allows users to resample redistricting plans using
#' inverse probability weighting techniques described in Rubin (1987). This
#' techniques reweights and resamples redistricting plans so that the resulting
#' sample is representative of a random sample from the uniform distribution.
#'
#' @return \code{redist.ipw} returns an object of class "redist". The object
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
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain
#' Monte Carlo." Working Paper.
#' Available at \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' Rubin, Donald. (1987) "Comment: A Noniterative Sampling/Importance Resampling
#' Alternative to the Data Augmentation Algorithm for Creating a Few Imputations
#' when Fractions of Missing Information are Modest: the SIR Algorithm."
#' Journal of the American Statistical Association.
#'
#' @examples \dontrun{
#' data(algdat.p20)
#'
#' ## Code to run the simulations in Figure 4 of Fifield, Higgins,
#' ## Imai and Tarr (2015)
#'
#' ## Get an initial partition
#' set.seed(1)
#' initcds <- algdat.p20$cdmat[,sample(1:ncol(algdat.p20$cdmat), 1)]
#'
#' ## Vector of beta weights
#' betaweights <- rep(NA, 10); for(i in 1:10){betaweights[i] <- 4^i}
#'
#' ## Run simulations - tempering population constraint
#' alg_253_20_st <- redist.mcmc(adjobj = algdat.p20$adjlist,
#' popvec = algdat.p20$precinct.data$pop,
#' initcds = initcds, nsims = 10000, betapop = -5.4,
#' betaweights = betaweights, temperbetapop = 1)
#'
#' ## Resample using inverse probability weighting.
#' ## Target distance from parity is 20%
#' alg_253_20_st <- redist.ipw(alg_253_20_st, targetpop = .2)
#'
#' }
#' @export
redist.ipw <- function(algout, targetpop = NULL){

  ## Warnings:
  if(!inherits(algout, "redist")){
    stop("Please provide a proper redist object")
  }

  ## Get indices drawn under target beta if tempering
  indbeta <- which(algout$beta_sequence == 1)

  ## Get indices of draws that meet target population
  if(!is.null(targetpop)){
    indpop <- which(algout$distance_parity <= targetpop)
  }else{
    indpop <- 1:ncol(algout$partitions)
  }

  ## Get intersection of indices
  inds <- intersect(indpop, indbeta)

  ## Construct weights
  psi <- algout[["energy_psi"]][inds]
  weights <- 1 / exp(-1 * psi)

  ## Resample indices
  inds <- sample(inds, length(inds), replace = TRUE, prob = weights)

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
