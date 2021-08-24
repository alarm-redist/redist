redist.preproc <- function(adj, total_pop, init_plan = NULL, ndists = NULL,
                           pop_tol = NULL,
                           counties = NULL,
                           group_pop = NULL,
                           areasvec = NULL,
                           borderlength_mat = NULL, ssdmat = NULL,
                           compactness_metric = NULL,
                           partisan_metric = NULL,
                           temper = NULL, constraint = NULL,
                           constraintweights = NULL,
                           betaseq = NULL, betaseqlength = NULL,
                           betaweights = NULL, adjswaps = TRUE, maxiterrsg = NULL,
                           contiguitymap = "rooks", tgt_min = 0.55, tgt_other = 0.25,
                           rvote,
                           dvote,
                           minorityprop = NULL,
                           verbose = TRUE
){

  #########################
  ## Inputs to function: ##
  #########################
  ## adj - adjacency object of geographic units. Accepts adjlist or adjmat
  ## total_pop - population of each of the units
  ## init_plan - initial congressional units. Must be contiguous partitions. Default is NULL
  ## ndists - number of desired congressional units. Default is NULL
  ## pop_tol - strength of hard population constraint. Defaulted to no
  ##           constraint. pop_tol = 0.01 implies a 1% population constraint.
  ## group_pop - vector of populations for a minority group. To be used
  ##               in conjunction with the segregation and vra M-H constraints
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
  ## tgt_min - vra constraint larger grouppop decimal
  ## tgt_other - vra constraint smaller grouppop decimal
  ## rvote - republican (or Party A) votes by precinct
  ## dvote - democratic (or Party B!=A) votes by precinct
  ## partisan_metric - either efficiency-gap or proportional-representation
  ##
  #######################
  ## Check missingness ##
  #######################
  if(missing(adj)){
    stop("Please supply adjacency matrix or list")
  }
  if(missing(total_pop)){
    stop("Please supply vector of geographic unit populations")
  }
  if(!is.null(constraintweights)){
    if((any(constraintweights == 0) & !is.null(constraint))){
      stop("If applying constraints or using simulated tempering, please set non-zero constraint by specifying the 'constraintweight' argument, and specify the names of the constraints in 'constraint'.")
    }
    if(any(!(constraint %in% c("compact", "vra", "segregation", "population",
                               "similarity", "countysplit", "partisan", "minority", 'hinge')))){
      stop("Please specify any combination of `compact`, `segregation`, vra`, `population`, `countysplit`, `similarity`, `partisan`, `minority`, `hinge` for constraint")
    }
  }

  ############################################
  ## If not a list, convert adjlist to list ##
  ############################################
  if(!is.list(adj)){

    ## If a matrix, check to see if adjacency matrix
    if(is.matrix(adj)){

      ## Is it square?
      squaremat <- (nrow(adj) == ncol(adj))
      ## All binary entries?
      binary <- ((length(unique(c(adj))) == 2) &
                   (sum(unique(c(adj)) %in% c(0, 1)) == 2))
      ## Diagonal elements all 1?
      diag <- (sum(diag(adj)) == nrow(adj))
      ## Symmetric?
      symmetric <- isSymmetric(adj)

      ## If all are true, change to adjlist and automatically zero-index
      if(squaremat & binary & diag & symmetric){

        ## Initialize object
        adjlist <- vector("list", nrow(adj))

        ## Loop through rows in matrix
        for(i in 1:nrow(adj)){

          ## Extract row
          adjvec <- adj[,i]
          ## Find elements it is adjacent to
          inds <- which(adj == 1)
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
    }else if(class(adj) == "SpatialPolygonsDataFrame"){ ## shp object

      ## Convert shp object to adjacency list
      adjlist <- redist.adjacency(st_as_sf(adj))


    }else{ ## If neither list, matrix, or shp, throw error
      stop("Please input an adjacency list, adjacency matrix, or Spatial
                 Polygons shp file")
    }

  }else if('sf' %in% class(adj)){
    adjlist <- redist.adjacency(adj)
  }else{

    ## Rename adjacency object as list
    adjlist <- adj

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






  if(is.null(init_plan) || isTRUE( 'smc' %in% init_plan)){
      map = redist_map(pop=total_pop,  pop_tol=ifelse(is.null(pop_tol), 0.05, pop_tol),
                       ndists=ndists,  adj=adj)
    invisible(capture.output(
        init_plan <- redist_smc(map, nsims = 1, silent = TRUE),
        type = 'message'))
    init_plan <- as.matrix(init_plan)[, 1]
  } else if(!is.null(init_plan) && 'rsg' %in% init_plan) {
    ##############################################################################
    ## If no init_plan == rsg, use Random Seed and Grow                         ##
    ## (Chen and Rodden 2013) algorithm                                         ##
    ##############################################################################

    ## Set up target pop, strength of constraint (5%)
    if(is.null(pop_tol)){
      pop_tol_rsg <- .05
    }else{
      pop_tol_rsg <- pop_tol
    }

    ## Print start
    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    if(verbose){
    cat("\n", append = TRUE)
    cat(divider, append = TRUE)
    cat("Using redist.rsg() to generate starting values.\n\n", append= TRUE)
    }
    ## Run the algorithm
    initout <- redist.rsg(adj = adjlist,
                          total_pop = total_pop,
                          ndists = ndists,
                          pop_tol = pop_tol_rsg,
                          verbose = FALSE,
                          maxiter = maxiterrsg)
    ## Get initial cds
    init_plan <- initout$plan

  } else {
    ###################################################################
    ## Check whether initial partitions (if provided) are contiguous ##
    ###################################################################
    if(!is.na(init_plan)[1]){
      if(sum(is.na(init_plan)) > 0){
        stop("You have NA's in your congressional districts. Please check the provided init_plan vector for NA entries.")
      }

      ndists <- length(unique(init_plan))
      divlist <- genAlConn(adjlist, init_plan)
      ncontig <- countpartitions(divlist)

      if(ncontig != ndists){
        stop(paste("Your initial congressional districts have ", ndists,
                   " unique districts but ",
                   ncontig, " contigous connected components. Please provide a starting map with contigous districts.", sep = ""))
      }
    }
  }

  ###########################################################
  ## Check other inputs to make sure they are right length ##
  ###########################################################
  if((length(total_pop) != length(adjlist)) | (sum(is.na(total_pop)) > 0)){
    stop("Each entry in adjacency list must have a corresponding entry
              in vector of populations")
  }
  if((length(init_plan) != length(adjlist)) | (sum(is.na(init_plan)) > 0)){
    stop("Each entry in adjacency list must have an initial congressional
             district assignment")
  }
  if("segregation" %in% constraint & is.null(group_pop)){
    stop("If applying the segregation constraint, please provide a vector
             of subgroup populations")
  }
  if("vra" %in% constraint & is.null(group_pop)){
    stop("If applying the vra constraint, please provide a vector
             of subgroup populations")
  }
  if("countysplit" %in% constraint & is.null(counties)){
    stop("If applying the county split constraint, please provide a numeric vector indicating county membership.")
  }
  if("segregation" %in% constraint & !(is.null(group_pop))){
    if((length(group_pop) != length(adjlist)) |
       (sum(is.na(group_pop)) > 0)){
      stop("If applying the segregation constraint, each entry in adjacency
              list must have corresponding entry in vector of group populations")
    }
  }
  if("vra" %in% constraint & !(is.null(group_pop))){
    if((length(group_pop) != length(adjlist)) |
       (sum(is.na(group_pop)) > 0)){
      stop("If applying the vra constraint, each entry in adjacency
              list must have corresponding entry in vector of group populations")
    }
  }
  if("countysplit" %in% constraint & !is.null(counties)){
    if(length(counties) != length(adjlist) | sum(is.na(counties)) > 0){
      stop("You do not have a county membership assigned for every unit.")
    }
  }

  if("partisan" %in% constraint){
    if(is.null(rvote)){
      stop('You must provide an integer vector to rvote when using partisan constraint.')
    }
    if(is.null(dvote)){
      stop('You must provide an integer vector to dvote when using partisan constraint.')
    }
    if(length(rvote) != length(adjlist)){
      stop("rvote must be an integer vector with corresponding entry for each entry in adjacency list")
    }
    if(length(dvote) != length(adjlist)){
      stop("dvote must be an integer vector with corresponding entry for each entry in adjacency list")
    }
    if(!'integer' %in% class(rvote) & !'numeric' %in% class(rvote)){
      stop('rvote must be an integer vector.')
    }
    if(!'integer' %in% class(dvote) & !'numeric' %in% class(dvote)){
      stop('dvote must be an integer vector.')
    }
    if(is.null(partisan_metric)){
      # defaults to efficiency-gap
      partisan_metric <- 'efficiency-gap'
    } else{
      if(!partisan_metric %in% c("efficiency-gap", "proportional-representation")){
        stop("partisan_metric is only implemented for choices efficiency-gap and proportional-representation.")
      }
    }
  } else{
    rvote <- c(0L,1L)
    dvote <- c(1L,0L)
  }

  if("minority" %in% constraint){
    if(!"numeric" %in% class(minorityprop)){
      stop('"minorityprop" must be of type numeric.')
    }
    if(length(minorityprop) > ndists){
      stop('"minorityprop" has more entries than there will be districts.')
    }
  }

  if("hinge" %in% constraint){
    if(!"numeric" %in% class(minorityprop)){
      stop('"minorityprop" must be of type numeric.')
    }
    if(length(minorityprop) > ndists){
      stop('"minorityprop" has more entries than there will be districts.')
    }
  }


  if (!any(c('hinge', 'minority') %in% constraint)) {
    minorityprop = 0 #init so it won't get mad in swMH input
  }

  ####################
  ## Zero-index cds ##
  ####################
  if(min(init_plan) != 0){
    init_plan <- redist.sink.plan(init_plan) - 1
  }
  if(length(unique(init_plan)) != (max(init_plan) + 1)){
    stop("The district numbers in init_plan must be consecutive. The input to `init_plan` could not be transformed using `redist.sink.plan()`.")
  }

  ## ------------------------------------
  ## Check VRA targets if necessary
  ## ------------------------------------
  if("vra" %in% constraint){
    if(! 'numeric' %in% class(tgt_min)){
      stop("Need vra constraint tgt_min to be numeric.")
    }
    if(! 'numeric' %in% class(tgt_other)){
      stop("Need vra constraint tgt_other to be numeric.")
    }
    if(tgt_min < 0 | tgt_min > 1){
      stop("Need vra constraint 0 <= tgt_min <= 1.")
    }
    if(tgt_other < 0 | tgt_other > 1){
      stop("Need vra constraint 0 <= tgt_other <= 1.")
    }
  }

  ####################################################
  ## Calculate parity and population margin allowed ##
  ####################################################
  dists <- length(unique(init_plan))
  if(is.null(pop_tol)){
    pop_tol <- 100
  }

  #####################################
  ## Set group_pop if not provided ##
  #####################################
  if(is.null(group_pop)){
    group_pop <- total_pop
  }

  ## -------------------------------------
  ## Set county membership if not provided
  ## -------------------------------------
  if(is.null(counties)){
    counties <- c(0, 0, 0, 0)
  }else{
    if(is.factor(counties)){
      counties <- as.numeric(counties)
    } else {
      counties <- redist.county.id(counties)
    }
    if(is.character(counties)){
     stop("Vector of counties must be numeric or factor.")
    }
    counties <- counties - min(counties)
  }

  ################################
  ## Set ssdmat if not provided ##
  ################################
  if(is.null(ssdmat) & "compact" %in% constraint & "fryer-holden" %in% compactness_metric){
    if(class(adj) == "SpatialPolygonsDataFrame"){
      adj <- sf::st_as_sf(adj)
      centroids <- sf::st_coordinates(sf::st_centroid(adj))
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
  if("compact" %in% constraint & "polsby-popper" %in% compactness_metric){
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
  beta <- ifelse(is.null(constraint) | (temper %in% c(TRUE)), 0, 1)
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
  if("vra" %in% constraint){
    weightvra <- constraintweights[which(constraint == "vra")]
  }else{
    weightvra <- 0
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

  if("partisan" %in% constraint){
    weightpartisan <- constraintweights[which(constraint == "partisan")]
  }else{
    weightpartisan <- 0
  }
  if("minority" %in% constraint){
    weightminority <- constraintweights[which(constraint == "minority")]
  }else{
    weightminority <- 0
  }
  if("hinge" %in% constraint){
    weighthinge <- constraintweights[which(constraint == "hinge")]
  }else{
    weighthinge <- 0
  }

  ###################################
  ## Check if betaspacing provided ##
  ###################################
  if("tempering" %in% temperbeta){
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
      total_pop = total_pop,
      init_plan = init_plan,
      group_pop = group_pop,
      areasvec = areasvec,
      borderlength_mat = borderlength_mat,
      ssdmat = ssdmat,
      counties = counties
    ),
    params = list(
      pctdistparity = pop_tol,
      dists = dists,
      beta = beta,
      temperbeta = temperbeta,
      betaseq = betaseq,
      betaweights = betaweights,
      adjswaps = adjswaps,
      weightpop = weightpop,
      weightcompact = weightcompact,
      weightseg = weightseg,
      weightvra = weightvra,
      weightsimilar = weightsimilar,
      weightcountysplit = weightcountysplit,
      weightpartisan = weightpartisan,
      weightminority = weightminority,
      weighthinge = weighthinge,
      tgt_min = tgt_min,
      tgt_other = tgt_other,
      rvote = rvote,
      dvote = dvote,
      partisan_metric = partisan_metric,
      minorityprop = minorityprop
    )
  )

  class(preprocout) <- "redist"

  return(preprocout)

}

