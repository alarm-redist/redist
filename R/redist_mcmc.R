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

#' MCMC Redistricting Simulator using Simulated Annealing
#'
#' \code{redist.mcmc.anneal} simulates congressional redistricting plans
#' using Markov chain Monte Carlo methods coupled with simulated annealing.
#'
#' @param adj adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param total_pop A vector containing the populations of each geographic
#' unit
#' @param ndists The number of congressional districts. The default is
#' \code{NULL}.
#' @param init_plan A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided,
#' random and contiguous congressional district assignments will be generated
#' using \code{redist.rsg}.
#' @param num_hot_steps The number of steps to run the simulator at beta = 0.
#' Default is 40000.
#' @param num_annealing_steps The number of steps to run the simulator with
#' linearly changing beta schedule. Default is 60000
#' @param num_cold_steps The number of steps to run the simulator at beta = 1.
#' Default is 20000.
#' @param eprob The probability of keeping an edge connected. The
#' default is \code{0.05}.
#' @param lambda The parameter determining the number of swaps to attempt
#' each iteration of the algorihtm. The number of swaps each iteration is
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
#' @param constraint Which constraint to apply. Accepts any combination of \code{compact},
#' \code{segregation}, \code{vra}, \code{population}, \code{similarity}, \code{partisan},
#' \code{minority}, \code{hinge}, \code{countysplit}, or \code{none}
#' (no constraint applied). The default is NULL.
#' @param constraintweights The weights to apply to each constraint. Should be a vector
#' the same length as constraint. Default is NULL.
#' @param compactness_metric The compactness metric to use when constraining on
#' compactness. Default is \code{fryer-holden}, the other implemented options
#' are \code{polsby-popper} and \code{edges-removed}.
#' @param partisan_metric The partisan metric to use when constraining on partisan metrics.
#' Only implemented are "efficiency-gap" (default) and "proportional-representation".
#' @param rngseed Allows the user to set the seed for the
#' simulations. Default is \code{NULL}.
#' @param maxiterrsg Maximum number of iterations for random seed-and-grow
#' algorithm to generate starting values. Default is 5000.
#' @param adapt_lambda Whether to adaptively tune the lambda parameter so that the Metropolis-Hastings
#' acceptance probability falls between 20\% and 40\%. Default is FALSE.
#' @param adapt_eprob Whether to adaptively tune the edgecut probability parameter so that the
#' Metropolis-Hastings acceptance probability falls between 20\% and 40\%. Default is
#' FALSE.
#' @param exact_mh Whether to use the approximate (0) or exact (1)
#' Metropolis-Hastings ratio calculation for accept-reject rule. Default is FALSE.
#' @param savename Filename to save simulations. Default is \code{NULL}.
#' @param verbose Whether to print initialization statement.
#' Default is \code{TRUE}.
#' @param ncores The number of cores available to parallelize over. Default is 1.
#' @param tgt_min The majority minority target percent as a decimal. Default is 0.55.
#' @param tgt_other The remaining target percent as a decimal. Default is 0.25.
#' @param rvote integer vector of votes for Republicans by precinct
#' @param dvote integer vector of votes for Democrats by precinct
#' @param minorityprop numeric vector of targeted minority proportions for the top
#' districts with that proportion
#' @param adjobj Deprecated, use adj. An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param popvec A vector containing the populations of each geographic
#' unit
#' @param initcds Deprecated, use init_plan. A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided,
#' random and contiguous congressional district assignments will be generated
#' using \code{redist.rsg}.
#' @param popcons Deprecated, use pop_tol. The strength of the hard population
#' constraint. \code{popcons} = 0.05 means that any proposed swap that
#' brings a district more than 5\% away from population parity will be
#' rejected. The default is \code{NULL}.
#' @param grouppopvec Deprecated, use group_pop. A vector of populations for some sub-group of
#' interest. The default is \code{NULL}.
#' @param countymembership Deprecated, use counties. A vector of county membership assignments. The default is \code{NULL}.
#' @param contiguitymap Deprecated. Use queens or rooks distance criteria for generating an
#' adjacency list from a "SpatialPolygonsDataFrame" data type.
#' Default is "rooks".
#'
#' @concept simulate
#' @export
redist.mcmc.anneal <- function(adj,
                               total_pop,
                               ndists = NULL,
                               init_plan = NULL,
                               num_hot_steps = 40000, num_annealing_steps = 60000,
                               num_cold_steps = 20000,
                               eprob = 0.05,
                               lambda = 0,
                               pop_tol = NULL,
                               group_pop = NULL,
                               areasvec = NULL,
                               counties = NULL,
                               borderlength_mat = NULL,
                               ssdmat = NULL,
                               constraint = NULL, constraintweights = NULL,
                               compactness_metric = "fryer-holden",
                               partisan_metric = "efficiency-gap",
                               rngseed = NULL, maxiterrsg = 5000,
                               adapt_lambda = FALSE, adapt_eprob = FALSE,
                                exact_mh = FALSE,
                               savename = NULL, verbose = TRUE,
                               ncores = 1, tgt_min = 0.55, tgt_other = 0.25, rvote = NULL,
                               dvote = NULL, minorityprop = NULL,
                               adjobj, popvec, initcds, popcons, grouppopvec, countymembership,
                               contiguitymap){

    if(!missing(adjobj)){
        .Deprecated(new = 'adj', old = 'adjobj')
        adj <- adjobj
    }
    if(!missing(popvec)){
        .Deprecated(new = 'total_pop', old = 'popvec')
        total_pop <- popvec
    }
    if(!missing(initcds)){
        .Deprecated(new = 'init_plan', old = 'initcds')
        init_plan <- initcds
    }
    if(!missing(popcons)){
        .Deprecated(new = 'pop_tol', old = 'popcons')
        pop_tol <- popcons
    }
    if(!missing(grouppopvec)){
        .Deprecated(new = 'group_pop', old = grouppopvec)
        group_pop <- grouppopvec
    }
    if(!missing(countymembership)){
        .Deprecated(new = 'counties', old = 'countymembership')
        counties <- countymembership
    }

    if(!missing(contiguitymap)){
        .Deprecated(msg = 'contiguitymap has been deprecated. Rooks adjacency built with redist.adjacency().')
    }
    contiguitymap <- 'rooks'

    if(verbose){
        ## Initialize ##
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        cat("redist.mcmc.anneal(): Automated Redistricting Simulation Using
         Markov Chain Monte Carlo\n\n", append = TRUE)
    }

    ## --------------
    ## Initial checks
    ## --------------
    if(missing(adj)){
        stop("Please supply adjacency matrix or list")
    }
    if(missing(total_pop)){
        stop("Please supply vector of geographic unit populations")
    }
    if(is.null(ndists) & is.null(init_plan)){
        stop("Please provide either the desired number of congressional districts
              or an initial set of congressional district assignments")
    }
    if(!(contiguitymap %in% c("queens", "rooks"))){
        stop("Please supply `queens` or `rooks` for a distance criteria")
    }
    if(!is.null(constraint) & is.null(constraintweights)){
        stop("Please provide a weight value in 'constraintweights' for each constraint specified in 'constraint'.")
    }
    if(!(compactness_metric %in% c("fryer-holden", "polsby-popper", "edges-removed"))){
        stop("We only support 'fryer-holden', 'polsby-popper', or 'edges-removed' as compactness metrics.")
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

    ## ------------------
    ## Preprocessing data
    ## ------------------
    cat("Preprocessing data.\n\n")
    preprocout <- redist.preproc(adj = adj, total_pop = total_pop,
                                 init_plan = init_plan, ndists = ndists,
                                 pop_tol = pop_tol,
                                 counties = counties,
                                 group_pop = group_pop,
                                 areasvec = areasvec,
                                 borderlength_mat = borderlength_mat,
                                 ssdmat = ssdmat,
                                 compactness_metric = compactness_metric,
                                 temper = FALSE,
                                 constraint = constraint,
                                 constraintweights = constraintweights,
                                 betaseq = "powerlaw", betaseqlength = 10,
                                 betaweights = NULL,
                                 adjswaps = TRUE, maxiterrsg = maxiterrsg,
                                 contiguitymap = contiguitymap,
                                 tgt_min = tgt_min,
                                 tgt_other = tgt_other,
                                 rvote = rvote,
                                 dvote = dvote,
                                 minorityprop = minorityprop,
                                 partisan_metric = partisan_metric)


    ## Set betas - if tempering, modified later
    weightpop <- preprocout$params$weightpop
    weightcompact <- preprocout$params$weightcompact
    weightseg <- preprocout$params$weightseg
    weightvra <- preprocout$params$weightvra
    weightsimilar <- preprocout$params$weightsimilar
    weightcountysplit <- preprocout$params$weightcountysplit
    weightpartisan <- preprocout$params$weightpartisan
    weightminority <- preprocout$params$weightminority
    weighthinge <- preprocout$params$weighthinge

    cat("Starting swMH().\n")
    algout <- swMH(aList = preprocout$data$adjlist,
                   cdvec = preprocout$data$init_plan,
                   cdorigvec = preprocout$data$init_plan,
                   popvec = preprocout$data$total_pop,
                   grouppopvec = preprocout$data$group_pop,
                   areas_vec = preprocout$data$areasvec,
                   county_membership = preprocout$data$counties,
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
                   weight_vra = weightvra,
                   weight_similar = weightsimilar,
                   weight_countysplit = weightcountysplit,
                   weight_partisan = weightpartisan,
                   weight_minority = weightminority,
                   weight_hinge = weighthinge,
                   adapt_beta = "annealing",
                   adjswap = preprocout$params$adjswaps,
                   exact_mh = exact_mh,
                   adapt_lambda = adapt_lambda,
                   adapt_eprob = adapt_eprob,
                   compactness_measure = compactness_metric,
                   partisan_measure = preprocout$params$partisan_metricpartisan_metric,
                   tgt_min = tgt_min,
                   tgt_other = tgt_other,
                   rvote = preprocout$params$rvote,
                   dvote = preprocout$params$dvote,
                   minorityprop = preprocout$params$minorityprop,
                   num_hot_steps = num_hot_steps,
                   num_annealing_steps = num_annealing_steps,
                   num_cold_steps = num_cold_steps)
    class(algout) <- "redist"

    ## -------------------------
    ## Combine and save the data
    ## -------------------------
    if(!is.null(savename)){
        save(algout, file = paste(savename, ".RData", sep = ""))
    }

    ## Examine the data
    return(algout)

}

#' redist.combine.anneal
#'
#' Combine files generated by redist.mcmc.anneal()
#'
#' @usage redist.combine.anneal(file_name)
#'
#' @param file_name The file name to search for in current working directory.
#'
#' @concept post
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
#' \item{constraint_segregation}{A vector containing the value of the
#' segregation constraint for each accepted redistricting plan.}
#' \item{constraint_vra}{A vector containing the value of the
#' vra constraint for each accepted redistricting plan.}
#' \item{constraint_similar}{A vector containing the value of the similarity
#' constraint for each accepted redistricting plan.}
#' \item{constraint_partisan}{A vector containing the value of the
#' partisan constraint for each accepted redistricting plan.}
#' \item{constraint_minority}{A vector containing the value of the
#' minority constraint for each accepted redistricting plan.}
#' \item{constraint_hinge}{A vector containing the value of the
#' hinge constraint for each accepted redistricting plan.}
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
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins,Imai and
#' Tarr (2015)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## Run the algorithm
#' alg_253 <- redist.mcmc(adj = fl25_adj, total_pop = fl25$pop,
#'                        init_plan = init_plan, nsims = 10000,
#'                        loops = 2, savename = "test")
#' out <- redist.combine(savename = "test", nloop = 2, nthin = 10)
#' }
#' @concept post
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

#' MCMC Redistricting Simulator
#'
#' \code{redist.mcmc} is used to simulate Congressional redistricting
#' plans using Markov Chain Monte Carlo methods.
#'
#' @param adj adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param total_pop A vector containing the populations of each geographic
#' unit
#'
#' @param nsims The number of simulations run before a save point.
#' @param ndists The number of congressional districts. The default is
#' \code{NULL}.
#' @param init_plan A vector containing the congressional district labels
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
#' @param temper Whether to use simulated tempering algorithm. Default is FALSE.
#' @param constraint Which constraint to apply. Accepts any combination of \code{compact},
#' \code{segregation}, \code{vra}, \code{population}, \code{similarity}, \code{partisan},
#' \code{minority}, \code{hinge}, \code{countysplit}, or \code{none}
#' (no constraint applied). The default is NULL.
#' @param constraintweights The weights to apply to each constraint. Should be a vector
#' the same length as constraint. Default is NULL.
#' @param compactness_metric The compactness metric to use when constraining on
#' compactness. Default is \code{fryer-holden}, the other implemented options
#' are \code{polsby-popper} and \code{edges-removed}.
#' @param partisan_metric The partisan metric to use when constraining on partisan metrics.
#' Only implemented is "efficiency-gap", the default.
#' @param ssd_denom The normalizing constant for the sum-of-squared distance Fryer-Holden metric.
#' Default is 1.0 (unnormalized).
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
#' @param maxiterrsg Maximum number of iterations for random seed-and-grow
#' algorithm to generate starting values. Default is 5000.
#' @param adapt_lambda Whether to adaptively tune the lambda parameter so that the Metropolis-Hastings
#' acceptance probability falls between 20\% and 40\%. Default is FALSE.
#' @param adapt_eprob Whether to adaptively tune the edgecut probability parameter so that the
#' Metropolis-Hastings acceptance probability falls between 20\% and 40\%. Default is
#' FALSE.
#' @param exact_mh Whether to use the approximate (0) or exact (1)
#' Metropolis-Hastings ratio calculation for accept-reject rule. Default is FALSE.
#' @param savename Filename to save simulations. Default is \code{NULL}.
#' @param verbose Whether to print initialization statement.
#' Default is \code{TRUE}.
#' @param tgt_min The majority minority target percent as a decimal. Default is 0.55.
#' @param tgt_other The remaining target percent as a decimal. Default is 0.25.
#' @param rvote integer vector of votes for Republicans by precinct
#' @param dvote integer vector of votes for Democrats by precinct
#' @param minorityprop numeric vector of targeted minority proportions for the top
#' districts with that proportion
#' @param adjobj Deprecated, use adj. An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param initcds Deprecated, use init_plan. A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided,
#' random and contiguous congressional district assignments will be generated
#' using \code{redist.rsg}.
#' @param popvec Deprecated, use total_pop. A vector containing the populations of each geographic
#' unit
#' @param popcons Deprecated, use pop_tol. The strength of the hard population
#' constraint. \code{popcons} = 0.05 means that any proposed swap that
#' brings a district more than 5\% away from population parity will be
#' rejected. The default is \code{NULL}.
#' @param grouppopvec Deprecated, use group_pop. A vector of populations for some sub-group of
#' interest. The default is \code{NULL}.
#' @param countymembership Deprecated, use counties. A vector of county membership assignments. The default is \code{NULL}.
#' @param contiguitymap Deprecated. Use queens or rooks distance criteria for generating an
#' adjacency list from a "SpatialPolygonsDataFrame" data type.
#' Default is "rooks".
#'
#' @details This function allows users to simulate redistricting plans
#' using Markov Chain Monte Carlo methods. Several constraints
#' corresponding to substantive requirements in the redistricting process
#' are implemented, including population parity and geographic
#' compactness. In addition, the function includes multiple-swap and
#' simulated tempering functionality to improve the mixing of the Markov
#' Chain.
#'
#' @return \code{redist.mcmc} returns an object of class "redist". The object
#' \code{redist} is a list that contains the following components (the
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
#' \item{constraint_segregation}{A vector containing the value of the
#' segregation constraint for each accepted redistricting plan.}
#' \item{constraint_vra}{A vector containing the value of the
#' vra constraint for each accepted redistricting plan.}
#' \item{constraint_similar}{A vector containing the value of the similarity
#' constraint for each accepted redistricting plan.}
#' \item{constraint_partisan}{A vector containing the value of the
#' partisan constraint for each accepted redistricting plan.}
#' \item{constraint_minority}{A vector containing the value of the
#' minority constraint for each accepted redistricting plan.}
#' \item{constraint_hinge}{A vector containing the value of the
#' hinge constraint for each accepted redistricting plan.}
#' \item{beta_sequence}{A vector containing the value of beta for each iteration
#' of the algorithm. Returned when tempering is being used.}
#' \item{mhdecisions_beta}{A vector specifying whether a proposed beta value was
#' accepted (1) or rejected (0) in a given iteration of the algorithm. Returned
#' when tempering is being used.}
#' \item{mhprob_beta}{A vector containing the Metropolis-Hastings acceptance
#' probability for each iteration of the algorithm. Returned when tempering
#' is being used.}
#'
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' @examples
#' \dontrun{
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins,Imai and
#' Tarr (2015)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## Run the algorithm
#' alg_253 <- redist.mcmc(adj = fl25_adj, total_pop = fl25$pop,
#'                        init_plan = init_plan, nsims = 10000,
#'                        loops = 2, savename = "test")
#' }
#' @concept simulate
#' @export
redist.mcmc <- function(adj,
                        total_pop,  nsims, ndists = NULL,
                        init_plan = NULL,
                        loopscompleted = 0, nloop = 1, nthin = 1, eprob = 0.05,
                        lambda = 0,
                        pop_tol = NULL,
                        group_pop = NULL,
                        areasvec = NULL,
                        counties = NULL,
                        borderlength_mat = NULL, ssdmat = NULL, temper = FALSE,
                        constraint = NULL, constraintweights = NULL,
                        compactness_metric = "fryer-holden",
                        partisan_metric = "efficiency-gap",
                        ssd_denom = 1.0,
                        betaseq = "powerlaw", betaseqlength = 10,
                        betaweights = NULL,
                        adjswaps = TRUE, rngseed = NULL, maxiterrsg = 5000,
                        adapt_lambda = FALSE, adapt_eprob = FALSE,
                        exact_mh = FALSE, savename = NULL,
                        verbose = TRUE, tgt_min = 0.55, tgt_other = 0.25,
                        rvote = NULL, dvote = NULL, minorityprop = NULL,
                        adjobj, popvec, initcds, popcons, grouppopvec, countymembership,
                        contiguitymap){

    if(!missing(adjobj)){
        .Deprecated(new = 'adj', old = 'adjobj')
        adj <- adjobj
    }
    if(!missing(popvec)){
        .Deprecated(new = 'total_pop', old = 'popvec')
        total_pop <- popvec
    }
    if(!missing(initcds)){
        .Deprecated(new = 'init_plan', old = 'initcds')
        init_plan <- initcds
    }
    if(!missing(popcons)){
        .Deprecated(new = 'pop_tol', old = 'popcons')
        pop_tol <- popcons
    }
    if(!missing(grouppopvec)){
        .Deprecated(new = 'group_pop', old = grouppopvec)
        group_pop <- grouppopvec
    }
    if(!missing(countymembership)){
        .Deprecated(new = 'counties', old = 'countymembership')
        counties <- countymembership
    }

    if(!missing(contiguitymap)){
        .Deprecated(msg = 'contiguitymap has been deprecated. Rooks adjacency built with redist.adjacency().')
    }
    contiguitymap <- 'rooks'



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
    if(missing(adj)){
        stop("Please supply adjacency matrix or list")
    }
    if(missing(total_pop)){
        stop("Please supply vector of geographic unit populations")
    }
    if(missing(nsims)){
        stop("Please supply number of simulations to run algorithm")
    }
    if(is.null(ndists) & is.null(init_plan)){
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
    if(!(compactness_metric %in% c("fryer-holden", "polsby-popper", 'edges-removed'))){
        stop("We only support 'fryer-holden', 'polsby-popper', or 'edges-removed' as compactness metrics.")
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
    cat("Preprocessing data.\n\n")
    preprocout <- redist.preproc(adj = adj, 
                                 total_pop = total_pop,
                                 init_plan = init_plan, 
                                 ndists = ndists,
                                 pop_tol = pop_tol,
                                 counties = counties,
                                 group_pop = group_pop,
                                 areasvec = areasvec,
                                 borderlength_mat = borderlength_mat,
                                 ssdmat = ssdmat,
                                 compactness_metric = compactness_metric,
                                 partisan_metric = partisan_metric,
                                 temper = temper,
                                 constraint = constraint, 
                                 constraintweights = constraintweights,
                                 betaseq = betaseq, 
                                 betaseqlength = betaseqlength,
                                 betaweights = betaweights,
                                 adjswaps = adjswaps, 
                                 maxiterrsg = maxiterrsg,
                                 contiguitymap = contiguitymap,
                                 tgt_min = tgt_min,
                                 tgt_other = tgt_other,
                                 rvote = rvote,
                                 dvote = dvote,
                                 minorityprop = minorityprop
                                 )

    ## Set betas - if tempering, modified later
    weightpop <- preprocout$params$weightpop
    weightcompact <- preprocout$params$weightcompact
    weightseg <- preprocout$params$weightseg
    weightvra <- preprocout$params$weightvra
    weightsimilar <- preprocout$params$weightsimilar
    weightcountysplit <- preprocout$params$weightcountysplit
    weightpartisan <- preprocout$params$weightpartisan
    weightminority <- preprocout$params$weightminority
    weighthinge <- preprocout$params$weighthinge

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
                cds <- preprocout$data$init_plan
            }

        }

        ## Run algorithm
        cat("Starting swMH().\n")
        algout <- swMH(aList = preprocout$data$adjlist,
                       cdvec = cds,
                       cdorigvec = preprocout$data$init_plan,
                       popvec = preprocout$data$total_pop,
                       grouppopvec = preprocout$data$group_pop,
                       areas_vec = preprocout$data$areasvec,
                       county_membership = preprocout$data$counties,
                       borderlength_mat = preprocout$data$borderlength_mat,
                       nsims = nsims,
                       eprob = eprob,
                       pct_dist_parity = preprocout$params$pctdistparity,
                       beta_sequence = preprocout$params$betaseq,
                       beta_weights = preprocout$params$betaweights,
                       ssdmat = preprocout$data$ssdmat,
                       lambda = lambda,
                       beta = preprocout$params$beta,
                       weight_population = weightpop,
                       weight_compact = weightcompact,
                       weight_segregation = weightseg,
                       weight_vra = weightvra,
                       weight_similar = weightsimilar,
                       weight_countysplit = weightcountysplit,
                       weight_partisan = weightpartisan,
                       weight_minority = weightminority,
                       weight_hinge = weighthinge,
                       adapt_beta = preprocout$params$temperbeta,
                       adjswap = preprocout$params$adjswaps,
                       exact_mh = exact_mh,
                       adapt_lambda = adapt_lambda,
                       adapt_eprob = adapt_eprob,
                       compactness_measure = compactness_metric,
                       partisan_measure = preprocout$params$partisan_metric,
                       ssd_denom = ssd_denom,
                       tgt_min = tgt_min,
                       tgt_other = tgt_other,
                       rvote = preprocout$params$rvote,
                       dvote = preprocout$params$dvote,
                       minorityprop = preprocout$params$minorityprop)

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
    if(nloop == 1){
        return(algout)
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
#' @usage redist.ipw(algout,
#' resampleconstraint = c("pop", "compact", "segregation", "similar"),
#' targetbeta, targetpop = NULL, temper = 0)
#'
#' @param algout An object of class "redist".
#' @param resampleconstraint The constraint implemented in the simulations: one
#' of "pop", "compact", "segregation", or "similar".
#' @param targetbeta The target value of the constraint.
#' @param targetpop The desired level of population parity. \code{targetpop} =
#' 0.01 means that the desired distance from population parity is 1\%. The
#' default is \code{NULL}.
#' @param temper A flag for whether simulated tempering was used to improve the
#' mixing of the Markov Chain. The default is \code{1}.
#'
#' @details This function allows users to resample redistricting plans using
#' inverse probability weighting techniques described in Rubin (1987). This
#' techniques reweights and resamples redistricting plans so that the resulting
#' sample is representative of a random sample from the uniform distribution.
#'
#' @return \code{redist.ipw} returns an object of class "redist". The object
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
#' \item{constraint_segregation}{A vector containing the value of the
#' segregation constraint for each accepted redistricting plan.}
#' \item{constraint_similar}{A vector containing the value of the similarity
#' constraint for each accepted redistricting plan.}
#' \item{constraint_vra}{A vector containing the value of the
#' vra constraint for each accepted redistricting plan.}
#' \item{constraint_partisan}{A vector containing the value of the
#' partisan constraint for each accepted redistricting plan.}
#' \item{constraint_minority}{A vector containing the value of the
#' minority constraint for each accepted redistricting plan.}
#' \item{constraint_hinge}{A vector containing the value of the
#' hinge constraint for each accepted redistricting plan.}
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
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Code to run the simulations in Figure 4 of Fifield, Higgins,
#' ## Imai and Tarr (2015)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## Vector of beta weights
#' betaweights <- rep(NA, 10); for(i in 1:10){betaweights[i] <- 4^i}
#'
#' ## Run simulations - tempering population constraint
#' alg_253_20_st <- redist.mcmc(adj = fl25_adj, total_pop = fl25$pop,
#'                              init_plan = init_plan, nsims = 10000,
#'                              constraint = 'population', constraintweights = 5.4,
#'                              betaweights = betaweights, temper = 1)
#'
#' ## Resample using inverse probability weighting.
#' ## Target distance from parity is 20%
#' alg_253_20_st <- redist.ipw(alg_253_20_st, resampleconstraint = "pop",
#'                             targetbeta = 1, targetpop = .2, temper = 1)
#' }
#' @concept post
#' @export
redist.ipw <- function(algout,
                       resampleconstraint = c("pop", "compact",
                                              "segregation", "similar"),
                       targetbeta,
                       targetpop = NULL,
                       temper = 0)
{

    ## Warnings:
    if(missing(algout) | class(algout) != "redist"){
        stop("Please provide a proper redist object")
    }
    if(length(resampleconstraint)!=1){
        stop("We currently only support one resamplingconstraint at a time.")
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

####################################
# S3 generics

#' Extract the redistricting matrix from a \code{redist} object
#' @method as.matrix redist
#' @param x redist object
#' @param \dots additional arguments
#' @export
as.matrix.redist = function(x, ...) {
    x$plans
}

#' @method print redist
#' @importFrom utils str
#' @export
print.redist = function(x, ...) {
    cat(x$nsims, " sampled plans with ",
        ifelse(x$algorithm == 'mcmc', max(x$plans[,1]) + 1, max(x$plans[,1])),
        " districts from a ",
        length(x$adj), "-unit map, drawn\n using ",
        c(mcmc="Markov chain Monte Carlo",
          smc="Sequential Monte Carlo")[x$algorithm], sep="")
    if (x$pct_dist_parity < 1)
        cat(" and a ", 100*x$pct_dist_parity, "% population constraint.\n", sep="")
    else
        cat(".\n")
    cat(str(x$plans))
}
