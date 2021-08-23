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

    ## Combine plans
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

#' Flip MCMC Redistricting Simulator using Simulated Annealing
#'
#' \code{redist.flip.anneal} simulates congressional redistricting plans
#' using Markov chain Monte Carlo methods coupled with simulated annealing.
#'
#' @param adj adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param total_pop A vector containing the populations of each geographic
#' unit
#' @param ndists The number of congressional districts. The default is
#' \code{NULL}.
#' @param init_plan A vector containing the congressional district labels
#' of each geographic unit. If not provided, random and contiguous congressional
#' district assignments will be generated using \code{redist_smc}. To use the old
#' behavior of generating with \code{redist.rsg}, provide init_plan = 'rsg'.
#' @param num_hot_steps The number of steps to run the simulator at beta = 0.
#' Default is 40000.
#' @param num_annealing_steps The number of steps to run the simulator with
#' linearly changing beta schedule. Default is 60000
#' @param num_cold_steps The number of steps to run the simulator at beta = 1.
#' Default is 20000.
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
#'
#' @return list of class redist
#'
#' @concept simulate
#' @export
redist.flip.anneal <- function(adj,
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
                               dvote = NULL, minorityprop = NULL){

    contiguitymap <- 'rooks'

    if(verbose){
        ## Initialize ##
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        cat("redist.flip.anneal(): Automated Redistricting Simulation Using
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
    if(is.null(ndists) & is.null(init_plan) || is.null(ndists) & is.character(init_plan) ){
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
    if(verbose){
        cat("Preprocessing data.\n\n")
    }
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
                                 partisan_metric = partisan_metric,
                                 verbose = verbose)


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

    if(verbose){
        cat("Starting swMH().\n")
    }

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
                   num_cold_steps = num_cold_steps,
                   verbose = as.logical(verbose))
    class(algout) <- "redist"

    ## -------------------------
    ## Combine and save the data
    ## -------------------------
    if(!is.null(savename)){
        saveRDS(algout, file = paste0(savename, ".rds"))
    }

    ## Examine the data
    algout$plans <- algout$plans + 1
    return(algout)

}

#' redist.combine.anneal
#'
#' Combine files generated by redist.flip.anneal()
#'
#' @usage redist.combine.anneal(file_name)
#'
#' @param file_name The file name to search for in current working directory.
#'
#' @return \code{redist.combine.anneal} returns an object of class "redist". The object
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
    nr <- length(algout$plans)
    nc <- length(fn)
    plans <- matrix(NA, nrow = nr, ncol = nc)

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
                plans[1:nr, i] <- algout$plans
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
            algout[[i]] <- plans
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


#' Combine successive runs of \code{redist.flip}
#'
#' \code{redist.combine} is used to combine successive runs of \code{redist.flip}
#' into a single data object
#'
#' @usage redist.combine(savename, nloop, nthin, temper)
#'
#' @param savename The name (without the loop or \code{.rds} suffix)
#' of the saved simulations.
#' @param nloop The number of loops being combined. Savename must be non-null.
#' @param nthin How much to thin the simulations being combined.
#' @param temper Wheterh simulated tempering was used (1) or not (0)
#' in the simulations. Default is 0.
#'
#' @details This function allows users to combine multiple successive runs of
#' \code{redist.flip} into a single \code{redist} object for analysis.
#'
#' @return \code{redist.combine} returns an object of class "redist". The object
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
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' @return  a redist object with entries combined
#'
#' @examples
#' \donttest{
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins,Imai and Tarr (2015)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## Run the algorithm
#' set.seed(1)
#' temp <- tempdir()
#' alg_253 <- redist.flip(adj = fl25_adj, total_pop = fl25$pop,
#'                        init_plan = init_plan, nsims = 10000,
#'                        nloop = 2, savename = paste0(temp, "/test"))
#' out <- redist.combine(savename = paste0(temp, "/test"), nloop = 2, nthin = 10)
#' }
#' @concept post
#' @export
redist.combine <- function(savename, nloop, nthin, temper = 0){
    ##############################
    ## Set up container objects ##
    ##############################
    algout <- readRDS(paste0(savename, "_loop1.rds"))
    names_obj <- names(algout)

    ## Create containers
    nr <- nrow(algout$plans)
    nc <- ncol(algout$plans)
    plans <- matrix(NA, nrow = nr,
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
        algout <- readRDS(paste0(savename, "_loop", i, ".rds"))

        ind <- ((i - 1) * (nc / nthin) + 1):(i * (nc / nthin))

        ## Store objects together
        for(j in 1:length(algout)){
            if(j == 1){
                plans[1:nr, ind] <- algout$plans[,indthin]
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
            algout[[i]] <- plans
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
    saveRDS(algout, file = paste0(savename, ".rds"))

    return(algout)
}

#' Flip MCMC Redistricting Simulator
#'
#' \code{redist.mcmc} is used to simulate Congressional redistricting
#' plans using Markov Chain Monte Carlo methods.
#'
#' @param adj adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param total_pop A vector containing the populations of each geographic
#' unit
#' @param nsims The number of simulations run before a save point.
#' @param ndists The number of congressional districts. The default is
#' \code{NULL}.
#' @param init_plan A vector containing the congressional district labels
#' of each geographic unit. If not provided, random and contiguous congressional
#' district assignments will be generated using \code{redist_smc}. To use the old
#' behavior of generating with \code{redist.rsg}, provide init_plan = 'rsg'.
#' @param loopscompleted Number of save points reached by the
#' algorithm. The default is \code{0}.
#' @param nloop The total number of save points for the algorithm. The
#' default is \code{1}. Note that the total number of simulations run
#' will be \code{nsims} * \code{nloop}. \code{savename} must be non-null.
#' @param warmup The number of warmup samples to discard. The default is 0.
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
#' \donttest{
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins, Imai and Tarr (2015)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## Run the algorithm
#' alg_253 <- redist.flip(adj = fl25_adj, total_pop = fl25$pop,
#'                        init_plan = init_plan, nsims = 10000)
#'
#'  ## You can also let it find a plan on its own!
#'  sims <- redist.flip(adj = fl25_adj, total_pop = fl25$pop,
#'                        ndists = 3, nsims = 10000)
#'
#'
#' }
#'
#' @concept simulate
#' @export
redist.flip <- function(adj,
                        total_pop,  nsims, ndists = NULL,
                        init_plan = NULL,
                        loopscompleted = 0, nloop = 1,
                        warmup = 0, nthin = 1, eprob = 0.05,
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
                        rvote = NULL, dvote = NULL, minorityprop = NULL){

    contiguitymap <- 'rooks'

    if(verbose){
        ## Initialize ##
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

        cat("\n", append = TRUE)
        cat(divider, append = TRUE)
        cat("redist.flip(): Automated Redistricting Simulation Using
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
    if(is.null(ndists) & is.null(init_plan) || is.null(ndists) & is.character(init_plan) ){
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
    if(verbose){
    cat("Preprocessing data.\n\n")
    }
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
                                 minorityprop = minorityprop,
                                 verbose = verbose
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
            cds <- algout$plans[,nsims]

            if(temper){
                beta <- algout$beta_sequence[nsims]
            }

            if(!is.null(rngseed) & is.numeric(rngseed)){
                set.seed(algout$randseed)
            }

            rm(list = "algout")

        } else{

            ## Reload the data if re-startomg
            if(loopstart > 1){

                ## Load the data
                algout <- readRDS(paste0(savename, "_loop", i - 1, ".rds"))

                ## Stop if number of simulations per loop is different
                if(nsims != ncol(algout[[1]])){
                    stop("Please specify the same number of simulations per
                     loop across all loops")
                }

                cds <- algout$plans[,nsims]

                if(temper){
                    beta <- algout$beta_sequence[nsims]
                }

                if(!is.null(rngseed) & is.numeric(rngseed)){
                    set.seed(algout$randseed)
                }

                rm(list = "algout")

            }else{
                cds <- preprocout$data$init_plan
            }

        }

        ## Run algorithm
        if(verbose){
            cat("Starting swMH().\n")
        }
        algout <- swMH(aList = preprocout$data$adjlist,
                       cdvec = cds,
                       cdorigvec = preprocout$data$init_plan,
                       popvec = preprocout$data$total_pop,
                       grouppopvec = preprocout$data$group_pop,
                       areas_vec = preprocout$data$areasvec,
                       county_membership = preprocout$data$counties,
                       borderlength_mat = preprocout$data$borderlength_mat,
                       nsims = nsims * nthin + warmup,
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
                       minorityprop = preprocout$params$minorityprop,
                       verbose = as.logical(verbose))

        class(algout) <- "redist"

        ## Save random number state if setting the seed
        if(!is.null(rngseed)){
            algout$randseed <- .Random.seed[3]
        }

        ## Save output
        if(nloop > 1){
            saveRDS(algout, file = paste0(savename, "_loop", i, ".rds"))
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
        saveRDS(algout, file = paste0(savename, ".rds"))
    }

    ## Examine the data
    if(nloop == 1){
        algout <- redist.warmup.chain(algout = algout, warmup = warmup)
        algout <- redist.thin.chain(algout, thin = nthin)
    }

    algout$plans <- algout$plans + 1

    return(algout)
}

#' Inverse probability reweighting for MCMC Redistricting
#'
#' \code{redist.ipw} properly weights and resamples simulated redistricting plans
#' so that the set of simulated plans resemble a random sample from the
#' underlying distribution. \code{redist.ipw} is used to correct the sample when
#' population parity, geographic compactness, or other constraints are
#' implemented.
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
#' @examples
#' \donttest{
#' data(iowa)
#' adj <- redist.adjacency(iowa)
#' init_plan <- iowa$cd_2010
#'alg <- redist.flip(adj = adj, total_pop = iowa$pop,
#'                   init_plan = init_plan, nsims = 1000,
#'                   constraint = 'population', constraintweights = 5.4)
#'
#' alg_ipw <- redist.ipw(algout = alg,
#'                      resampleconstraint = 'pop',
#'                      targetbeta = 1,
#'                      targetpop = 0.05)
#' }
#'
#' @concept post
#' @export
redist.ipw <- function(algout,
                       resampleconstraint = c("pop", "compact",
                                              "segregation", "similar"),
                       targetbeta,
                       targetpop = NULL,
                       temper = 0){

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
        indbeta <- 1:ncol(algout$plans)
    }

    ## Get indices of draws that meet target population
    if(!is.null(targetpop)){
        indpop <- which(algout$distance_parity <= targetpop)
    }else{
        indpop <- 1:ncol(algout$plans)
    }

    ## Get intersection of indices
    inds <- intersect(indpop, indbeta)
    ## Construct weights
    psi <- algout[[paste0("constraint_", resampleconstraint)]][inds]
    weights <- 1 / exp(targetbeta * psi)

    ## Resample indices
    inds <- sample(inds, length(inds), replace = TRUE, prob = weights)

    ## Subset the entire list
    algout_new <- vector(mode = "list", length = length(algout))
    for(i in 1:length(algout_new)){

        ## Subset the matrix first, then the vectors
        if(i == 1){
            algout_new[[i]] <- algout[[i]][,inds]
        } else if(length(algout[[i]]) == 1) {
            algout_new[[i]] <- algout[[i]]
        } else if(all(names(algout[[i]]) == 'adj')) {
            algout_new[[i]] <- algout[[i]]
        } else {
            algout_new[[i]] <- algout[[i]][inds]
        }
    }
    names(algout_new) <- names(algout)

    ## Change class
    class(algout_new) <- "redist"

    return(algout_new)

}

redist.warmup.chain <- function(algout, warmup = 1){
    if(warmup <= 0){
        return(algout)
    }
    inds <- 1:warmup
    algout_new <- vector(mode = "list", length = length(algout))
    for(i in 1:length(algout)){

        ## Subset the matrix first, then the vectors
        if(i == 1){
            algout_new[[i]] <- algout[[i]][,-inds]
        } else if(length(algout[[i]]) == 1) {
            algout_new[[i]] <- algout[[i]]
        } else if(all(names(algout[[i]]) == 'adj')) {
            algout_new[[i]] <- algout[[i]]
        } else {
            algout_new[[i]] <- algout[[i]][-inds]
        }

    }
    names(algout_new) <- names(algout)
    class(algout_new) <- "redist"
    return(algout_new)
}


redist.thin.chain <- function(algout, thin = 100){
    if(thin <= 1){
        return(algout)
    }

    inds <- seq(1, ncol(algout$plans), by = thin)
    algout_new <- vector(mode = "list", length = length(algout))
    for(i in 1:length(algout)){

        ## Subset the matrix first, then the vectors
        if(i == 1){
            algout_new[[i]] <- algout[[i]][,inds]
        } else if(length(algout[[i]]) == 1) {
            algout_new[[i]] <- algout[[i]]
        } else if(!is.null(names(algout[[i]])) & all(names(algout[[i]]) == 'adj')) {
            algout_new[[i]] <- algout[[i]]
        } else {
            algout_new[[i]] <- algout[[i]][inds]
        }

    }
    names(algout_new) <- names(algout)
    class(algout_new) <- "redist"
    return(algout_new)
}

