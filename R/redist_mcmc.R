###########################################
## Author: Ben Fifield
## Institution: Princeton University
## Date Created: 2015/02/04
## Date Modified: 2015/03/09
## Purpose: R wrapper to run swMH() code (non-mpi)
###########################################


combine.par.anneal <- function(a, b) {

    .Deprecated("rbind.redist_plans()")
    ## Names of object
    name_out <- names(a)

    ## Create output object
    output_obj <- vector(mode = "list", length = length(a))

    ## Combine plans
    for (i in 1:length(a)) {
        if (i == i) {
            output_obj[[i]] <- cbind(a[[i]], b[[i]])
        } else {
            output_obj[[i]] <- c(a[[i]], b[[i]])
        }
    }

    names(output_obj) <- name_out
    return(output_obj)

}

#' (Deprecated) Flip MCMC Redistricting Simulator using Simulated Annealing
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
#' @param constraints A `redist_constr` list of constraints
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
#'
#' @return list of class redist
#'
#' @concept simulate
#' @export
redist.flip.anneal <- function(adj,
                               total_pop,
                               ndists = NULL,
                               init_plan = NULL,
                               constraints = redist_constr(),
                               num_hot_steps = 40000, num_annealing_steps = 60000,
                               num_cold_steps = 20000,
                               eprob = 0.05,
                               lambda = 0,
                               pop_tol = NULL,
                               rngseed = NULL, maxiterrsg = 5000,
                               adapt_lambda = FALSE, adapt_eprob = FALSE,
                               exact_mh = FALSE,
                               savename = NULL, verbose = TRUE) {
    .Deprecated("redist_flip_anneal", msg = "Please use `redist_flip_anneal. This is gone as of 4.1.")

    if (verbose) {
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
    if (missing(adj)) {
        stop("Please supply adjacency matrix or list")
    }
    if (missing(total_pop)) {
        stop("Please supply vector of geographic unit populations")
    }
    if (is.null(ndists) & is.null(init_plan) || is.null(ndists) & is.character(init_plan)) {
        stop("Please provide either the desired number of congressional districts
              or an initial set of congressional district assignments")
    }

    ## Set seed before first iteration of algorithm if provided by user
    if (!is.null(rngseed) & is.numeric(rngseed)) {
        set.seed(rngseed)
    }

    if (adapt_lambda) {
        adapt_lambda <- 1
    } else {
        adapt_lambda <- 0
    }
    if (adapt_eprob) {
        adapt_eprob <- 1
    } else {
        adapt_eprob <- 0
    }
    if (exact_mh) {
        exact_mh <- 1
    } else {
        exact_mh <- 0
    }

    ## ------------------
    ## Preprocessing data
    ## ------------------
    if (verbose) {
        cat("Preprocessing data.\n\n")
    }
    preprocout <- redist.preproc(adj = adj, total_pop = total_pop,
        init_plan = init_plan, ndists = ndists,
        pop_tol = pop_tol,
        temper = FALSE,
        betaseq = "powerlaw", betaseqlength = 10,
        betaweights = NULL,
        adjswaps = TRUE, maxiterrsg = maxiterrsg,
        verbose = verbose)


    if (verbose) {
        cat("Starting swMH().\n")
    }

    algout <- swMH(aList = preprocout$data$adjlist,
        cdvec = preprocout$data$init_plan,
        popvec = preprocout$data$total_pop,
        constraints = as.list(constraints),
        nsims = 100,
        eprob = eprob,
        pct_dist_parity = preprocout$params$pctdistparity,
        beta_sequence = preprocout$params$betaseq,
        beta_weights = preprocout$params$betaweights,
        lambda = lambda,
        beta = 0,
        adapt_beta = "annealing",
        adjswap = preprocout$params$adjswaps,
        exact_mh = exact_mh,
        adapt_lambda = adapt_lambda,
        adapt_eprob = adapt_eprob,
        num_hot_steps = num_hot_steps,
        num_annealing_steps = num_annealing_steps,
        num_cold_steps = num_cold_steps,
        verbose = as.logical(verbose))
    class(algout) <- "redist"

    ## -------------------------
    ## Combine and save the data
    ## -------------------------
    if (!is.null(savename)) {
        saveRDS(algout, file = paste0(savename, ".rds"))
    }

    ## Examine the data
    algout$plans <- algout$plans + 1
    return(algout)

}

#' (Deprecated) redist.combine.anneal
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
#' \item{constraint_qps}{A vector containing the value of the
#' QPS constraint for each accepted redistricting plan.}
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
redist.combine.anneal <- function(file_name) {

    .Deprecated("rbind.redist_plans")
    ## List files
    fn <- list.files()[grep(file_name, list.files())]
    if (length(fn) == 0) {
        stop("Can't find any files in current working directory with that name.")
    }
    load(fn[1])
    names_obj <- names(algout)

    # Create containers
    nr <- length(algout$plans)
    nc <- length(fn)
    plans <- matrix(NA, nrow = nr, ncol = nc)

    veclist <- vector(mode = "list", length = length(algout) - 1)
    for (i in 1:length(veclist)) {
        veclist[[i]] <- rep(NA, nc)
    }

    ## ------------
    ## Combine data
    ## ------------
    for (i in 1:length(fn)) {
        ## Load data
        load(fn[i])

        ## Store objects together
        for (j in 1:length(algout)) {
            if (j == 1) {
                plans[1:nr, i] <- algout$plans
            } else {
                veclist[[j - 1]][i] <- algout[[j]]
            }
        }
    }

    ## ---------------------------
    ## Store data in algout object
    ## ---------------------------
    algout <- vector(mode = "list", length = length(algout))
    for (i in 1:length(algout)) {
        if (i == 1) {
            algout[[i]] <- plans
        } else {
            algout[[i]] <- veclist[[i - 1]]
        }
    }
    names(algout) <- names_obj

    ## -------------
    ## Output object
    ## -------------
    class(algout) <- "redist"
    return(algout)

}


#' (Deprecated) Combine successive runs of \code{redist.flip}
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
#' \item{constraint_qps}{A vector containing the value of the
#' QPS constraint for each accepted redistricting plan.}
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
#' # alg_253 <- redist.flip(adj = fl25_adj, total_pop = fl25$pop,
#' # init_plan = init_plan, nsims = 10000,
#' # nloop = 2, savename = paste0(temp, "/test"))
#' # out <- redist.combine(savename = paste0(temp, "/test"), nloop = 2, nthin = 10)
#' }
#' @concept post
#' @export
redist.combine <- function(savename, nloop, nthin, temper = 0) {
    .Deprecated("rbind.redist_plans")
    ##############################
    ## Set up container objects ##
    ##############################
    algout <- readRDS(paste0(savename, "_loop1.rds"))
    names_obj <- names(algout)

    ## Create containers
    nr <- nrow(algout$plans)
    nc <- ncol(algout$plans)
    plans <- matrix(NA, nrow = nr,
        ncol = (nc*nloop/nthin))

    veclist <- vector(mode = "list", length = length(algout) - 1)
    for (i in 1:length(veclist)) {
        veclist[[i]] <- rep(NA, (nc*nloop/nthin))
    }

    ## Indices for thinning
    indthin <- which((1:nc) %% nthin == 0)

    ####################################
    ## Combine data in multiple loops ##
    ####################################
    for (i in 1:nloop) {

        ## Load data
        algout <- readRDS(paste0(savename, "_loop", i, ".rds"))

        ind <- ((i - 1)*(nc/nthin) + 1):(i*(nc/nthin))

        ## Store objects together
        for (j in 1:length(algout)) {
            if (j == 1) {
                plans[1:nr, ind] <- algout$plans[, indthin]
            } else {
                veclist[[j - 1]][ind] <- algout[[j]][indthin]
            }
        }

    }

    #################################
    ## Store data in algout object ##
    #################################
    algout <- vector(mode = "list", length = length(algout))
    for (i in 1:length(algout)) {
        if (i == 1) {
            algout[[i]] <- plans
        } else {
            algout[[i]] <- veclist[[i - 1]]
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

#' (Deprecated) Flip MCMC Redistricting Simulator
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
#' @param constraints A `redist_constr` list.
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
#' @param temper Whether to use simulated tempering algorithm. Default is FALSE.
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
#' \item{constraint_qps}{A vector containing the value of the
#' QPS constraint for each accepted redistricting plan.}
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
#'     init_plan = init_plan, nsims = 10000)
#'
#' ## You can also let it find a plan on its own!
#' sims <- redist.flip(adj = fl25_adj, total_pop = fl25$pop,
#'     ndists = 3, nsims = 10000)
#' }
#'
#' @concept simulate
#' @export
redist.flip <- function(adj,
                        total_pop,  nsims, ndists = NULL,
                        init_plan = NULL, constraints = redist_constr(),
                        loopscompleted = 0, nloop = 1,
                        warmup = 0, nthin = 1, eprob = 0.05,
                        lambda = 0,
                        pop_tol = NULL,
                        temper = FALSE,
                        betaseq = "powerlaw", betaseqlength = 10,
                        betaweights = NULL,
                        adjswaps = TRUE, rngseed = NULL, maxiterrsg = 5000,
                        adapt_lambda = FALSE, adapt_eprob = FALSE,
                        exact_mh = FALSE, savename = NULL,
                        verbose = TRUE) {
    .Deprecated("redist_flip", msg = "Please use `redist_flip`. This will be gone in 4.1.")

    if (verbose) {
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
    if (missing(adj)) {
        stop("Please supply adjacency matrix or list")
    }
    if (missing(total_pop)) {
        stop("Please supply vector of geographic unit populations")
    }
    if (missing(nsims)) {
        stop("Please supply number of simulations to run algorithm")
    }
    if (is.null(ndists) & is.null(init_plan) || is.null(ndists) & is.character(init_plan)) {
        stop("Please provide either the desired number of congressional districts
              or an initial set of congressional district assignments")
    }
    if (nloop > 1 & missing(savename)) {
        stop("Please supply save directory if saving simulations at checkpoints")
    }

    ## Set seed before first iteration of algorithm if provided by user
    if (!is.null(rngseed) & is.numeric(rngseed)) {
        set.seed(rngseed)
    }

    if (adapt_lambda) {
        adapt_lambda <- 1
    } else {
        adapt_lambda <- 0
    }
    if (adapt_eprob) {
        adapt_eprob <- 1
    } else {
        adapt_eprob <- 0
    }
    if (exact_mh) {
        exact_mh <- 1
    } else {
        exact_mh <- 0
    }

    #####################
    ## Preprocess data ##
    #####################
    if (verbose) {
        cat("Preprocessing data.\n\n")
    }
    preprocout <- redist.preproc(adj = adj,
        total_pop = total_pop,
        init_plan = init_plan,
        ndists = ndists,
        pop_tol = pop_tol,
        temper = temper,
        betaseq = betaseq,
        betaseqlength = betaseqlength,
        betaweights = betaweights,
        adjswaps = adjswaps,
        maxiterrsg = maxiterrsg,
        verbose = verbose
    )


    ## Get starting loop value
    loopstart <- loopscompleted + 1

    #######################
    ## Run the algorithm ##
    #######################
    for (i in loopstart:nloop) {
        ## Get congressional districts, tempered beta values
        if (i > loopstart) {
            cds <- algout$plans[, nsims]

            if (temper) {
                beta <- algout$beta_sequence[nsims]
            }

            if (!is.null(rngseed) & is.numeric(rngseed)) {
                set.seed(algout$randseed)
            }

            rm(list = "algout")

        } else {

            ## Reload the data if re-startomg
            if (loopstart > 1) {

                ## Load the data
                algout <- readRDS(paste0(savename, "_loop", i - 1, ".rds"))

                ## Stop if number of simulations per loop is different
                if (nsims != ncol(algout[[1]])) {
                    stop("Please specify the same number of simulations per
                     loop across all loops")
                }

                cds <- algout$plans[, nsims]

                if (temper) {
                    beta <- algout$beta_sequence[nsims]
                }

                if (!is.null(rngseed) & is.numeric(rngseed)) {
                    set.seed(algout$randseed)
                }

                rm(list = "algout")

            } else {
                cds <- preprocout$data$init_plan
            }

        }

        ## Run algorithm
        if (verbose) {
            cat("Starting swMH().\n")
        }
        algout <- swMH(aList = preprocout$data$adjlist,
            cdvec = cds,
            popvec = preprocout$data$total_pop,
            constraints = constraints,
            nsims = nsims*nthin + warmup,
            eprob = eprob,
            pct_dist_parity = preprocout$params$pctdistparity,
            beta_sequence = preprocout$params$betaseq,
            beta_weights = preprocout$params$betaweights,
            lambda = lambda,
            beta = preprocout$params$beta,
            adapt_beta = preprocout$params$temperbeta,
            adjswap = preprocout$params$adjswaps,
            exact_mh = exact_mh,
            adapt_lambda = adapt_lambda,
            adapt_eprob = adapt_eprob,
            verbose = as.logical(verbose))

        class(algout) <- "redist"

        ## Save random number state if setting the seed
        if (!is.null(rngseed)) {
            algout$randseed <- .Random.seed[3]
        }

        ## Save output
        if (nloop > 1) {
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
    if (nloop > 1) {
        redist.combine(savename = savename, nloop = nloop,
            nthin = nthin,
            temper = temperflag)
    } else if (!is.null(savename)) {
        saveRDS(algout, file = paste0(savename, ".rds"))
    }

    ## Examine the data
    if (nloop == 1) {
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
#' @param plans An object of class `redist_plans` from `redist_flip()`.
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
#' \item{constraint_qps}{A vector containing the value of the
#' QPS constraint for each accepted redistricting plan.}
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
#' map_ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
#' cons <- redist_constr(map_ia)
#' cons <- add_constr_pop_dev(cons, strength = 5.4)
#' alg <- redist_flip(map_ia, nsims = 500, constraints = cons)
#'
#' alg_ipw <- redist.ipw(plans = alg,
#'     resampleconstraint = "pop_dev",
#'     targetbeta = 1,
#'     targetpop = 0.05)
#' }
#'
#' @concept post
#' @export
redist.ipw <- function(plans,
                       resampleconstraint = c("pop_dev", "edges_removed",
                           "segregation", "status_quo"),
                       targetbeta,
                       targetpop = NULL,
                       temper = 0) {

    ## Warnings:
    if (missing(plans) | !inherits(plans, "redist_plans")) {
        cli_abort("Please provide {.arg plans} as a {.cls redist_plans}.")
    }

    plans_ref <- subset_ref(plans)
    plans <- subset_sampled(plans)

    if (length(resampleconstraint) != 1) {
        cli_abort("We currently only support one resamplingconstraint at a time.")
    }
    if (!(resampleconstraint %in% c("pop_dev", "edges_removed", "segregation", "status_quo"))) {
        cli_abort("We do not provide support for that constraint at this time")
    }
    if (missing(targetbeta)) {
        cli_abort("Please specify the target beta value")
    }

    ## Get indices drawn under target beta if tempering
    if (temper == 1) {
        indbeta <- which(plans$beta_sequence == targetbeta)
    } else {
        indbeta <- seq_len(ncol(get_plans_matrix(plans)))
    }

    ## Get indices of draws that meet target population
    if (!is.null(targetpop)) {
        indpop <- which(plans$distance_parity <= targetpop)
    } else {
        indpop <- seq_len(ncol(get_plans_matrix(plans)))
    }

    ## Get intersection of indices
    inds <- intersect(indpop, indbeta)
    ## Construct weights
    psi <- plans[[paste0("constraint_", resampleconstraint)]][inds]
    weights <- 1/exp(targetbeta*psi)

    ## Resample indices
    inds <- sample(inds, length(inds), replace = TRUE, prob = weights)
    ndists <- max(plans$district)
    indx <- unlist(lapply(inds, function(x) {seq(ndists*(x - 1) + 1, ndists*x, by = 1)}))

    ## Subset the entire list
    plans %>% slice(indx)
}

redist.warmup.chain <- function(algout, warmup = 1) {
    if (warmup <= 0) {
        return(algout)
    }
    inds <- 1:warmup
    algout_new <- vector(mode = "list", length = length(algout))
    for (i in 1:length(algout)) {

        ## Subset the matrix first, then the vectors
        if (i == 1) {
            algout_new[[i]] <- algout[[i]][, -inds]
        } else if (length(algout[[i]]) == 1) {
            algout_new[[i]] <- algout[[i]]
        } else if (all(names(algout[[i]]) == "adj")) {
            algout_new[[i]] <- algout[[i]]
        } else {
            algout_new[[i]] <- algout[[i]][-inds]
        }

    }
    names(algout_new) <- names(algout)
    class(algout_new) <- "redist"
    return(algout_new)
}


redist.thin.chain <- function(algout, thin = 100) {
    if (thin <= 1) {
        return(algout)
    }

    inds <- seq(1, ncol(algout$plans), by = thin)
    algout_new <- vector(mode = "list", length = length(algout))
    for (i in 1:length(algout)) {

        ## Subset the matrix first, then the vectors
        if (i == 1) {
            algout_new[[i]] <- algout[[i]][, inds]
        } else if (length(algout[[i]]) == 1) {
            algout_new[[i]] <- algout[[i]]
        } else if (!is.null(names(algout[[i]])) & all(names(algout[[i]]) == "adj")) {
            algout_new[[i]] <- algout[[i]]
        } else {
            algout_new[[i]] <- algout[[i]][inds]
        }

    }
    names(algout_new) <- names(algout)
    class(algout_new) <- "redist"
    return(algout_new)
}
