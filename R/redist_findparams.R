#############################################
## Author: Ben Fifield
## Institution: Princeton University
## Date Created: 2015/08/26
## Date Modified: 2015/08/26
## Purpose: Function to try parameters and
##          get estimates of performance
#############################################

run_sims <- function(i, params, adjobj, popvec, nsims, ndists, initcds,
                     ssdmat, grouppopvec, countymembership, names, maxiterrsg, report_all,
                     adapt_lambda, adapt_eprob,
                     nstartval_store, maxdist_startval, logarg){
    
    ## Get this iteration
    p_sub <- params %>% dplyr::slice(i)
    if(logarg){
        sink(paste0("log_", i, ".txt"))
        cat("Parameter Values:\n")
        cat(c(p_sub))
        cat("\n")
    }
    
    ## Set parameter values
    if("eprob" %in% names){
        eprob <- p_sub$eprob
    }else{
        eprob <- 0.05
    }

    if("lambda" %in% names){
        lambda <- p_sub$eprob
    }else{
        lambda <- 0
    }

    if(!("popcons" %in% names)){
        popcons <- 100
    }else{
        popcons <- p_sub$popcons
    }

    ## Set constraints
    constraintvec <- c()
    weightvec <- c()
    addweights <- FALSE
    if("weight_population" %in% names){
        constraintvec <- c(constraintvec, "population")
        weightvec <- c(weightvec, p_sub$weight_population)
        addweights <- TRUE
    }
    if("weight_compact" %in% names){
        constraintvec <- c(constraintvec, "compact")
        weightvec <- c(weightvec, p_sub$weight_compact)
        addweights <- TRUE
    }
    if("weight_segregation" %in% names){
        constraintvec <- c(constraintvec, "segregation")
        weightvec <- c(weightvec, p_sub$weight_segregation)
        addweights <- TRUE
    }
    if("weight_similarity" %in% names){
        constraintvec <- c(constraintvec, "similarity")
        weightvec <- c(weightvec, p_sub$weight_similarity)
        addweights <- TRUE
    }
    if("weight_countysplit" %in% names){
        constraintvec <- c(constraintvec, "countysplit")
        weightvec <- c(weightvec, p_sub$weight_countysplit)
        addweights <- TRUE
    }
    if(addweights == FALSE){
        constraintvec <- NULL
        weightvec <- NULL
    }

    ## Run siulations
    out <- redist.mcmc(adjobj = adjobj, popvec = popvec, nsims = nsims,
                       ndists = ndists, ssdmat = ssdmat,
                       grouppopvec = grouppopvec,
                       countymembership = countymembership, #ctk-cran-note
                       initcds = initcds, eprob = eprob, lambda = lambda,
                       popcons = popcons, 
                       constraint = constraintvec,
                       constraintweights = weightvec, #ctk-cran-note
                       maxiterrsg = maxiterrsg,
                       adapt_lambda = adapt_lambda,
                       adapt_eprob = adapt_eprob)
    if(adapt_eprob){
        final_eprob <- out$final_eprob
    }
    if(adapt_lambda){
        final_lambda <- out$final_lambda
    }

    ## Sample districts for use as starting values
    ## Divide equally by distance
    inds <- which(1 - out$distance_original < maxdist_startval)
    cuts <- c(0, round(quantile(1:length(inds), (1:nstartval_store)/nstartval_store)))
    if(length(inds) == 0){
        cat(paste0("No maps available under parameter set ", i, ".\n"))
        startval <- NULL
    }else{
        startval <- matrix(NA, nrow(out$partitions), nstartval_store)
        for(i in 1:nstartval_store){
            sub <- inds[inds > cuts[i] & inds <= cuts[i+1]]
            startval[,i] <- out$partitions[,sample(sub, 1)]
        }
        startval <- as.matrix(startval)
    }
    
    ## Get quantiles
    quant <- floor(nsims / 4)
    q1 <- 1:quant
    q2 <- (quant + 1):(2*quant)
    q3 <- (2*quant + 1):(3 * quant)
    q4 <- (3*quant + 1):nsims

    ## Check acceptance rate
    mh_acceptance <- round(sum(out$mhdecisions) / length(out$mhdecisions),
                           digits = 3)

    ## Check population parity
    pop_parity <- round(mean(out$distance_parity), digits = 3)
    med_pop_parity <- round(median(out$distance_parity), digits = 3)
    range_pop_parity <- round(range(out$distance_parity), digits = 3)
    q1_pop_median <- round(median(out$distance_parity[q1]), digits = 3)
    q2_pop_median <- round(median(out$distance_parity[q2]), digits = 3)
    q3_pop_median <- round(median(out$distance_parity[q3]), digits = 3)
    q4_pop_median <- round(median(out$distance_parity[q4]), digits = 3)

    ## Check distance to original
    dist_orig <- round(mean(out$distance_original), digits = 3)
    med_dist_orig <- round(median(out$distance_original), digits = 3)
    range_dist_orig <- round(range(out$distance_original), digits = 3)
    q1_dist_median <- round(median(out$distance_original[q1]), digits = 3)
    q2_dist_median <- round(median(out$distance_original[q2]), digits = 3)
    q3_dist_median <- round(median(out$distance_original[q3]), digits = 3)
    q4_dist_median <- round(median(out$distance_original[q4]), digits = 3)

    ## Share of counties split
    if(!is.null(countymembership)){
        ncounties_split <- unlist(lapply(1:nsims, function(x){
            cd_assign <- out$partitions[,x]
            return(sum(tapply(cd_assign, countymembership, function(y){ifelse(length(unique(y)) > 1, 1, 0)})))
        }))
        starting_county_split <- ncounties_split[1]
        q1_countysplit_median <- median(ncounties_split[q1])
        q2_countysplit_median <- median(ncounties_split[q2])
        q3_countysplit_median <- median(ncounties_split[q3])
        q4_countysplit_median <- median(ncounties_split[q4])
    }

    ## -----------------
    ## Report statistics
    ## -----------------
    out <- paste("## -------------------------------------\n",
                 "## -------------------------------------\n",
                 "## Parameter Values for Simulation", i, "\n")
    if(!adapt_eprob){
        out <- paste0(out, "## Edgecut probability = ", eprob, "\n")
    }else{
        out <- paste0(out, "## Final adaptive edgecut probability = ",
                      final_eprob, "\n")
    }
    if(!adapt_lambda){
        out <- paste0(out, "## Lambda = ", lambda, "\n")
    }else{
        out <- paste0(out, "## Final adaptive lambda = ", final_lambda, "\n")
    }
    if(popcons != 100){
        out <- paste0(out, "## Hard population constraint = ", popcons, "\n")
    }else{
        out <- paste0(out, "## No hard population constraint applied\n")
    }
    if(!is.null(constraintvec)){
        out <- paste0(out, "## Setting constraints on ", as.character(constraintvec), "\n",
                      "## Weights  = ", weightvec, "\n")
    }else{
        out <- paste0(out, "## Not setting any soft constraints\n")
    }
    out <- paste0(out,
                 "## -------------------------------------\n",
                 "## Diagnostics:\n",
                 "## Metropolis-Hastings Acceptance Ratio = ", mh_acceptance, "\n")
    if("population" %in% constraintvec | report_all == TRUE){
        out <- paste0(
            out, "## Mean population parity distance = ",
            pop_parity, "\n",
            "## Median population parity distance = ",
            med_pop_parity, "\n",
            "## Population parity range = ",
            paste(range_pop_parity, collapse = " "),
            "\n",
            "## MCMC Iteration quantiles of population parity median = ",
            paste(q1_pop_median, q2_pop_median,
                  q3_pop_median, q4_pop_median, sep = " "),
            "\n")
    }
    if("similarity" %in% constraintvec | report_all == TRUE){
        out <- paste0(
            out,
            "\n## Mean share of geographies equal to initial assignment = ",
            dist_orig, "\n",
            "## Median share of geographies equal to initial assignment = ",
            med_dist_orig, "\n",
            "## Range of share of geographies equal to initial assignment = ",
            paste(range_dist_orig, collapse = " "), "\n",
            "## MCMC Iteration quantiles of geography distance to initial assignment = ",
            paste(q1_dist_median, q2_dist_median,
                  q3_dist_median, q4_dist_median, sep = " "),
            "\n")
    }
    if(!is.null(countymembership)){
        if("countysplit" %in% constraintvec | report_all == TRUE){
            out <- paste0(
                out,
                "\n## Median number of counties split = ",
                mean(ncounties_split), "\n",
                "## Median number of counties split = ",
                median(ncounties_split), "\n",
                "## Range of number of counties split = ",
                paste(range(ncounties_split), collapse = " "), "\n",
                "## Initial number of counties split = ",
                starting_county_split, "\n",
                "## MCMC Iteration quantiles of number of counties split = ",
                paste(q1_countysplit_median, q2_countysplit_median, q3_countysplit_median, q4_countysplit_median, sep = " "), "\n"
            )
        }
    }
    out <- paste0(out,
                  "## -------------------------------------\n",
                  "## -------------------------------------\n\n")
    if(logarg){
        sink()
    }
    
    return(list(printout = out, startval = startval))

}

#' Run parameter testing for \code{redist.mcmc}
#'
#' \code{redist.findparams} is used to find optimal parameter values of
#' \code{redist.mcmc} for a given map.
#'
#' @usage redist.findparams(adjobj, popvec, nsims, ndists = NULL, initcds = NULL,
#' adapt_lambda = FALSE, adapt_eprob = FALSE,
#' params, ssdmat = NULL, grouppopvec = NULL, countymembership = NULL,
#' nstartval_store, maxdist_startval,
#' maxiterrsg = 5000, report_all = TRUE,
#' parallel = FALSE, nthreads = NULL, log = FALSE, verbose = TRUE)
#'
#' @param adjobj An adjacency matrix, list, or object of class
#' "SpatialPolygonsDataFrame."
#' @param popvec A vector containing the populations of each
#' geographic unit.
#' @param nsims The number of simulations run before a save point.
#' @param ndists The number of congressional districts.
#' The default is \code{NULL}.
#' @param initcds A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided, random
#' and contiguous congressional district assignments will be generated using \code{redist.rsg}.
#' @param adapt_lambda Whether to adaptively tune the lambda parameter so that the Metropolis-Hastings
#' acceptance probability falls between 20\% and 40\%. Default is FALSE.
#' @param adapt_eprob Whether to adaptively tune the edgecut probability parameter so that the
#' Metropolis-Hastings acceptance probability falls between 20\% and 40\%. Default is
#' FALSE.
#' @param params A matrix of parameter values to test, such as the output of
#' \code{expand.grid}. Parameters accepted for \code{params} include \code{eprob},
#' \code{lambda}, \code{popcons}, \code{beta}, and \code{constraint}.
#' @param ssdmat A matrix of squared distances between geographic
#' units. The default is \code{NULL}.
#' @param grouppopvec A vector of populations for some sub-group of
#' interest. The default is \code{NULL}.
#' @param countymembership A vector of county membership assignments. The default is \code{NULL}.
#' @param nstartval_store The number of maps to sample from the preprocessing chain
#' for use as starting values in future simulations. Default is 1.
#' @param maxdist_startval The maximum distance from the starting map that
#' sampled maps should be. Default is 100 (no restriction). 
#' @param maxiterrsg Maximum number of iterations for random seed-and-grow
#' algorithm to generate starting values. Default is 5000.
#' @param report_all Whether to report all summary statistics for each set of
#' parameter values. Default is \code{TRUE}.
#' @param parallel Whether to run separate parameter settings in parallel.
#' Default is \code{FALSE}.
#' @param nthreads Number of parallel tasks to run, declared outside of the
#' function. Default is \code{NULL}.
#' @param log Whether to open a log to track progress for each parameter combination
#' being tested. Default is FALSE.
#' @param verbose Whether to print additional information about the tests.
#' Default is \code{TRUE}.
#'
#' @details This function allows users to test multiple parameter settings of
#' \code{redist.mcmc} in preparation for a longer run for analysis.
#'
#' @return \code{redist.findparams} returns a print-out of summary statistics
#' about each parameter setting.
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#' 
#' @importFrom dplyr slice
#' 
#' @examples \dontrun{
#' data(algdat.pfull)
#'
#' ## Code to run the simulations in Figure 4 in Fifield, Higgins, Imai and
#' ## Tarr (2015)
#'
#' ## Get an initial partition
#' set.seed(1)
#' initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]
#'
#' params <- expand.grid(eprob = c(.01, .05, .1))
#'
#' ## Run the algorithm
#' redist.findparams(adjobj = algdat.pfull$adjlist,
#' popvec = algdat.pfull$precinct.data$pop,
#' initcds = initcds, nsims = 10000, params = params)
#' }
#' @export
redist.findparams <- function(adjobj, popvec, nsims, ndists = NULL, initcds = NULL,
                              adapt_lambda = FALSE, adapt_eprob = FALSE,
                              params, ssdmat = NULL, grouppopvec = NULL,
                              countymembership = NULL, 
                              nstartval_store = 1, maxdist_startval = 100,
                              maxiterrsg = 5000, report_all = TRUE,
                              parallel = FALSE, nthreads = NULL, log = FALSE, verbose = TRUE){

    ## Get number of trial parameter values to test
    trials <- nrow(params)

    ## Starting statement
    if(verbose){
        cat(paste("## ------------------------------\n",
            "## redist.findparams(): Parameter tuning for redist.mcmc()\n",
            "## Searching over", trials, "parameter combinations\n",
            "## ------------------------------\n\n", sep = " "))
    }

    ## Get parameters in params
    valid_names <- c("eprob", "lambda", "popcons", "weight_compact", "weight_population", "weight_segregation", "weight_similarity", "weight_countysplit")
    names <- names(params)
    if(sum(names %in% valid_names) < length(names)){
        invalid_name <- names[!(names %in% valid_names)]
        stop(paste(invalid_name, "is not a valid params input. Please see documentation.\n", sep = " "))
    }

    ## Check ndists, initcds
    if(is.null(ndists) & is.null(initcds)){
        stop("Please either supply a vector of starting congressional district assignments in `initcds' or a target number of congressional districts in `ndists`")
    }
    if(!is.null(initcds)){
        ndists <- length(unique(initcds))
    }
    if(sum("lambda" %in% names & adapt_lambda) > 0){
        warning("You have specified a grid of lambda values to search and set `adapt_lambda` to TRUE. Setting `adapt_lambda` to FALSE.")
        adapt_lambda <- FALSE
    }
    if(sum("eprob" %in% names & adapt_eprob) > 0){
        warning("You have specified a grid of eprob values to search and set `adapt_eprob` to TRUE. Setting `adapt_eprob` to FALSE.")
        adapt_eprob <- FALSE
    }    
    if("weight_segregation" %in% names & is.null(grouppopvec)){
        stop("If constraining on segregation, please provide a vector of group population.")
    }
    if("weight_compact" %in% names & is.null(ssdmat)){
        stop("If constraining on compactness, please provide a distances matrix.")
    }
    if("weight_similarity" %in% names & is.null(initcds)){
        stop("If constraining on similarity, please provide a vector of initial congressional district assignments.")
    }
    if("weight_countysplit" %in% names & is.null(countymembership)){
        stop("If constraining the number of county splits, please provide a vector of county assignments.")
    }

    if(parallel){ ## Parallel
        
        ## Check to see if threads declared
        if(is.null(nthreads)){
            stop("If parallelizing, please declare the number of threads")
        }

        ## Statement initializing parallelization
        if(verbose){
            cat(paste("## -----------------------------\n",
                "## Parallelizing over", nthreads, "processors\n",
                "## -----------------------------\n\n", sep = " "))
        }
        
        ## Set parallel environment
        if(verbose){
            cl <- makeCluster(nthreads, outfile = "")
        }else{
            cl <- makeCluster(nthreads)
        }
        registerDoParallel(cl)
        
        ## Execute foreach loop
        ret <- foreach(i = 1:trials, .verbose = verbose) %dopar% {

            ## Run simulations
            out <- run_sims(i, params, adjobj, popvec, nsims, ndists, initcds,
                            ssdmat, grouppopvec, countymembership, names, maxiterrsg, report_all,
                            adapt_lambda, adapt_eprob,
                            nstartval_store, maxdist_startval, log)
            
            ## Return values
            return(out)
            
        }

        printout <- c()
        startval <- vector(mode = "list", length = trials)
        for(i in 1:trials){
            printout <- paste(printout, ret[[i]]$printout)
            startval[[i]] <- ret[[i]]$startval
        }

    }else{ ## Sequential

        ## Create container for report
        printout <- c()
        startval <- vector(mode = "list", length = trials)
        
        ## Start loop over parameter values
        for(i in 1:trials){

            ## Run simulations
            out <- run_sims(i = i, params = params, adjobj, popvec, nsims, ndists, initcds,
                            ssdmat, grouppopvec, countymembership, names, maxiterrsg, report_all,
                            adapt_lambda, adapt_eprob,
                            nstartval_store, maxdist_startval, log)
            
            ## Add to printout
            printout <- paste(printout, out$printout)
            startval[[i]] <- out$startval
            
        }
        
    }

    if(parallel){
        stopCluster(cl)
    }

    cat(paste(printout, collapse = ""))
    return(list(diagnostics = printout, startvals = startval))

}

