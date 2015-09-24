#############################################
## Author: Ben Fifield
## Institution: Princeton University
## Date Created: 2015/08/26
## Date Modified: 2015/08/26
## Purpose: Function to try parameters and
##          get estimates of performance
#############################################

run_sims <- function(i, params, adjobj, popvec, nsims, ndists, initcds,
                     ssdmat, grouppopvec, names, maxiterrsg, report_all){
    
    ## Get this iteration
    p_sub <- params[i,]

    ## Set parameter values
    if(!("eprob" %in% names)){
        eprob <- 0.05
    }else{
        eprob <- p_sub$eprob
    }

    if(!("lambda" %in% names)){
        lambda <- 0
    }else{
        lambda <- p_sub$lambda
    }

    if(!("popcons" %in% names)){
        popcons <- 100
    }else{
        popcons <- p_sub$popcons
    }

    if(!("beta" %in% names)){
        beta <- 0
    }else{
        beta <- p_sub$beta
    }

    if(!("constraint" %in% names)){
        constraint <- "none"
    }else{
        constraint <- p_sub$constraint
        if(!(constraint %in% c("compact", "segregation", "population",
                               "similarity"))){
            stop("Please select either `compact,` `segregation,` `population,` or `similarity` for constraint type")
        }
    }

    ## Warnings
    if(constraint == "segregation" & is.null(grouppopvec)){
        stop("If constraining on segregation, please provide a vector of group population")
    }
    if(constraint == "compact" & is.null(ssdmat)){
        stop("If constraining on compactness, please provide a distances matrix")
    }
    if(constraint == "similarity" & is.null(initcds)){
        stop("If constraining on similarity, please provide a vector of initial congressional district assignments")
    }

    ## Run siulations
    out <- redist.mcmc(adjobj = adjobj, popvec = popvec, nsims = nsims,
                       ndists = ndists, ssdmat = ssdmat,
                       grouppopvec = grouppopvec,
                       initcds = initcds, eprob = eprob, lambda = lambda,
                       popcons = popcons, beta = beta,
                       constraint = constraint,
                       maxiterrsg = maxiterrsg)

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

    ## Report statistics
    out <- paste("########################################\n",
                 "## Parameter Values for Simulation", i, "\n",
                 "## Edgecut probability =", eprob, "\n",
                 "## Lambda =", lambda, "\n")
    if(popcons != 100){
        out <- paste(out, "## Hard population constraint =", popcons, "\n",
                     sep = " ")
    }else{
        out <- paste(out, "## No hard population constraint applied\n",
                     sep = " ")
    }
    out <- paste(out, "## Soft constraint is", as.character(constraint), "\n",
                 "## Target beta =", beta, "\n",
                 "########################################\n",
                 "## Diagnostics:\n",
                 "## Metropolis-Hastings Acceptance Ratio =", mh_acceptance,
                 "\n", sep = " ")
    if(constraint == "population" | report_all == TRUE){
        out <- paste(out, "## Mean population parity distance =",
                     pop_parity, "\n",
                     "## Median population parity distance =",
                     med_pop_parity, "\n",
                     "## Population parity range =",
                     paste(range_pop_parity, collapse = " "),
                     "\n",
                     "## MCMC Iteration quantiles of population parity median =",
                     q1_pop_median, q2_pop_median, q3_pop_median, q4_pop_median, "\n",
                     sep = " ")
    }
    if(constraint == "similarity" | report_all == TRUE){
        out <- paste(out, "\n## Mean share of geographies equal to initial assignment =", dist_orig, "\n",
                     "## Median share of geographies equal to initial assignment =", med_dist_orig, "\n",
                     "## Range of share of geographies equal to initial assignment =", paste(range_dist_orig, collapse = " "), "\n", sep = " ")
    }
    out <- paste(out, "########################################\n\n",
                 sep = " ")

    return(out)

}

redist.findparams <- function(adjobj, popvec, nsims, ndists = NULL, initcds = NULL,
                              params, ssdmat = NULL, grouppopvec = NULL,
                              maxiterrsg = 5000, report_all = TRUE,
                              parallel = FALSE, nthreads = NULL, verbose = TRUE){

    ## Get number of trial parameter values to test
    trials <- nrow(params)

    ## Starting statement
    if(verbose){
        cat(paste("########################################\n",
            "## redist.findparams(): Parameter tuning for redist.mcmc()\n",
            "## Searching over", trials, "parameter combinations\n",
            "########################################\n\n", sep = " "))
    }

    ## Get parameters in params
    valid_names <- c("eprob", "lambda", "popcons", "beta", "constraint")
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

    if(parallel){ ## Parallel

        ## Check to see if threads declared
        if(is.null(nthreads)){
            stop("If parallelizing, please declare the number of threads")
        }

        ## Statement initializing parallelization
        if(verbose){
            cat(paste("########################################\n",
                "## Parallelizing over", nthreads, "processors\n",
                "########################################\n\n", sep = " "))
        }
        
        ## Set parallel environment
        if(verbose){
            cl <- makeCluster(nthreads, outfile = "")
        }else{
            cl <- makeCluster(nthreads)
        }
        registerDoParallel(cl)
        
        ## Execute foreach loop
        printout <- foreach(i = 1:trials, .combine = paste,
                            .verbose = verbose,
                            .export = c("params", "adjobj", "popvec", "nsims",
                                "ndists", "initcds", "ssdmat",
                                "grouppopvec", "names",
                                "maxiterrsg", "report_all", "run_sims")) %dopar% {

            ## Run simulations
            out <- run_sims(i, params, adjobj, popvec, nsims, ndists, initcds,
                            ssdmat, grouppopvec, names, maxiterrsg, report_all)
            
            ## Return values
            return(out)
            
        }

    }else{ ## Sequential

        ## Create container for report
        printout <- c()
        
        ## Start loop over parameter values
        for(i in 1:trials){

            ## Run simulations
            out <- run_sims(i, params, adjobj, popvec, nsims, ndists, initcds,
                            ssdmat, grouppopvec, names, maxiterrsg, report_all)
            
            ## Add to printout
            printout <- paste(printout, out)
            
        }
        
    }

    cat(printout)

    if(parallel){
        stopCluster(cl)
    }

}

