redist.rsg <- function(adj.list,
                       population,
                       ndists,
                       thresh,
                       verbose = TRUE,
                       maxiter=5000
                       )
{
    if(verbose){
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")
        
        cat("\n")
        cat(divider)
        cat("redist.rsg(): Automated Redistricting Starts\n\n")
    }

    target.pop <- sum(population) / ndists
    
    ## Main Call to Computation - if returning NA, break.
    ## If returning districts but not contiguous, repeat
    ## First attempt
    time <- system.time(ret <- .Call('redist_rsg',
                                     PACKAGE = 'redist',
                                     adj.list,
                                     population,
                                     ndists,
                                     target.pop,
                                     thresh,
                                     as.integer(maxiter)
                                     ))
    ## Make another call if stuck, but only do one more try
    ## because maxiter might be too low
    if(is.na(ret$district_membership[1])){
        time <- system.time(ret <- .Call('redist_rsg',
                                         PACKAGE = 'redist',
                                         adj.list,
                                         population,
                                         ndists,
                                         target.pop,
                                         thresh,
                                         as.integer(maxiter)
                                         ))
    }

    if(is.na(ret$district_membership[1])){

        stop("redist.rsg() failed to return a valid partition. Try increasing maxiterrsg")
        
    }

    if(verbose){
        cat(paste("\n\t", ndists, " districts built using ",
                  length(adj.list), " precincts in ",
                  round(time[3], digits=2), " seconds...\n\n", sep = ""), append = TRUE)
    }
    
    return(ret)

}

