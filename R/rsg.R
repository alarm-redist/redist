redist.rsg <- function(adj.list,
                       population,
                       ndists,
                       thresh,
                       verbose = TRUE
                       )
{

    if(verbose){
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")
        
        cat("\n")
        cat(divider)
        cat("redist.rsg(): Automated Redistricting Starts\n\n")
    }
    
    target.pop <- sum(population) / ndists
    
    ## Main Call to Computation
    time <- system.time(ret <- .Call('redist_rsg',
                                     PACKAGE = 'redist',
                                     adj.list,
                                     population,
                                     ndists,
                                     target.pop,
                                     thresh
                                     ))

    if(verbose){
        cat(paste("\n\t", ndists, " districts built using ", length(adj.list), " precincts in ",
                  round(time[3], digits=2), " seconds...\n\n", sep = ""))
    }
    
    return(ret)

}

