redist.rsg <- function(adj.list,
                    population,
                    Ndistrict,
                    target.pop,
                    thresh
                    ) {

    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    cat("\n")
    cat(divider)
    cat("redist.rsg(): Automated Redistricting Starts\n\n")

    ## Main Call to Computation
    time <- system.time(ret <- .Call('redist_rsg',
                 PACKAGE = 'redist',
				adj.list,
				population,
				Ndistrict,
				target.pop,
				thresh
                 ))

    cat(paste("\n\t", Ndistrict, " districts built using ", length(adj.list), " precincts in ",
    	round(time[3], digits=2), " seconds...\n\n", sep = ""))

    return(ret)
}
