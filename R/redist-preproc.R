redist.preproc <- function(adj, total_pop, init_plan = NULL, ndists = NULL,
                           pop_tol = NULL,
                           temper = NULL,
                           betaseq = NULL, betaseqlength = NULL,
                           betaweights = NULL, adjswaps = TRUE, maxiterrsg = NULL,
                           verbose = TRUE
) {

    #########################
    ## Inputs to function: ##
    #########################
    ## adj - adjacency object of geographic units. Accepts adjlist or adjmat
    ## total_pop - population of each of the units
    ## init_plan - initial congressional units. Must be contiguous partitions. Default is NULL
    ## ndists - number of desired congressional units. Default is NULL
    ## pop_tol - strength of hard population constraint. Defaulted to no
    ## constraint. pop_tol = 0.01 implies a 1% population constraint.
    ## group_pop - vector of populations for a minority group. To be used
    ## in conjunction with the segregation and vra M-H constraints
    ## ssdmat - matrix of squared distances between population units.
    ## To be used when applying the compactness constraint.
    ## beta - target strength of constraint in MH ratio. Defaults to 0.
    ## temper - whether to use tempering (parallel or simulated) algorithms.
    ## Defaults to `none` (no tempering)
    ## betaseq - Spacing for beta sequence if tempering. Default is power law
    ## spacing, but can also be provided by user
    ## betaseqlength - Number of temperatures in the beta sequence. Default is
    ## ten
    ## betaweights - Vector of weights for beta sequence. Provided by user
    ## adjswaps - Flag for adjacent swaps for geyer-thompson tempering or MPI
    ## parallel tempering. Default to TRUE
    ## maxiterrsg - Maximum number of iterations for RSG algorithm
    ##
    #######################
    ## Check missingness ##
    #######################
    if (missing(adj)) {
        stop("Please supply adjacency matrix or list")
    }
    if (missing(total_pop)) {
        stop("Please supply vector of geographic unit populations")
    }

    ############################################
    ## If not a list, convert adjlist to list ##
    ############################################
    if (!is.list(adj)) {

        ## If a matrix, check to see if adjacency matrix
        if (is.matrix(adj)) {

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
            if (squaremat & binary & diag & symmetric) {

                ## Initialize object
                adjlist <- vector("list", nrow(adj))

                ## Loop through rows in matrix
                for (i in 1:nrow(adj)) {

                    ## Extract row
                    adjvec <- adj[, i]
                    ## Find elements it is adjacent to
                    inds <- which(adj == 1)
                    ## Remove self-adjacency
                    inds <- inds[inds != i, ]
                    ## Zero-index
                    inds <- inds - 1
                    ## Put in adjlist
                    adjlist[[i]] <- inds

                }

            } else { ## If not valid adjacency matrix, throw error
                stop("Please input valid adjacency matrix")
            }
        } else if (inherits(adj, "SpatialPolygonsDataFrame")) { ## shp object

            ## Convert shp object to adjacency list
            adjlist <- redist.adjacency(st_as_sf(adj))


        } else { ## If neither list, matrix, or shp, throw error
            stop("Please input an adjacency list, adjacency matrix, or Spatial
                 Polygons shp file")
        }

    } else if ("sf" %in% class(adj)) {
        adjlist <- redist.adjacency(adj)
    } else {

        ## Rename adjacency object as list
        adjlist <- adj

        ## Is list zero-indexed?
        minlist <- min(unlist(adjlist))
        maxlist <- max(unlist(adjlist))
        oneind <- (sum(minlist == 1, maxlist == length(adjlist)) == 2)
        zeroind <- (sum(minlist == 0, maxlist == (length(adjlist) - 1)) == 2)

        if (oneind) {
            ## Zero-index list
            for (i in 1:length(adjlist)) {
                adjlist[[i]] <- adjlist[[i]] - 1
            }
        } else if (!(oneind | zeroind)) {
            ## if neither oneind or zeroind, then stop
            stop("Adjacency list must be one-indexed or zero-indexed")
        }

    }

    if (is.null(init_plan) || isTRUE("smc" %in% init_plan)) {
        map <- redist_map(pop = total_pop,  pop_tol = ifelse(is.null(pop_tol), 0.05, pop_tol),
            ndists = ndists,  adj = adj)
        invisible(capture.output(
            init_plan <- redist_smc(map, nsims = 1, silent = TRUE),
            type = "message"))
        init_plan <- as.matrix(init_plan)[, 1]
    } else if (!is.null(init_plan) && "rsg" %in% init_plan) {
        ##############################################################################
        ## If no init_plan == rsg, use Random Seed and Grow                         ##
        ## (Chen and Rodden 2013) algorithm                                         ##
        ##############################################################################

        ## Set up target pop, strength of constraint (5%)
        if (is.null(pop_tol)) {
            pop_tol_rsg <- .05
        } else {
            pop_tol_rsg <- pop_tol
        }

        ## Print start
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

        if (verbose) {
            cat("\n", append = TRUE)
            cat(divider, append = TRUE)
            cat("Using redist.rsg() to generate starting values.\n\n", append = TRUE)
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
        if (!is.na(init_plan)[1]) {
            if (sum(is.na(init_plan)) > 0) {
                stop("You have NA's in your congressional districts. Please check the provided init_plan vector for NA entries.")
            }

            ndists <- length(unique(init_plan))
            divlist <- genAlConn(adjlist, init_plan)
            ncontig <- countpartitions(divlist)

            if (ncontig != ndists) {
                stop(paste("Your initial congressional districts have ", ndists,
                    " unique districts but ",
                    ncontig, " contigous connected components. Please provide a starting map with contigous districts.", sep = ""))
            }
        }
    }

    ###########################################################
    ## Check other inputs to make sure they are right length ##
    ###########################################################
    if ((length(total_pop) != length(adjlist)) | (sum(is.na(total_pop)) > 0)) {
        stop("Each entry in adjacency list must have a corresponding entry
              in vector of populations")
    }
    if ((length(init_plan) != length(adjlist)) | (sum(is.na(init_plan)) > 0)) {
        stop("Each entry in adjacency list must have an initial congressional
             district assignment")
    }

    ####################
    ## Zero-index cds ##
    ####################
    if (min(init_plan) != 0) {
        init_plan <- vctrs::vec_group_id(init_plan) - 1
    }
    if (length(unique(init_plan)) != (max(init_plan) + 1)) {
        stop("The district numbers in init_plan must be consecutive. The input to `init_plan` could not be transformed using `redist.sink.plan()`.")
    }

    ####################################################
    ## Calculate parity and population margin allowed ##
    ####################################################
    dists <- length(unique(init_plan))
    if (is.null(pop_tol)) {
        pop_tol <- 100
    }

    ########################
    ## Set up constraints ##
    ########################
    beta <- ifelse((temper %in% c(TRUE)), 0, 1)
    temperbeta <- ifelse(temper, "tempering", "none")

    ###################################
    ## Check if betaspacing provided ##
    ###################################
    if ("tempering" %in% temperbeta) {
        if (betaseq[1] == "powerlaw") {

            ## Generate power law sequence
            betaseq <- rep(NA, betaseqlength)
            for (i in 1:length(betaseq)) {
                betaseq[i] <- (0.1^((i - 1)/(length(betaseq) - 1)) - .1)/.9
            }

        } else if (is.vector(betaseq)) {
            betaseq <- betaseq
        } else if (!is.vector(betaseq) & betaseq[1] != "powerlaw") {
            stop("Please provide valid sequence of betas")
        }
        if (is.null(betaweights)) {
            betaweights <- rep(1, length(betaseq))
        }
    } else {
        betaseq <- c(1, 1, 1, 1)
        betaweights <- c(1, 1, 1, 1)
    }

    ## Reverse beta sequence
    betaseq <- rev(betaseq)

    ########################################
    ## Convert adjacent swaps flag to 0/1 ##
    ########################################
    adjswaps <- adjswaps*1

    #################
    ## Return list ##
    #################
    preprocout <- list(
        data = list(
            adjlist = adjlist,
            total_pop = total_pop,
            init_plan = init_plan
        ),
        params = list(
            pctdistparity = pop_tol,
            dists = dists,
            beta = beta,
            temperbeta = temperbeta,
            betaseq = betaseq,
            betaweights = betaweights,
            adjswaps = adjswaps
        )
    )

    class(preprocout) <- "redist"
    preprocout
}
