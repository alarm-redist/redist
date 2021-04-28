###########################################
## Author: Michael Higgins
## Institution: Princeton University
## Purpose: R wrapper to run full enumeration code
###########################################

list_to_mat <- function(A){
    AM <- matrix(0,length(A),length(A));
    for(i in 1:length(A)){
        for(j in 1:length(A[[i]])){
            AM[i,A[[i]][j]+1] <- 1
        }
    }
    return(AM)
}

#' Deprecated: Exact Redistricting Plan Enumerator
#'
#' \code{redist.enumerate} uses a spanning-tree method to fully enumerate all
#' valid redistricting plans with $n$ districts given a set of geographic units.
#' \code{redist.enumerate} also allows users to implement minimum and maximum
#' numbers of geographic units per district, as well as population parity
#' requirements.
#'
#' @param adj An adjacency list, matrix, or object of class
#' \code{SpatialPolygonsDataFrame}.
#' @param ndists The desired number of congressional districts. The default is 2.
#' @param total_pop A vector of geographic unit populations. The default is
#' \code{NULL}.
#' @param nconstraintlow Lower bound for number of geographic units to include in
#' a district. The default is \code{NULL}.
#' @param nconstrainthigh Lower bound for number of geographic units to include
#' in a district. The default is \code{NULL}.
#' @param pop_tol The strength of the hard population constraint.
#' \code{pop_tol} = 0.05 means that any proposed swap that brings a district more
#' than 5\% away from population parity will be rejected. The default is
#' \code{NULL}.
#'
#' @param adjobj Deprecated, use adj. An adjacency list, matrix, or object of class
#' \code{SpatialPolygonsDataFrame}.
#' @param popvec Deprecated, use total_pop. A vector of geographic unit populations.
#' @param popcons Deprecated, use pop_tol. The strength of the hard population constraint.
#' \code{popcons} = 0.05 means that any proposed swap that brings a district more
#' than 5\% away from population parity will be rejected. The default is
#' \code{NULL}.
#' @param contiguitymap Use queens or rooks distance criteria for generating an
#' adjacency list from a "SpatialPolygonsDataFrame" data type. Default is "rooks".
#'
#' @details This function allows users to input a set of geographic units to
#' generate all valid partitions of $n$ congressional districts. The function
#' uses a set of spanning-tree methods to generate all valid, contiguous
#' partitions, which makes it more efficient than brute-force methods. However,
#' even with these methods, full redistricting problems quickly become
#' intractable, necessitating the use of the MCMC-based methods implemented in
#' \code{redist.mcmc}.
#'
#' @return \code{redist.enumerate} returns an object of class "list". Each entry
#' in the list is a vector of congressional district assignments, where the first
#' entry in the vector corresponds to the congressional district assignment of
#' the first geographic unit.
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' @concept enumerate
#' @export
redist.enumerate <- function(adj,
                             ndists = 2,
                             total_pop = NULL,
                             nconstraintlow = NULL,
                             nconstrainthigh = NULL,
                             pop_tol = NULL,
                             adjobj,
                             popvec,
                             popcons,
                             contiguitymap = "rooks"){
    .Deprecated(new = 'redist.enumpart')


    if(!missing(adjobj)){
        adj <- adjobj
        .Deprecated(new = 'adj', old = 'adjobj')
    }
    if(!missing(popvec)){
        total_pop <- popvec
        .Deprecated(new = 'total_pop', old = 'popvec')
    }
    if(!missing(popcons)){
        pop_tol <- popcons
        .Deprecated(old = 'popcons', new = 'pop_tol')
    }


    ## Warnings
    if(is.null(total_pop) & !is.null(pop_tol)){
        stop("If constraining on population, please provide a vector of populations for geographic units.")
    }
    if(!(contiguitymap %in% c("queens", "rooks"))){
        stop("Please supply `queens` or `rooks` for a distance criteria")
    }

    ############################################
    ## If not a list, convert adjlist to list ##
    ############################################
    ## NOTE - FOR ENUMERATION, WE WANT ONE-INDEXING VERSUS ZERO-INDEXING
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
                    ## Put in adjlist
                    adjlist[[i]] <- inds

                }

            }else { ## If not valid adjacency matrix, throw error
                stop("Please input valid adjacency matrix")
            }
        }else if(class(adj) == "SpatialPolygonsDataFrame"){ ## shp object

            ## Distance criterion
            queens <- ifelse(contiguitymap == "rooks", FALSE, TRUE)

            ## Convert shp object to adjacency list
            adjlist <- redist.adjacency(st_as_sf(adj))

            ## Change class to list
            class(adjlist) <- "list"

        }else{ ## If neither list, matrix, or shp, throw error
            stop("Please input an adjacency list, adjacency matrix, or Spatial
                 Polygons shp file")
        }

    }else{

        ## Rename adjacency object as list
        adjlist <- adj

        ## Is list zero-indexed?
        minlist <- min(unlist(adjlist))
        maxlist <- max(unlist(adjlist))
        oneind <- (sum(minlist == 1, maxlist == length(adjlist)) == 2)
        zeroind <- (sum(minlist == 0, maxlist == (length(adjlist) - 1)) == 2)

        if(zeroind){
            ## Zero-index list
            for(i in 1:length(adjlist)){
                adjlist[[i]] <- adjlist[[i]] + 1
            }
        }else if(!(oneind | zeroind)){
            ## if neither oneind or zeroind, then stop
            stop("Adjacency list must be one-indexed or zero-indexed")
        }

    }

    #########################
    ## Other preprocessing ##
    #########################

    ## Clean the input ##
    ## Store the number of nodes
    adjListLength <- length(adjlist)

    ## Set population constraint ##
    if(!is.null(pop_tol)){
        parity <- sum(total_pop) / ndists
        popConstraintLow <- parity - total_pop * parity
        popConstraintHigh <- parity + total_pop * parity
    }
    if(is.null(total_pop)){
        if(!is.null(total_pop)){
            popConstraintLow <- min(total_pop)
            popConstraintHigh <- sum(total_pop)
        }else{
            total_pop <- rep(1,adjListLength)
            popConstraintLow <- 1
            popConstraintHigh <- adjListLength
        }
    }

    ## If there is no pop vector,
    ## Default total_pop to vector of 1's,
    ## and default popConstraintLow and popConstraintHigh to 1 and adjListLength respectively
    if(is.null(total_pop)){
        total_pop <- rep(1,adjListLength)
        popConstraintLow <- 1
        popConstraintHigh <- adjListLength
    }

    ## De-facto minimums and maximum on number of units (if not specified)
    if(is.null(nconstraintlow)){
        nconstraintlow <- 1
    }

    ## The most amount of nodes to be contained within a block is
    ## numNodes - nconstraintlow*(ndists-1)
    if(is.null(nconstrainthigh)){
        nconstrainthigh <- adjListLength - nconstraintlow*(ndists-1)
    }

    ## Run the cpp function and return the output
    out <- cppGeneratePartitions(adjlist,
                                 ndists,
                                 total_pop,
                                 nconstraintlow,
                                 nconstrainthigh,
                                 popConstraintLow,
                                 popConstraintHigh)
    return(out)

}

#' DEFUNCT: Sample partitions using spanning trees
#'
#' \code{redist.samplepart} uses a spanning tree method to randomly sample
#' redistricting plans. DEFUNCT, use \code{\link{redist_smc}} instead.
#'
#'
#' @param adj An adjacency list, matrix, or object of class
#' \code{SpatialPolygonsDataFrame}.
#' @param ndists The desired number of congressional districts
#' @param total_pop Population vector for adjacency object. Provide if
#' filtering by population
#' @param pop_filter Boolean. Whether or not to filter on population parity.
#' Default is FALSE.
#' @param pop_tol Strength of population filter if filtering on
#' distance to parity.
#' @param pop_constraint Deprecated, use pop_tol. Strength of population filter if filtering on
#' distance to parity.
#' @param nsims Number of samples to draw. Default is 1000.
#' @param ncores Number of cores to parallelize over for parity calculation and
#' compactness calculation. Default is 1.
#'
#' @param adjobj Deprecated, use adj. An adjacency list, matrix, or object of class
#' \code{SpatialPolygonsDataFrame}.
#' @param popvec Deprecated, use total_pop. Population vector for adjacency object. Provide if
#' filtering by population
#' @param contiguitymap Use queens or rooks distance criteria for generating an
#' adjacency list from a "SpatialPolygonsDataFrame" data type.
#' Default is "rooks".
#' @param nsamp Deprecated, use nsims. Number of samples to draw. Default is 1000.
#' @param n_cores Deprecated, use ncores. Number of cores to parallelize over for parity calculation and
#' compactness calculation. Default is 1.
#'
#' @return \code{redist.samplepart} returns a list where the first entry is the
#' randomly sampled redistricting plan, and the second entry is the number of
#' possible redistricting plans from the implied spanning tree.
#' @export
#' @importFrom parallel mclapply
redist.samplepart <- function(adj,  ndists, total_pop = NULL,
                              pop_filter = FALSE, pop_tol = 0.5,
                               nsims = 1000, ncores = 1,
                              adjobj, popvec, pop_constraint,
                              contiguitymap = "rooks", nsamp, n_cores){
    .Defunct("redist_smc", "redist",
             "To sample plans using spanning trees, use `redist_smc` or `redist_mergesplit`.")
}
