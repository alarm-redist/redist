###########################################
## Author: Michael Higgins
## Institution: Princeton University
## Purpose: R wrapper to run full enumeration code
###########################################

#' Exact Redistricting Plan Enumerator
#'
#' \code{redist.enumerate} uses a spanning-tree method to fully enumerate all
#' valid redistricting plans with $n$ districts given a set of geographic units.
#' \code{redist.enumerate} also allows suers to implement minimum and maximum
#' numbers of geographic units per district, as well as population parity
#' requirements.
#'
#' @usage redist.enumerate(adjobj,
#' ndists = 2, popvec = NULL, nconstraintlow = NULL,
#' nconstrainthigh = NULL, popcons = NULL, contiguitymap = "rooks")
#'
#' @param adjobj An adjacency list, matrix, or object of class
#' \code{SpatialPolygonsDataFrame}.
#' @param ndists The desired number of congressional districts. The default is 2.
#' @param popvec A vector of geographic unit populations. The default is
#' \code{NULL}.
#' @param nconstraintlow Lower bound for number of geographic units to include in
#' a district. The default is \code{NULL}.
#' @param nconstrainthigh Lower bound for number of geographic units to include
#' in a district. The default is \code{NULL}.
#' @param popcons The strength of the hard population constraint.
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
#' @examples
#' \dontrun{
#' data(algdat.pfull)
#' test <- redist.enumerate(adjobj = algdat.pfull$adjlist)
#' }
#' @export
redist.enumerate <- function(adjobj,
                             ndists = 2,
                             popvec = NULL,
                             nconstraintlow = NULL,
                             nconstrainthigh = NULL,
                             popcons = NULL,
                             contiguitymap = "rooks"){

    ## Warnings
    if(is.null(popvec) & !is.null(popcons)){
        stop("If constraining on population, please provide a vector of populations for geographic units.")
    }
    if(!(contiguitymap %in% c("queens", "rooks"))){
        stop("Please supply `queens` or `rooks` for a distance criteria")
    }

    ############################################
    ## If not a list, convert adjlist to list ##
    ############################################
    ## NOTE - FOR ENUMERATION, WE WANT ONE-INDEXING VERSUS ZERO-INDEXING
    if(!is.list(adjobj)){

        ## If a matrix, check to see if adjacency matrix
        if(is.matrix(adjobj)){

            ## Is it square?
            squaremat <- (nrow(adjobj) == ncol(adjobj))
            ## All binary entries?
            binary <- ((length(unique(c(adjobj))) == 2) &
                           (sum(unique(c(adjobj)) %in% c(0, 1)) == 2))
            ## Diagonal elements all 1?
            diag <- (sum(diag(adjobj)) == nrow(adjobj))
            ## Symmetric?
            symmetric <- isSymmetric(adjobj)

            ## If all are true, change to adjlist and automatically zero-index
            if(squaremat & binary & diag & symmetric){
                
                ## Initialize object
                adjlist <- vector("list", nrow(adjobj))

                ## Loop through rows in matrix
                for(i in 1:nrow(adjobj)){

                    ## Extract row
                    adjvec <- adjobj[,i]
                    ## Find elements it is adjacent to
                    inds <- which(adjobj == 1)
                    ## Remove self-adjacency
                    inds <- inds[inds != i,]
                    ## Put in adjlist
                    adjlist[[i]] <- inds
                    
                }
                
            }else { ## If not valid adjacency matrix, throw error
                stop("Please input valid adjacency matrix")
            }
        }else if(class(adjobj) == "SpatialPolygonsDataFrame"){ ## shp object

            ## Distance criterion
            queens <- ifelse(contiguitymap == "rooks", FALSE, TRUE)
            
            ## Convert shp object to adjacency list
            adjlist <- poly2nb(adjobj, queen = queens)
            
            ## Change class to list
            class(adjlist) <- "list"
            
        }else{ ## If neither list, matrix, or shp, throw error
            stop("Please input an adjacency list, adjacency matrix, or Spatial
                 Polygons shp file")
        }

    }else{

        ## Rename adjacency object as list
        adjlist <- adjobj

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
    if(!is.null(popcons)){
        parity <- sum(popvec) / ndists
        popConstraintLow <- parity - popcons * parity
        popConstraintHigh <- parity + popcons * parity
    }
    if(is.null(popcons)){
        if(!is.null(popvec)){
            popConstraintLow <- min(popvec)
            popConstraintHigh <- sum(popvec)
        }else{
            popVec <- rep(1,adjListLength)
            popConstraintLow <- 1
            popConstraintHigh <- adjListLength
        }
    }
    
    ## If there is no pop vector, 
    ## Default popvec to vector of 1's,
    ## and default popConstraintLow and popConstraintHigh to 1 and adjListLength respectively
    if(is.null(popvec)){
        popvec <- rep(1,adjListLength)
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
                                 popvec,
                                 nconstraintlow,
                                 nconstrainthigh,
                                 popConstraintLow,
                                 popConstraintHigh)
    return(out)
	
}

