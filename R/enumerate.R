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

#' Sample partitions using spanning trees
#'
#' \code{redist.samplepart} uses a spanning tree method to randomly sample
#' redistricting plans.
#'
#' @usage redist.samplepart(adjobj, ndists, popvec, pop_filter, pop_constraint,
#' contiguitymap, nsamp, n_cores)
#'
#' @param adjobj An adjacency list, matrix, or object of class
#' \code{SpatialPolygonsDataFrame}.
#' @param ndists The desired number of congressional districts
#' @param popvec Population vector for adjacency object. Provide if
#' filtering by population
#' @param pop_filter Boolean. Whether or not to filter on population parity.
#' Default is FALSE.
#' @param pop_constraint Strength of population filter if filtering on
#' distance to parity.
#' @param contiguitymap Use queens or rooks distance criteria for generating an
#' adjacency list from a "SpatialPolygonsDataFrame" data type.
#' Default is "rooks".
#' @param nsamp Number of samples to draw. Default is 1000.
#' @param n_cores Number of cores to parallelize over for parity calculation and
#' compactness calculation. Default is 1.
#'
#' @return \code{redist.samplepart} returns a list where the first entry is the
#' randomly sampled redistricting plan, and the second entry is the number of
#' possible redistricting plans from the implied spanning tree.
#' @export
#' @importFrom parallel mclapply
redist.samplepart <- function(adjobj, ndists, popvec = NULL,
                              pop_filter = FALSE, pop_constraint = .5,
                              contiguitymap = "rooks", nsamp = 1000,
                              n_cores = 1){

    ## if(compact_filter & !inherits(adjobj, "SpatialPolygonsDataFrame")){
    ##     stop("If filtering on compactness, adjobj must be of class SpatialPolygonsDataFrame.")
    ## }
    if(pop_filter & is.null(popvec)){
        stop("If filtering on population, you must provide a vector of populations
for each geographic unit.")
    }
    if(!is.null(popvec)){
        if(pop_filter & sum(is.na(popvec)) > 0){
            stop("You have NAs in your vector of geographic unit populations.")
        }
    }
    
    ## --------------------------------------
    ## If not a list, convert adjlist to list
    ## --------------------------------------
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
                    ## Zero-index
                    inds <- inds - 1
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
            
            ## Zero-index list
            for(i in 1:length(adjlist)){
                adjlist[[i]] <- adjlist[[i]] - 1
            }
            
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
        
        if(oneind){
            ## Zero-index list
            for(i in 1:length(adjlist)){
                adjlist[[i]] <- adjlist[[i]] - 1
            }
        }else if(!(oneind | zeroind)){
            ## if neither oneind or zeroind, then stop
            stop("Adjacency list must be one-indexed or zero-indexed")
        }
        
    }
    if(pop_filter & length(popvec) != length(adjlist)){
        stop("Your population vector does not contain an entry for every unit.")
    }

    ## -------------------------
    ## Run enumeration algorithm
    ## -------------------------
    ## Run spanning tree
    cat("Sampling partitions using spanning tree method.\n")
    enum_out <- sample_partition(
        aList = adjlist, aMat = list_to_mat(adjlist),
        num_partitions = ndists, num_samples = nsamp,
        threads = n_cores
    )

    ## --------------
    ## Filtering step
    ## --------------
    if(pop_filter){
        cat("Calculating distance from population parity for sampled partitions.\n")
        popdist <- unlist(
            mclapply(1:ncol(enum_out$partitions), function(x, pops = popvec){
                part <- enum_out$partitions[,x]
                targ <- sum(pops) / length(unique(part))
                tab <- tapply(pops, part, sum)
                return(max(abs((tab/targ - 1))))
            }, mc.cores = n_cores)
        )
        inds_pop <- which(popdist <= pop_constraint)
        ## if(!compact_filter){
            inds_sub <- inds_pop
        ## }
    }
    ## if(compact_filter){
    ##     cat("Calculating Polsby-Popper score for sampled partitions.\n")
    ##     shp_sf <- st_as_sf(adjobj)
    ##     cpct <- unlist(
    ##         mclapply(1:ncol(enum_out$partitions), function(x, sf_obj = shp_sf){
    ##             part <- enum_out$partitions[,x]
    ##             dists <- unique(part)
    ##             pp <- rep(NA, length(dists))
    ##             for(i in 1:length(dists)){
    ##                 sf_sub <- st_combine(sf_obj[part == dists[i],])
    ##                 ar <- st_area(sf_sub)
    ##                 per <- st_length(st_boundary(st_union(sf_sub)))
    ##                 pp[i] <- 4 * pi * ar / per^2
    ##             }
    ##             return(min(pp))
    ##         }, mc.cores = n_cores)
    ##     )
    ##     inds_compact <- which(cpct >= polsbypopper_constraint)
    ##     if(!pop_filter){
    ##         inds_sub <- inds_compact
    ##     }
    ## }
    ## if(pop_filter & compact_filter){
    ##     inds_sub <- intersect(inds_pop, inds_compact)
    ## }

    ## Subset
    ## if(compact_filter | pop_filter){
    if(pop_filter){
        if(length(inds_sub) == 0){
            cat("No draws found within target parameters. Returning all partitions.\n")
            inds_sub <- 1:ncol(enum_out$partitions)
        }
        enum_out$partitions <- enum_out$partitions[,inds_sub]
        enum_out$prob_partitions <- enum_out$prob_partitions[inds_sub]
        if(pop_filter){
            enum_out$distance_parity <- popdist[inds_sub]
        }
        ## if(compact_filter){
        ##     enum_out$polsby_popper <- cpct[inds_sub]
        ## }
    }
    
    return(enum_out)
    
}
