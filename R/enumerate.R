###########################################
## Author: Michael Higgins
## Institution: Princeton University
## Purpose: R wrapper to run full enumeration code
###########################################

redist.enumerate <- function(adjobj,
                             ndists = 2,
                             popvec = NULL,
                             nconstraintlow = NULL,
                             nconstrainthigh = NULL,
                             popcons = NULL,
                             contiguitymap = "rooks"){

    #########################
    ## Inputs to function: ##
    #########################
    ## adjobj - adjacency object of geographic units. Accepts adjlist, adjmat, or spatialPolygonsDataFrame
    ## ndists - number of partitions wanted
    ## popvec - vector of populations for geographic units
    ## nconstraintlow - minimum number of precincts contained in district
    ## nconstrainthigh - maximum number of precincts contained in district
    ## popcons - strenght of population constraint. Decimal
    ## contiguitymap - Use queens distance or rooks distance for creating an adjlist
    ##                 from map. Defaults to rooks.

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

