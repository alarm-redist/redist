########################################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/07/15
## Date Modified: 2020/07/16
## Purpose: R function to compute gerrymandering metrics
########################################################

#' Calculate gerrymandering metrics for a set of districts
#' 
#' \code{redist.metrics} is used to compute different gerrymandering metrics for a
#' set of maps.
#' 
#' @param district_membership A numeric vector (if only one map) or matrix with one row 
#' for each precinct and one column for each map. Required.
#' @param measure A vector with a string for each measure desired.
#' @param rvote A numeric vector with the Republican vote for each precinct.
#' @param dvote A numeric vector with the Democratic vote for each precinct.
#' @param nloop A numeric to specify loop number. Defaults to 1 if only one map provided 
#' and the column number if multiple maps given.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' 
#' @details This function computes specified compactness scores for a map.  If 
#' there is more than one precinct specified for a map, it aggregates to the district level
#' and computes one score.
#' 
#' Seats is computed as the expected number of seats with no change in votes.
#' 
#' 
#' 
#' @export
redist.metrics <- function(district_membership, 
                           measure = "Seats", 
                           rvote, dvote, 
                           nloop = 1, 
                           ncores = 1){
  
  # Check Inputs
  if(class(district_membership) == 'redist'){
    district_membership <- district_membership$partitions
  }  
  
  if(!any(class(district_membership) %in% c('numeric', 'integer', 'matrix'))){
    stop('Please provide "district_membership" as a numeric vector or matrix.')
  }
  if(!is.matrix(district_membership)){
    district_membership <- as.matrix(district_membership)
  } 
  
  if(measure == "all"){
    measure <-  c("Seats")
  }
  match.arg(arg = measure,several.ok = TRUE, choices = c("Seats"))
  
 # TODO add checks for rvote, dvote class + that all dimensions conform for n prec
  
  
  if(class(nloop) != 'numeric'){
    stop('Please provide "nloop" as a numeric.')
  }
  
  if(class(ncores) != 'numeric'){
    stop('Please provide "ncores" as a numeric.')
  }
  
  # Precompute a few useful variables
  nd <- length(unique(district_membership[,1]))
  np <- nrow(district_membership)
  
  
  
  
}