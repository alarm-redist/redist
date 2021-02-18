#' Reorders district numbers
#' 
#' Ensures that for each column in the district_membership object, the first 
#' district listed is 1, the second is 2, up to n districts. Assumes that all
#' columns have the same number of districts as the first.
#'
#' @param district_membership A numeric vector (if only one map) or 
#' matrix with one row for each precinct and one column for each map.
#'
#' @return integer matrix 
#' @export
#'
#' @examples \dontrun{
#' cds <- matrix(c(rep(c(4L,5L,2L,1L,3L),5), 
#' rep(c(5L,4L,3L,2L,1L),2), rep(c(4L,5L,2L,1L,3L),3)), nrow = 25)
#' redist.reorder(cds)
#' }
redist.reindex <- function(district_membership){
  # Check inputs
  if(missing(district_membership)){
    stop('"district_membership" is required.')
  }
  if(any(class(district_membership)%in% c('numeric', 'integer'))){
    district_membership <- as.matrix(district_membership)
  }
  if(!('matrix' %in% class(district_membership))){
    stop('Please provide "district_membership" as a matrix.')
  }
  
  # Prep objects for Rcpp
  nd <- length(unique(district_membership[,1]))
  
  # reindex!
  return(reindex(dm = district_membership, nd = nd))
}
