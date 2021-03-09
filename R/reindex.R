#' Reorders district numbers
#' 
#' Ensures that for each column in the plans object, the first 
#' district listed is 1, the second is 2, up to n districts. Assumes that all
#' columns have the same number of districts as the first.
#' @param plans A numeric vector (if only one map) or 
#' matrix with one row for each precinct and one column for each map.
#' @param district_membership Deprecated, use plans. A numeric vector (if only one map) or 
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
redist.reindex <- function(plans, district_membership){
  if(!missing(district_membership)){
    .Deprecated(new = 'plan', old = 'district_membership')
    plan <- district_membership
  }
  
  # Check inputs
  if(missing(plans)){
    stop('"plans" is required.')
  }
  if(any(class(plans)%in% c('numeric', 'integer'))){
    plans <- as.matrix(plans)
  }
  if(!('matrix' %in% class(plans))){
    stop('Please provide "plans" as a matrix.')
  }
  
  # Prep objects for Rcpp
  nd <- length(unique(plans[,1]))
  
  # reindex!
  return(reindex(dm = plans, nd = nd))
}
