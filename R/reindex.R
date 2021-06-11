#' Reorders district numbers
#'
#' Ensures that for each column in the plans object, the first
#' district listed is 1, the second is 2, up to n districts. Assumes that all
#' columns have the same number of districts as the first.
#' @param plans A numeric vector (if only one map) or
#' matrix with one row for each precinct and one column for each map.
#'
#' @return integer matrix
#' @export
#'
#' @examples
#' cds <- matrix(c(rep(c(4L,5L,2L,1L,3L),5),
#' rep(c(5L,4L,3L,2L,1L),2), rep(c(4L,5L,2L,1L,3L),3)), nrow = 25)
#' redist.reorder(cds)
#' 
redist.reorder <- function(plans){
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


#' Sink Plans to 1:ndists
#'
#' Takes a plan and renumbers it to be from 1:ndists
#'
#' @param plan vector of assignments, required.
#'
#' @return A vector with an ID that corresponds from 1:ndists
#'
#' @concept prepare
#' @export
#' @examples
#' data(fl25_enum)
#' plan <- fl25_enum$plans[, 5118]
#' # Subset based on something:
#' plan <- plan[plan!=2]
#' plan <- redist.sink.plan(plan)
#' # Now plan can be used with redist.flip()
#' plan
#'
#'
redist.sink.plan <- function(plan){
  if(class(plan) %in% c('character', 'numeric','integer')){
    uc <- unique(sort(plan))
    plan_id <- rep(0, length(plan))
    for(i in 1:length(plan)){
      plan_id[i] <- which(uc == plan[i])
    }
  } else{
    stop('Please provide "plan" as a  numeric or integer vector.')
  }

  return(plan_id)
}
