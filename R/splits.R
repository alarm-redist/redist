#' Count County Splits
#'
#' @param district_membership A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param counties A vector of county names or county ids.
#'
#' @return integer vector with one number for each map
#' @export
#'
redist.splits <- function(district_membership, counties){

  if(missing(district_membership)){
    stop('Please provide an argument to district_membership.')
  }
  if(missing(counties)){
    stop('Please provide an argument to counties.')
  }
  if(class(counties) %in% c('character', 'numeric','integer')){
    uc <- unique(sort(counties))
    county_id <- rep(0, nrow(district_membership))
    for(i in 1:nrow(district_membership)){
      county_id[i] <- which(uc == counties[i])
    }
  } else{
    stop('Please provide "counties" as a character, numeric, or integer vector.')
  }

  county_splits <- splits(district_membership, community = county_id)

  return(county_splits)

}

#' Identify which counties re split by a plan
#'
#' @param plan A vector of precinct/unit assignments
#' @param counties A vector of county names or county ids.
#'
#' @return A logical vector which is \code{TRUE} for precincts belonging to
#'   counties which are split
is_county_split = function(plan, counties) {
    (tapply(plan, counties, FUN=function(y) length(unique(y))) > 1)[counties]
}
