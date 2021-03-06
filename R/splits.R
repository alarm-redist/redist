#' Count County Splits
#' 
#' @param plans A numeric vector (if only one map) or matrix with one row 
#' for each precinct and one column for each map. Required.
#' @param district_membership Deprecated, use plans. A numeric vector (if only one map) or matrix with one row 
#' for each precinct and one column for each map. Required.
#' @param counties A vector of county names or county ids.
#'
#' @return integer vector with one number for each map
#' @export
#' 
redist.splits <- function(plans, district_membership, counties){
   
  if(!missing(district_membership)){
    .Deprecated(new = 'plans', old = 'district_membership')
    plans <- district_membership
  }
  
  if(missing(plans)){
    stop('Please provide an argument to plans.')
  } 
  if(missing(counties)){
    stop('Please provide an argument to counties.')
  }
  if(class(counties) %in% c('character', 'numeric','integer')){
    uc <- unique(sort(counties))
    county_id <- rep(0, nrow(plans))
    for(i in 1:nrow(plans)){
      county_id[i] <- which(uc == counties[i])
    }
  } else{
    stop('Please provide "counties" as a character, numeric, or integer vector.')
  }
  
  county_splits <- splits(plans, community = county_id)
  
  return(county_splits)

}
