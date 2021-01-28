
#' Calculate Group Percent by District
#' 
#' \code{redist.group.percent} computes the percentage that a group makes up in
#' each district across a matrix of maps.
#'
#' @param district_membership A matrix with one row 
#' for each precinct and one column for each map. Required.
#' @param grouppop A numeric vector with the population of the group for every precinct.
#' @param fullpop A numeric vector with the population for every precinct.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' 
#' @importFrom foreach %do% %dopar% foreach
#' @return matrix with percent for each district
#' 
#' @export
#' 
#' @examples \dontrun{
#' data("fl25")
#' data(algdat)
#' cd <- algdat.p10$cdmat[,1:5]
#'
#' redist.group.percent(district_membership = cd, 
#'                     grouppop = fl25$BlackPop, 
#'                     fullpop = fl25$TotPop)
#' }
redist.group.percent <- function(district_membership, grouppop, fullpop, ncores = 1){
  
  if(!any(class(fullpop) %in% c('numeric', 'integer'))){
    stop('Please provide "fullpop" as a numeric vector.')
  }
  if(!any(class(grouppop) %in% c('numeric', 'integer'))){
    stop('Please provide "groupop" as a numeric vector.')
  }
  
  if(!any(class(district_membership) %in% c('integer', 'numeric', 'matrix'))){
    stop('Please provide "district_membership" as a matrix.')
  }
  
  if(!is.matrix(district_membership)){
    district_membership <- as.matrix(district_membership)
  }
  
  if(length(fullpop) != nrow(district_membership)){
    stop('Arguments "district_membership" and "fullpop" do not have same number of precincts.')
  }
  if(length(grouppop) != nrow(district_membership)){
    stop('Arguments "district_membership" and "groupop" do not have same number of precincts.')
  }
  
  apply(district_membership, 2, function(x){
    group_pop <- tapply(grouppop, x, sum)
    dist_pop <- tapply(fullpop, x, sum)
    return(group_pop/dist_pop)
  })
}
