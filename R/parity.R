#' Calculates Population Parity
#' 
#' \code{redist.parity} computes the population parity of a matrix of maps.
#'
#' @param district_membership A matrix with one row 
#' for each precinct and one column for each map. Required.
#' @param population A numeric vector with the population for every precinct.
#' @return numeric vector with the population parity for each column
#' @export
#'
#' @examples
redist.parity <- function(district_membership, population){
  
  if(!any(class(population %in% c('numeric', 'integer')))){
    stop('Please provide "population" as a numeric vector.')
  }
  
  if(!any(class(district_membership) %in% c('numeric', 'matrix'))){
    stop('Please provide "district_membership" as a matrix.')
  }
  
  if(!is.matrix(district_membership)){
    district_membership <- as.matrix(district_membership)
  }
  
  if(length(population) != nrow(district_membership)){
    stop('Arguments "district_membership" and "population" do not have same number of precincts.')
  }
  
  out <- apply(district_membership, 2, function(x){
    distpop <- tapply(population, x, sum)
    parpop <- sum(distpop) / length(distpop)
    max(abs(distpop / parpop - 1))
  })
  return(out)
}