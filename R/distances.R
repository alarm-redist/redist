#' Compute Distance between Partitions
#'
#' @param district_membership A matrix with one row for each precinct and one 
#' column for each map. Required.
#' @param measure String vector indicating which distances to compute. Implemented 
#' currently are "Hamming", "Manhattan", and "Euclidean". Use all to return all implemented 
#' measures.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' 
#' @return list of matrices of distances with one for each distance measure selected
#'
#' @examples \dontrun{
#' data("algdat.p10")
#' distances <- redist.distances(district_membership = algdat.p10$cdmat)
#' distances$Hamming[1:5,1:5]
#' }
#' @export
redist.distances <- function(district_membership, measure = "Hamming", ncores = 1){
  
  #check inputs
  if(measure == "all"){
    measure <- c("Hamming", "Manhattan", "Euclidean")
  }
  
  # init vars
  distances <- list()
  name <- c()
  done <- 0
  
  # Compute Hamming Distance Metric
  if("Hamming" %in% measure){
  nc <- min(ncores, ncol(district_membership))
  if (nc == 1){
    `%oper%` <- `%do%`
  } else {
    `%oper%` <- `%dopar%`
    cl <- makeCluster(nc, , setup_strategy = 'sequential')
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  ham <- foreach(map = 1:ncol(district_membership), .combine = "cbind") %oper% {
    hamming(v = district_membership[,map], m = district_membership)
  }
  colnames(ham) <- NULL
  
  done = done + 1
  distances[[done]] <- ham
  names(distances)[done] <- "Hamming"
  }
  
  
  # Compute Manhattan Distance Metric
  if("Manhattan" %in% measure){
    nc <- min(ncores, ncol(district_membership))
    if (nc == 1){
      `%oper%` <- `%do%`
    } else {
      `%oper%` <- `%dopar%`
      cl <- makeCluster(nc, , setup_strategy = 'sequential')
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
    }
    
    man <- foreach(map = 1:ncol(district_membership), .combine = "cbind") %oper% {
      minkowski(v = district_membership[,map], m = district_membership, p = 1)
    }
    colnames(man) <- NULL
    
    done = done + 1
    distances[[done]] <- man
    names(distances)[done] <- "Manhattan"
  }
  
  # Compute Euclidean Distance Metric
  if("Euclidean" %in% measure){
    nc <- min(ncores, ncol(district_membership))
    if (nc == 1){
      `%oper%` <- `%do%`
    } else {
      `%oper%` <- `%dopar%`
      cl <- makeCluster(nc, , setup_strategy = 'sequential')
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
    }
    
    euc <- foreach(map = 1:ncol(district_membership), .combine = "cbind") %oper% {
      minkowski(v = district_membership[,map], m = district_membership, p = 2)
    }
    colnames(euc) <- NULL
    
    done = done + 1
    distances[[done]] <- euc
    names(distances)[done] <- "Euclidean"
  }
  
  return(distances)
}

utils::globalVariables(names = "map")