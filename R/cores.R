#' Identify Cores of a District (Heuristic)
#'
#' Creates a grouping ID to unite geographies and perform analysis on a smaller set of
#' precincts. Given a k, it identifies all precincts within k edges of another district
#' on a graph. Each precinct within k of another district gets it own group. Connected
#' components more than k away are givent the same grouping variable.
#'
#' @param adjacency A zero indexed adjacency list.
#' @param district_membership An integer vector or matrix column of district assignments.
#' @param k Number of steps to check for. Defaults to 1.
#' @param focus Optional. Integer. A single district to focus on.
#' @param simplify Optional. Logical. Whether to return extra information or just grouping ID. 
#'
#' @return integer vector if simplify is false. Otherwise returns a tibble with the grouping 
#' variable and additional information.
#' 
#' @importFrom dplyr row_number cur_group_id
#' 
#' @export
redist.identify.cores <- function(adjacency, district_membership, k = 1, focus = NULL, simplify = TRUE){
  if(missing(adjacency)){
    stop('Please provide an object to adjacency.')
  }
  
  if(! 'list' %in% class(adjacency)){
    stop('adjacency must be an adjacency list.')
  }
  
  if(missing(district_membership)){
    stop('Please provide an object to district_membership.')
  }
  if('matrix' %in% class(district_membership)){
    district_membership <- district_membership[,1]
  }
  

  # init a nice empty list
  cd_within_k <- lapply(1:length(district_membership), FUN = function(x){integer(0)})

  core <- cores(adj = adjacency, dm = district_membership, k = k, cd_within_k = cd_within_k)

  if(!is.null(focus)){
    idx <- unlist(lapply(core$cd_within_k, FUN = function(x){focus %in% x})) | district_membership == focus
    
    core$k <- ifelse(idx, core$k, 0)
    
    conncomp <- update_conncomp(dm  = core$dm, kvec = core$k, adj = adjacency)
    core$conncomp <- conncomp
  } 
  
  tb <- tibble(dm = district_membership, k = core$k, cc = core$conncomp) %>% 
    group_by(dm, k) %>% 
    mutate(gid = row_number()) %>% 
    ungroup() %>% 
    mutate(gid = ifelse(k == 0, cc, gid)) %>% 
    mutate(gid = paste0(dm, '-', k, '-', gid)) %>% group_by(gid) %>% 
    mutate(group = cur_group_id()) %>% ungroup()
  
  gid <- tb$group

  
  if(simplify){
    return(gid)
  } else{
    core$gid <- gid
    return(core)
  }

  
}

globalVariables(c('dm', 'cc'))
