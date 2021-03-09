#' Identify Cores of a District (Heuristic)
#'
#' Creates a grouping ID to unite geographies and perform analysis on a smaller set of
#' precincts. Given a k, it identifies all precincts within k edges of another district
#' on a graph. Each precinct within k of another district gets it own group. Connected
#' components more than k away are givent the same grouping variable.
#'
#' @param adj zero indexed adjacency list.
#' @param plan An integer vector or matrix column of district assignments.
#' @param k Number of steps to check for. Defaults to 1.
#' @param focus Optional. Integer. A single district to focus on.
#' @param simplify Optional. Logical. Whether to return extra information or just grouping ID. 
#' @param adjacency Deprecated, use adj. A zero indexed adjacency list.
#' @param district_membership Deprecated, use plan. An integer vector or matrix column of district assignments.
#'
#' @return integer vector if simplify is false. Otherwise returns a tibble with the grouping 
#' variable and additional information.
#' 
#' @importFrom dplyr row_number cur_group_id
#' 
#' @export
redist.identify.cores <- function(adj, adjacency, plan, district_membership, k = 1, focus = NULL, simplify = TRUE){
  if(!missing(adjacency)){
    adj <- adjacency
    .Deprecated(new = 'adj', old = 'adjacency')
  }
  if(!missing(district_membership)){
    plan <- district_membership
    .Deprecated(new = 'plan', old = 'district_membership')
  }
  
  
  if(missing(adj)){
    stop('Please provide an object to adj.')
  }
  
  if(! 'list' %in% class(adj)){
    stop('adj must be an adj list.')
  }
  
  if(missing(plan)){
    stop('Please provide an object to plan.')
  }
  if('matrix' %in% class(plan)){
    plan <- plan[,1]
  }
  

  # init a nice empty list
  cd_within_k <- lapply(1:length(plan), FUN = function(x){integer(0)})

  core <- cores(adj = adj, dm = plan, k = k, cd_within_k = cd_within_k)

  if(!is.null(focus)){
    idx <- unlist(lapply(core$cd_within_k, FUN = function(x){focus %in% x})) | plan == focus
    
    core$k <- ifelse(idx, core$k, 0)
    
    conncomp <- update_conncomp(dm  = core$dm, kvec = core$k, adj = adj)
    core$conncomp <- conncomp
  } 
  
  tb <- tibble(dm = plan, k = core$k, cc = core$conncomp) %>% 
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

#' Uncoarsen a District Matrix
#'
#' After a cores analysis or other form of coarsening, sometimes you need
#' to be at the original geography level to be comparable. This takes in a
#' coarsened matrix and uncoarsens it to the original level
#' 
#' @param plans A coarsened matrix of plans.
#' @param district_membership Deprecated, use plans. A coarsened district membership matrix
#' @param group_index The index used to coarsen the shape.
#'
#' @return matrix
#' @export
redist.uncoarsen <- function(plans, district_membership, group_index){
  if(!missing(district_membership)){
    plan <- district_membership
    .Deprecated(new = 'plan', old = 'district_membership')
  }
  
  
  uncoarse <- matrix(nrow = length(group_index), 
                     ncol = ncol(plans))
  
  remain <- sort(unique(group_index))
  
  
  for(i in 1:length(group_index)){
    uncoarse[i,] <- plans[which(group_index[i] == remain),]
  }
  
  return(uncoarse)
}


globalVariables(c('dm', 'cc'))
