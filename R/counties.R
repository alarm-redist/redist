#' Relabel Discontinuous Counties
#'
#' @param adj adjacency list
#' @param counties character vector of county names
#' @param simplify boolean - TRUE returns a numeric vector of ids, while FALSE
#' appends a number when there are multiple connected components.
#'
#' @return character vector of county names
#'
#' @importFrom dplyr group_by mutate ungroup rowwise
#'
#' @concept prepare
#' @export
#' @examples
#' set.seed(2)
#' data(fl25)
#' data(fl25_adj)
#' counties <- sample(c(rep('a', 20), rep('b', 5)))
#' redist.county.relabel(fl25_adj, counties)
#'
redist.county.relabel <- function(adj, counties,  simplify = TRUE){

  if(missing(adj)){
    stop('Please provide an argument to adj.')
  }

  if(missing(counties)){
    stop('Please provide a character vector of counties.')
  }
  if(length(adj) != length(counties)){
    stop('adj and group are different lengths.')
  }

  if('numeric' %in% class(counties) | 'integer' %in% class(counties) ){
    counties <- as.character(counties)
  }

  groups <- rep(0, length(counties))
  sorted <- sort(unique(counties))
  for(i in 1:length(counties)){
    groups[i] <- which(sorted == counties[i])
  }

  component <- data.frame(counties = counties, comp = contiguity(adj, groups)) %>%
    group_by(counties) %>%
    mutate(comps = max(comp)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(countiescomp = ifelse(comps > 1, paste0(counties, '-', comp), counties)) %>%
    ungroup()
  if(simplify){
    return(redist.county.id(component$countiescomp))
  } else{
    return(component$countiescomp)
  }
}


#' Create County IDs
#'
#' @param counties vector of counties, required.
#'
#' @return A vector with an ID that corresponds from 1:n counties
#'
#' @concept prepare
#' @export
#' @examples
#' set.seed(2)
#' counties <- sample(c(rep('a', 20), rep('b', 5)))
#' redist.county.id(counties)
#'
redist.county.id <- function(counties){
  if(class(counties) %in% c('character', 'numeric','integer')){
    uc <- unique(sort(counties))
    county_id <- rep(0, length(counties))
    for(i in 1:length(counties)){
      county_id[i] <- which(uc == counties[i])
    }
  } else{
    stop('Please provide "counties" as a character, numeric, or integer vector.')
  }

  return(county_id)
}

globalVariables(c('comp', 'comps'))
