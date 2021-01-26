#' Relabel Discontinuous Counties
#'
#' @param adjacency adjacency list
#' @param counties character vector of county names
#'
#' @return character vector of county names
#' 
#' @importFrom dplyr group_by mutate ungroup rowwise
#' 
#' @export
#' @examples \dontrun{
#' 
#' set.seed(2)
#' counties <- sample(c(rep('a', 20), rep('b', 5)))
#' redist.relabel.counties(counties)
#' 
#' }
redist.county.relabel <- function(adjacency, counties){
  if(missing(adjacency)){
    stop('Please provide an argument to adjacency.')
  }
  
  if(missing(counties)){
    stop('Please provide a character vector of counties.')
  }
  if(length(adjacency) != length(counties)){
    stop('Adjacency and group are different lengths.')
  }
  
  groups <- rep(0, length(counties))
  sorted <- sort(unique(counties))
  for(i in 1:length(counties)){
    groups[i] <- which(sorted == counties[i])
  }

  component <- data.frame(counties = counties, comp = contiguity(adjacency, groups)) %>% 
    group_by(counties) %>% 
    mutate(comps = max(comp)) %>% 
    ungroup() %>% 
    rowwise() %>% 
    mutate(countiescomp = ifelse(comps > 1, paste0(counties, comp), counties)) %>% 
    ungroup()
  
  return(component$countiescomp) 
}


#' Create County IDs
#'
#' @param counties vector of counties, required.
#'
#' @return A vector with an ID that corresponds from 1:n counties
#' @export
#' @examples \dontrun{
#' set.seed(2)
#' counties <- sample(c(rep('a', 20), rep('b', 5)))
#' redist.county.id(counties)
#' }
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