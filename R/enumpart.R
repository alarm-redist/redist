#' Initialize enumpart
#'
#'This ensures that the enumerate partitions programs is prepared to run.
#'This must be run once per install of the redist package.
#'
#' @return 0 on success
#' @export
#' @references
#' Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
#' "The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
#' Forthcoming, Statistics and Public Policy.
#'
#' @concept enumerate
#' @examples \dontrun{
#' redist.init.enumpart()
#' }
redist.init.enumpart <- function(){
  # Update makefile to direct to library only if Windows
  if(Sys.info()[['sysname']] == 'Windows'){
    makecontent <- readLines(system.file('enumpart/Makefile', package = 'redist'))
    makecontent[7] <- "\tg++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG -lpsapi"
    writeLines(text = makecontent, con = system.file('enumpart/Makefile', package = 'redist'))
  }

  servr::make(dir = system.file('enumpart', package = 'redist'), verbose = FALSE)

  if(Sys.info()[['sysname']] == 'Windows'){
    sys::exec_wait('python', args= c('-m', 'pip', 'install', 'networkx', '--user'))
  } else {
    sys::exec_wait('python3', args= c('-m', 'pip', 'install', 'networkx', '--user'))
  }


  # Necessary to avoid bad CRAN submissions:
  if(Sys.info()[['sysname']] == 'Windows'){
    makecontent <- readLines(system.file('enumpart/Makefile', package = 'redist'))
    makecontent[7] <- "\tg++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG"
    writeLines(text = makecontent, con = system.file('enumpart/Makefile', package = 'redist'))
  }

  return(0)
}


#' Prepares a run of the enumpart algorithm by ordering edges
#'
#' @param adj zero indexed adjacency list
#' @param unordered_path valid path to output the unordered adjacency map to
#' @param ordered_path valid path to output the ordered adjacency map to
#' @param weight_path A path (not including ".dat") to store a space-delimited
#' file containing a vector of vertex weights. Only supply with total_pop.
#' @param total_pop the vector of precinct populations. Only supply with weight_path
#'
#'
#' @return 0 on success
#' @export
#' @importFrom sys exec_wait
#'
#' @references
#' Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
#' "The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
#' Forthcoming, Statistics and Public Policy.
#' @concept enumerate
#' @examples \dontrun{
#' temp <- tempdir()
#' data(fl25)
#' adj <- redist.adjacency(fl25)
#' redist.prep.enumpart(adj = adj, unordered_path = paste0(temp, '/unordered'),
#'                      ordered_path = paste0(temp, '/ordered'))
#' }
redist.prep.enumpart <- function(adj, unordered_path, ordered_path,
                                 weight_path = NULL, total_pop = NULL){

  if (is.null(weight_path) & !is.null(total_pop)) {
    stop('`weight_path` not null, but `total_pop` is. Provide both or none.')
  }
  if (!is.null(weight_path) & is.null(total_pop)) {
    stop('`total_pop` not null, but `weight_path` is. Provide both or none.')
  }


  # Return the list to 1 indexing
  adj <- lapply(adj, function(x){x+1})

  # Remove any duplicates:
  adj <- lapply(adj, unique)

  ## Sink
  adj_map <- c()
  for(k in 1:length(adj)){
    sub <- adj[[k]]
    sub <- sub[sub > k]
    if(length(sub) > 0){
      for(l in 1:length(sub)){
        adj_map <- rbind(adj_map, c(k, sub[l]))
      }
    }
  }

  utils::write.table(data.frame(adj_map), file = paste0(unordered_path,".dat"),
                     quote=FALSE, row.names=FALSE, col.names=FALSE)

  ## Order edges

  if(Sys.info()[['sysname']] == 'Windows'){
    res <- sys::exec_wait('python',
                          args = system.file('python/ndscut.py', package = 'redist'),
                          std_in = paste0(unordered_path, '.dat'),
                          std_out = paste0(ordered_path, '.dat'))
  } else {
    res <- sys::exec_wait('python3',
                          args = system.file('python/ndscut.py', package = 'redist'),
                          std_in = paste0(unordered_path, '.dat'),
                          std_out = paste0(ordered_path, '.dat'))
  }

  if (!is.null(weight_path)) {
    utils::write.table(t(total_pop), file = paste0(weight_path,".dat"),
                       quote=FALSE, row.names=FALSE, col.names=FALSE)
  }

  return(res)
}

#' Runs the enumpart algorithm
#'
#' @param ordered_path Path used in redist.prep.enumpart (not including ".dat")
#' @param out_path Valid path to output the enumerated districts
#' @param ndists number of districts to enumerate
#' @param all boolean. TRUE outputs all districts. FALSE samples n districts.
#' @param n integer. Number of districts to output if all is FALSE. Returns
#' districts selected from uniform random distribution.
#' @param weight_path A path (not including ".dat") to a space-delimited file containing a vector of
#' vertex weights, to be used along with \code{lower} and \code{upper}.
#' @param lower A lower bound on each partition's total weight, implemented by rejection sampling.
#' @param upper An upper bound on each partition's total weight.
#' @param options Additional enumpart arguments. Not recommended for use.
#'
#' @references
#' Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
#' "The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
#' Forthcoming, Statistics and Public Policy.
#'
#' @return 0 on success
#' @export
#' @concept enumerate
#'
#' @examples \dontrun{
#' temp <- tempdir()
#' redist.run.enumpart(ordered_path = paste0(temp, '/ordered'),
#' out_path = paste0(temp, '/enumerated'))
#' }
redist.run.enumpart <- function(ordered_path, out_path, ndists = 2,
                                all = TRUE, n  = NULL, weight_path = NULL,
                                lower = NULL, upper = NULL, options = NULL){
  ndists <- as.integer(ndists)
  n <- as.integer(n)

  # use args based on types
  if (is.null(options)) {
      if (all) {
          options <- c('-k', ndists, '-comp', '-allsols')
      } else{
          if (is.null(n)) {
            stop('n must be specified when all is FALSE.')
          }
          options <- c('-k', ndists, '-comp', '-sample', n)
      }
  }

  if (!is.null(lower)) {
    options <-  c(options, "-lower", as.character(lower))
  }
  if (!is.null(upper)) {
    options = c(options, "-upper", as.character(upper))
  }

  if (is.null(weight_path)) {
      options <- c(paste0(ordered_path, '.dat'), options)
  } else {
      options <- c(paste0(ordered_path, '.dat'), paste0(weight_path, ".dat"), options)
  }

  ## Run enumpart
  res <- sys::exec_wait(paste0(system.file('enumpart', package = 'redist'), '/enumpart'),
                 args = options,
                 std_out = paste0(out_path, '.dat'), std_err = TRUE)

  return(res)
}




#' Read Results from enumpart
#'
#' @param out_path out_path specified in redist.run.enumpart
#' @param skip number of lines to skip
#' @param n_max max number of lines to read
#'
#' @return district_membership matrix
#' @export
#' @references
#' Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
#' "The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
#' Forthcoming, Statistics and Public Policy.
#'
#' @importFrom readr read_lines
#' @concept enumerate
#' @examples \dontrun{
#' temp <- tempdir()
#' cds <- redist.read.enumpart(out_path = paste0(temp,'/enumerated'))
#' }
redist.read.enumpart <- function(out_path, skip = 0,  n_max = -1L){
  sols <- readr::read_lines(paste0(out_path, ".dat"), skip = skip,
                            n_max = n_max, lazy = FALSE)
  sols <- apply(do.call("cbind", strsplit(sols, " ")), 2, as.numeric)
  return(sols + 1L)
}


# check if last edge
#
# @param i integer, current frontier
# @param v integer, vertex to search for
# @param edges edgelist matrix
#
# @return bool
#
is_last <- function(i, v, edges){
  if(i == nrow(edges)){
    return(TRUE)
  }
  for(j in (i+1):nrow(edges)){
    if(v ==  edges[j, 1] | v == edges[j, 2]){
      return(FALSE)
    }
  }
  return(TRUE)
}


#' Calculate Frontier Size
#'
#' @param ordered_path path to ordered path created by redist.prep.enumpart
#'
#' @return List, four objects
#' \itemize{
#' \item{max}{numeric, maximum frontier size}
#' \item{average}{numeric, average frontier size}
#' \item{average_sq}{numeric, average((frontier size)^2)}
#' \item{sequence}{numeric vector, lists out all sizes for every frontier}
#' }
#' @export
#' @concept enumerate
#'
#' @importFrom stringr str_split
#' @examples \dontrun{
#' data(fl25)
#' adj <- redist.adjacency(fl25)
#' redist.prep.enumpart(adj, 'unordered', 'ordered')
#' redist.calc.frontier.size('ordered')
#' }
redist.calc.frontier.size <- function(ordered_path){
  lines_in <- readLines(paste0(ordered_path,'.dat'))
  n <- length(lines_in)

  edges_unsort <- apply(stringr::str_split(string = lines_in, pattern = ' ', simplify = T),2, as.integer)
  edges <- cbind(apply(edges_unsort,1,min), apply(edges_unsort,1,max))

  frontier_sizes <- rep(NA_real_, 1 +n)
  frontier <- rep(FALSE, n)
  frontier_sizes[1] <- 0

  for(i in 1:n){
    e1 <- edges[i,1]
    e2 <- edges[i,2]
    frontier[e1] <- TRUE
    frontier[e2] <- TRUE

    if(is_last(i, e1, edges)){
      frontier[e1] <- FALSE
    }
    if(is_last(i, e2, edges)){
      frontier[e2] <- FALSE
    }

    frontier_sizes[i+1] <- sum(frontier)
  }


  return(
    list(max = max(frontier_sizes),
         average = mean(frontier_sizes),
         average_sq = mean(frontier_sizes^2),
         sequence = frontier_sizes)
  )
}

#' Enumerate All Parititions
#'
#' Single function for standard enumeration analysis.
#'
#' @param adj zero indexed adjacency list.
#' @param unordered_path valid path to output the unordered adjacency map to
#' @param ordered_path valid path to output the ordered adjacency map to
#' @param out_path Valid path to output the enumerated districts
#' @param ndists number of districts to enumerate
#' @param all boolean. TRUE outputs all districts. FALSE samples n districts.
#' @param n integer. Number of districts to output if all is FALSE. Returns
#' districts selected from uniform random distribution.
#' @param weight_path A path (not including ".dat") to a space-delimited file containing a vector of
#' vertex weights, to be used along with \code{lower} and \code{upper}.
#' @param lower A lower bound on each partition's total weight, implemented by rejection sampling.
#' @param upper An upper bound on each partition's total weight.
#' @param init Runs redist.init.enumpart. Defaults to false. Should be run on first use.
#' @param read boolean. Defaults to TRUE. reads
#' @param total_pop the vector of precinct populations
#'
#' @return List with entries district_membership and parity.
#'
#' @concept enumerate
#' @export
redist.enumpart <- function(adj, unordered_path, ordered_path,
                            out_path, ndists = 2, all = TRUE, n = NULL,
                            weight_path=NULL, lower=NULL, upper=NULL,
                            init = FALSE, read = TRUE, total_pop = NULL){
  if(init){
    redist.init.enumpart()
  }

  prep <- redist.prep.enumpart(adj = adj,
                               unordered_path = unordered_path,
                               ordered_path = ordered_path,
                               weight_path = weight_path,
                               total_pop = total_pop)
  if(!prep){
    run <- redist.run.enumpart(ordered_path = ordered_path,
                               out_path = out_path,
                               ndists = ndists,
                               all = all,
                               n = n,
                               weight_path = weight_path,
                               lower = lower,
                               upper = upper)
  }

  if(read){
    cds <- redist.read.enumpart(out_path = out_path)
    if(!is.null(total_pop)){
      par <- redist.parity(plans = cds, total_pop = total_pop)
    } else{
      par <- rep(NA_real_, ncol(cds))
    }
    out <- list(plans = cds, parity = par)
  } else{
    return(0)
  }

  return(out)

}
