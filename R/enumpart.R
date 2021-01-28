#' Initialize enumpart
#'
#'This ensures that the enumerate partitions programs is prepared to run.
#'This must be run once per install of the redist package.
#'
#' @return 0 on success
#' @export
#' @importFrom servr make
#' @references
#' Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
#' "The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
#' Forthcoming, Statistics and Public Policy.
#'
#' @examples \dontrun{
#' redist.init.enumpart()
#' }
redist.init.enumpart <- function(){
  # Update makefile to direct to library only if Windows
  if(Sys.info()[['sysname']] == 'Windows'){
    makecontent <- readLines(system.file('enumpart/Makefile', package = 'redist'))
    makecontent[7] <-"\tg++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG -lpsapi"
    writeLines(text = makecontent, con = system.file('enumpart/Makefile', package = 'redist'))
  }

  servr::make(dir = system.file('enumpart', package = 'redist'), verbose = FALSE)
  sys::exec_wait('python3', args= c('-m', 'pip', 'install', 'networkx', '--user'))
  return(0)
}


#' Prepares a run of the enumpart algorithm by ordering edges
#'
#' @param adjlist zero indexed adjacency list
#' @param unordered_path valid path to output the unordered adjacency map to
#' @param ordered_path valid path to output the ordered adjacency map to
#'
#' @return 0 on success
#' @export
#' @importFrom readr write_delim
#' @importFrom sys exec_wait
#'
#' @references
#' Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
#' "The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
#' Forthcoming, Statistics and Public Policy.
#' @examples \dontrun{
#' data("algdat.p10")
#' redist.prep.enumpart(adjlist = algdat.p10$adjlist, unordered_path = '../unordered',
#' ordered_path = '../ordered')
#' }
redist.prep.enumpart <- function(adjlist, unordered_path, ordered_path){
  # Return the list to 1 indexing
  adjlist <- lapply(adjlist, function(x){x+1})

  ## Sink
  adjlist_map <- c()
  for(k in 1:length(adjlist)){
    sub <- adjlist[[k]]
    sub <- sub[sub > k]
    if(length(sub) > 0){
      for(l in 1:length(sub)){
        adjlist_map <- rbind(adjlist_map, c(k, sub[l]))
      }
    }
  }

  readr::write_delim(data.frame(adjlist_map),
                     file = paste0(unordered_path,".dat"),
                     col_names = FALSE)

  ## Order edges
  res <- sys::exec_wait('python3',
                        args = system.file('python/ndscut.py', package = 'redist'),
                        std_in = paste0(unordered_path, '.dat'),
                        std_out = paste0(ordered_path, '.dat'))

  return(res)
}

#' Runs the enumpart algorithm
#'
#' @param ordered_path Path used in redist.prep.enumpart (not including ".dat")
#' @param out_path Valid path to output the enumerated districts
#' @param ndist number of districts to enumerate
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
#'
#' @examples \dontrun{
#' redist.run.enumpart(ordered_path = '../ordered', out_path = '../enumerated')
#' }
redist.run.enumpart <- function(ordered_path, out_path, ndist=2, all=TRUE, n=NULL,
                                weight_path=NULL, lower=NULL, upper=NULL, options=NULL){
  ndist <- as.integer(ndist)
  n <- as.integer(n)

  # use args based on types
  if (all) {
      options <- c('-k', ndist, '-comp', '-allsols', options)
  } else {
      if(is.null(n))
          stop('n must be specified when all is FALSE.')
      options <- c('-k', ndist, '-comp', '-sample', n, options)
  }

  if (!is.null(lower)) options = c(options, "-lower", as.character(lower))
  if (!is.null(upper)) options = c(options, "-upper", as.character(upper))

  if (is.null(weight_path)) {
      options <- c(paste0(ordered_path, '.dat'), options)
  } else {
      options <- c(paste0(ordered_path, '.dat'), paste0(weight_path, ".dat"), options)
  }

  ## Run enumpart
  res <- sys::exec_wait(system.file('enumpart/enumpart', package = 'redist'),
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
#' @importFrom readr read_lines
#' @references
#' Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
#' "The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
#' Forthcoming, Statistics and Public Policy.
#'
#' @examples \dontrun{
#' cds <- redist.read.enumpart(out_path = '../enumerated')
#' }
redist.read.enumpart <- function(out_path, skip = 0,  n_max = -1L){
  sols <- readr::read_lines(paste0(out_path, ".dat"))#, skip = skip, n_max = n_max)
  sols <- apply(do.call("cbind", strsplit(sols, " ")), 2, as.numeric)
  return(sols)
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
#' @param adjlist zero indexed adjacency list
#' @param unordered_path valid path to output the unordered adjacency map to
#' @param ordered_path valid path to output the ordered adjacency map to
#' @param out_path Valid path to output the enumerated districts
#' @param ndist number of districts to enumerate
#' @param all boolean. TRUE outputs all districts. FALSE samples n districts.
#' @param n integer. Number of districts to output if all is FALSE. Returns
#' districts selected from uniform random distribution.
#' @param weight_path A path (not including ".dat") to a space-delimited file containing a vector of
#' vertex weights, to be used along with \code{lower} and \code{upper}.
#' @param lower A lower bound on each partition's total weight, implemented by rejection sampling.
#' @param upper An upper bound on each partition's total weight.
#' @param init Runs redist.init.enumpart. Defaults to false. Should be run on first use.
#' @param read boolean. Defaults to TRUE. reads
#' @param population Integer Vector. Defaults to NULL. If supplied, computes the parity.
#'
#' @return
#' @export
#'
redist.enumpart <- function(adjlist, unordered_path, ordered_path, out_path,
                            ndist=2, all=TRUE, n=NULL, init=FALSE,
                            weight_path=NULL, lower=NULL, upper=NULL,
                            read=TRUE, population=NULL) {
  if(init){
    redist.init.enumpart()
  }

  prep <- redist.prep.enumpart(adjlist = adjlist,
                               unordered_path = unordered_path,
                               ordered_path = ordered_path)
  if(!prep){
    run <- redist.run.enumpart(ordered_path = ordered_path,
                               out_path = out_path,
                               ndist = ndist,
                               all = all,
                               n = n,
                               weight_path = weight_path,
                               lower = lower,
                               upper = upper)
  }

  if(read){
    cds <- redist.read.enumpart(out_path = out)
    if(!is.null(population)){
      par <- redist.parity(district_membership = cds, population = population)
    } else{
      par <- rep(NA_real_, ncol(cds))
    }
    out <- list(district_membership = cds, parity = par)
  } else{
    return(0)
  }

  return(out)

}
