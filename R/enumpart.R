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
  servr::make(dir = system.file('enumpart', package = 'redist'))
  sys::exec_wait('python', args= c('-m', 'pip', 'install', 'networkx', '--user'))
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
                     path = paste0(unordered_path,".dat"),
                     col_names = FALSE)
  
  ## Order edges
  res <- sys::exec_wait('python', 
                        args = system.file('python/ndscut.py', package = 'redist'), 
                        std_in = paste0(unordered_path, '.dat'), 
                        std_out = paste0(ordered_path, '.dat'))  
  
  return(res)
}

#' Runs the enumpart algorithm
#'
#' @param ordered_path Path used in redist.prep.enumpart
#' @param out_path Valid path to output the enumerated districts
#' @param ndist number of districts to enumerate
#' @param all boolean. TRUE outputs all districts. FALSE samples n districts. 
#' @param n integer. Number of districts to output if all is FALSE. Returns 
#' districts selected from uniform random distribution.
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
redist.run.enumpart <- function(ordered_path, out_path, ndist = 2, all = TRUE, n  = NULL, options = NULL){
  ndist <- as.integer(ndist)
  n <- as.integer(n)
  
  # use args based on types
  if(!is.null(options)){
  if(all){
    options <- c('-k', ndist, '-comp', '-allsols')
  } else{
    if(is.null(n))
      stop('n must be specified when all is FALSE.')
    options <- c('-k', ndist, '-comp', '-sample', n)
  }
  options <- c(paste0(ordered_path, '.dat'), options)
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