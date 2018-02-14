## ## Deprecated code from the geiger r package
## area.between.curves <- function(x, f1, f2){
##     xrange = c(min(x), max(x))
##     a<-0.0;
##     for(i in 1:length(x)) {
##         if(x[i]>=xrange[1] & x[i]<=xrange[2]) {
##             if(i==1) {
##                 lhs<-0
##             } else if(x[i-1]<xrange[1]) {
##                 lhs<-xrange[1]
##             } else lhs<-x[i-1];
##             if(i==length(x)) {
##                 rhs<-x[i]
##             } else if(x[i+1]>xrange[2]) {
##                 rhs<-xrange[2];
##             } else rhs<-x[i+1];
##             a<-a+(f2[i]-f1[i])*(rhs-lhs)/2;
##         } else if(i!=1)
##             if(x[i-1]>=xrange[1] & x[i-1]<=xrange[2]) {
##                 y1<-f1[i-1]+(f1[i]-f1[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
##                 y2<-f2[i-1]+(f2[i]-f2[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
##                 a<-a+(y2-y1)*(xrange[2]-x[i-1])/2;
##             } else if(i!=length(x))
##                 if(x[i+1]>=xrange[1] & x[i+1]<=xrange[2]) {
##                     y1<-f1[i]+(f1[i+1]-f1[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
##                     y2<-f2[i]+(f2[i+1]-f2[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
##                     a<-a+(y2-y1)*(x[i+1]-xrange[1])/2;
##                 }
##     }
##     return(a)
## }

## #' Calculate standard redistricting diagnostic statistics
## #'
## #' \code{redist.calcstats} can calculate various statistics used to
## #' summarize ensembles of simulated redistricting plans.
## #'
## #' @usage redist.calcstats(algout, group1vote, group2vote, margin, swing, npoints, stats)
## #'
## #' @param algout An object of class "redist" or a matrix of congressional district
## #' assignments where columns are different redistricting plans.
## #' @param group1vote Vote counts for group 1 (for instance, Democrats)
## #' @param group2vote Vote counts for group 2 (for instance, Republicans)
## #' @param margin_mc The margin (50% - margin, 50% + margin) for counting "close" plans
## #' for the "mc" stat.
## #' @param margin_ps The margin (50% - swing, 50% + swing) over which to calculate the
## #' partisan symmetry statistic.
## #' @param npoints The number of points along which to estimate the integral for the
## #' partisan symmetry statistic.
## #' @param stats The stats to calculate. "ec" refers to electoral competitiveness,
## #' "mc" is a count of plans within a particular vote margin, "ps_simple" is a
## #' simple measure of partisan symmetry that counts the expected number of seats
## #' given a party vote share minus the expected number of seats that party would receive
## #' under the other party's vote share, "eg" is the efficiency gap, and "ps_full" is a
## #' full implementation of the partisan symmetry diagnostic that measures the area
## #' between the empirical seats-votes curve and a null, balanced seats-votes curve.
## #' @param nc Number of cores to parallelize calculation over. Default is detectCores() - 1.
## #'
## #' @return \code{redist.calcstats} returns a list of estimated statistics for each plan.
## #'
## #' @export
## #' @importFrom parallel mclapply detectCores
## redist.calcstats <- function(algout, group1vote, group2vote,
##                        margin_mc = .05, margin_ps = .05, npoints = 10,
##                        stats = c("ec", "mc", "ps_simple", "eg", "ps_full"),
##                        nc = NULL){

##     ## Check
##     if(sum(stats %in% c("ec", "mc", "ps_simple", "eg", "ps_full")) < length(stats)){
##         badstat <- stats[!(stats %in% c("ec", "mc", "ps_simple", "eg", "ps_full"))]
##         stop("The specified statistics", badstat, "are not supported.")
##     }
##     if(inherits(algout, "redist")){
##         algout <- algout$partitions
##     }
##     if(length(group1vote) != nrow(algout)){
##         stop("The number of geographic units in algout is not equal to the length of group1vote.")
##     }
##     if(length(group2vote) != nrow(algout)){
##         stop("The number of geographic units in algout is not equal to the length of group2vote.")
##     }
##     if(is.null(nc)){
##         nc <- detectCores()-1
##     }

##     ## Container
##     stat_out <- vector(mode = "list", length = length(stats))
##     ind <- 1

##     ## Electoral competitiveness
##     if("ec" %in% stats){
##         cat("Calculating electoral competitiveness.\n")
##         t_p <- mclapply(
##             1:ncol(algout), function(i){
##                 x <- algout[,i]
##                 unq <- length(unique(x))
##                 return(sum(abs(tapply(group1vote, x, sum) /
##                                (tapply(group1vote, x, sum) +
##                                 tapply(group2vote, x, sum)) - .5)) / unq)
##             }, mc.cores = nc
##         )
##         t_e <- mclapply(
##             1:ncol(algout), function(i){
##                 x <- algout[,i]
##                 unq <- length(unique(x))
##                 n_d <- tapply(group1vote, x, sum)
##                 n_r <- tapply(group2vote, x, sum)
##                 b_r <- sum(n_d > n_r)
##                 return(abs(b_r / unq - .5))
##             }, mc.cores = nc
##         )
##         alpha <- 1
##         beta <- 4/3
##         f <- unlist(t_p) * (1 + alpha * unlist(t_e)) * beta

##         stat_out[[ind]] <- f
##         ind <- ind + 1
        
##     }

##     ## Close seats
##     if("mc" %in% stats){
##         cat("Counting number of seats within specified margin.\n")
##         closeseats <- mclapply(
##             1:ncol(algout), function(i){
##                 x <- algout[,i]
##                 demvote <- tapply(group1vote, x, sum)
##                 repvote <- tapply(group2vote, x, sum)
##                 return(sum(demvote / (demvote + repvote) > .5 - margin_mc &
##                            demvote / (demvote + repvote) < .5 + margin_mc))
##             }, mc.cores = nc
##         )

##         stat_out[[ind]] <- unlist(closeseats)
##         ind <- ind + 1

##     }

##     ## Partisan Symmetry - simple
##     if("ps_simple" %in% stats){
##         cat("Calculating partisan symmetry by seats-votes inversion.\n")
        
##         ## Set up swing
##         statebase <- sum(group1vote) /
##             (sum(group1vote) + sum(group2vote))
##         swingshare <- statebase - (1 - statebase)
##         swing <- (group1vote + group2vote) * abs(swingshare)

##         ## Induce swing
##         group1vote_swing <- group1vote - swing
##         group2vote_swing <- group2vote + swing

##         ps <- mclapply(
##             1:ncol(algout), function(i){
##                 x <- algout[,i]
##                 demvote <- tapply(group1vote, x, sum)
##                 repvote <- tapply(group2vote, x, sum)
##                 f_v0 <- sum(demvote > repvote) / length(demvote)

##                 demvote_swing <- tapply(group1vote_swing, x, sum)
##                 repvote_swing <- tapply(group2vote_swing, x, sum)
##                 f_v0p <- sum(demvote_swing > repvote_swing) / length(demvote)

##                 return(abs(f_v0 - (1 - f_v0p)))
##             }, mc.cores = nc
##         )

##         stat_out[[ind]] <- unlist(ps)
##         ind <- ind + 1
        
##     }

##     ## Efficiency gap
##     if("eg" %in% stats){
##         cat("Calculating efficiency gap.\n")
##         egap <- mclapply(
##             1:ncol(algout), function(i){
##                 x <- algout[,i]
##                 demvote <- tapply(group1vote, x, sum)
##                 repvote <- tapply(group2vote, x, sum)
##                 seatmarg <- sum(demvote > repvote) - .5
##                 votemarg <- sum(demvote) / (sum(demvote) + sum(repvote)) - .5
##                 return(seatmarg - 2 * votemarg)
##             }, mc.cores = nc
##         )

##         stat_out[[ind]] <- unlist(egap)
##         ind <- ind + 1
        
##     }

##     ## Partisan symmetry - full distribution
##     if("ps_full" %in% stats){
##         cat("Calculating partisan symmetry with full seats-votes curve. May take a while.\n")
##         ## Set up swing
##         statebase <- sum(group1vote) /
##             (sum(group1vote) + sum(group2vote))
##         equal <- .5 - statebase

##         v_s <- margin_ps
##         inc <- seq(0, v_s, length.out = npoints)

##         bias <- unique(c(rev(equal - inc), equal + inc))

##         null_df <- data.frame(voteshare = c(0, .5, 1),
##                               seatshare = c(0, length(unique(algout[,1]))/2,
##                                             length(unique(algout[,1]))))
##         lm_null <- lm(seatshare ~ voteshare, null_df)

##         voteshare <- statebase + bias
##         pred_seatshare <- predict(lm_null, data.frame(voteshare))

##         ## Calculate
##         pbias <- mclapply(
##             1:ncol(algout), function(i){
##                 x <- algout[,i]

##                 ## Get seats
##                 nseats <- rep(NA, length(bias))
##                 for(j in 1:length(bias)){
##                     swing <- (group1vote + group2vote) * abs(bias[j])
##                     if(bias[j] < 0){
##                         demvote_swing <- group1vote - swing
##                         repvote_swing <- group2vote + swing
##                     }else{
##                         demvote_swing <- group1vote + swing
##                         repvote_swing <- group2vote - swing
##                     }
##                     nseats[j] <- sum(tapply(demvote_swing, x, sum) >
##                                      tapply(repvote_swing, x, sum))
##                 }

##                 ## Calc area
##                 gt0_area <- area.between.curves(
##                     voteshare[voteshare >= .5],
##                     nseats[voteshare >= .5],
##                     pred_seatshare[voteshare >= .5]
##                 )
##                 lt0_area <- area.between.curves(
##                     voteshare[voteshare <= .5],
##                     nseats[voteshare <= .5],
##                     pred_seatshare[voteshare <= .5]
##                 )

##                 return(gt0_area - lt0_area)
##             }, mc.cores = nc
##         )

##         stat_out[[ind]] <- unlist(pbias)

##     }

##     stats <- stats[match(c("ec", "mc", "ps_simple", "eg", "ps_full"), stats)]
##     stats <- stats[!is.na(stats)]
    
##     names(stat_out) <- stats
##     return(stat_out)
    
## }

## #' Polsby-Popper calculation for MCMC
## #' @export
## redist.polsbypopper <- function(algout, adj_list, areas_vec, borderlength_mat, pop_vec, discrete = FALSE){

##     ## Unpack objects
##     nsims <- ncol(algout$partitions)
##     if(min(unlist(adj_list)) == 1){
##         adj_list <- lapply(adj_list, function(x){x-1})
##     }
##     partitions <- algout$partitions    
##     cds_unique <- unique(partitions[,1])

##     store_pp <- matrix(NA, nsims, 3)
##     for(i in 1:nsims){
##         cds <- partitions[,i]
##         sub_al <- genAlConn(adj_list, cds)
##         boundary_indicator <- findBoundary(adj_list, sub_al)
##         store_sim_pp <- rep(NA, length(cds_unique))
##         for(j in 1:length(cds_unique)){
##             cd_ind <- which(cds == cds_unique[j]) - 1
##             pp_out <- calc_polsbypopper(cd_ind, areas_vec, boundary_indicator,
##                                         borderlength_mat, pop_vec, adj_list,
##                                         discrete)
##             store_sim_pp[j] <- pp_out
##         }
##         store_pp[i, 1] <- mean(store_sim_pp)
## 	store_pp[i, 2] <- max(store_sim_pp)
##         store_pp[i, 3] <- min(store_sim_pp)
##     }

##     return(store_pp)
    
## }

## #' Reock calculation for MCMC
## #' @export
## redist.reock <- function(algout, areas_vec, coord_mat){

##     ## Unpack objects
##     nsims <- ncol(algout$partitions)
##     partitions <- algout$partitions
##     pwdmat <- calcPWDh(coord_mat)
##     cds_unique <- unique(partitions[,1])

##     store_reock <- matrix(NA, nsims, 2)
##     for(i in 1:nsims){
##         cds <- partitions[,i]
##         store_sim_reock <- rep(NA, length(cds_unique))
##         for(j in 1:length(cds_unique)){
##             cd_ind <- which(cds == cds_unique[j])
##             area_dist <- sum(areas_vec[cd_ind])
##             max_dist <- max(coord_mat[cd_ind, cd_ind])
##             area_circle <- pi * (max_dist / 2) ^2
##             store_sim_reock[j] <- area_dist / area_circle
##         }
##         store_reock[i, 1] <- mean(store_sim_reock)
##         store_reock[i, 2] <- max(store_sim_reock)
##         store_reock[i, 3] <- min(store_sim_reock)
##     }

##     return(store_reock)
    
## }

## #' Diagnostic plotting functionality for MCMC redistricting.
## #'
## #' \code{redist.diagplot} generates several common MCMC diagnostic plots.
## #'
## #' @usage redist.diagplot(sumstat,
## #' plot = c("trace", "autocorr", "densplot", "mean", "gelmanrubin"),
## #' logit = FALSE, savename = NULL)
## #'
## #' @param sumstat A vector, list, \code{mcmc} or \code{mcmc.list} object
## #' containing a summary statistic of choice.
## #' @param plot The type of diagnostic plot to generate: one of "trace",
## #' "autocorr", "densplot", "mean", "gelmanrubin". If \code{plot = "gelmanrubin"},
## #' the input \code{sumstat} must be of class \code{mcmc.list} or \code{list}.
## #' @param logit Flag for whether to apply the logistic transformation for the
## #' summary statistic. The default is \code{FALSE}.
## #' @param savename Filename to save the plot. Default is \code{NULL}.
## #'
## #' @details This function allows users to generate several standard diagnostic
## #' plots from the MCMC literature, as implemented by Plummer et. al (2006).
## #' Diagnostic plots implemented include trace plots, autocorrelation plots,
## #' density plots, running means, and Gelman-Rubin convergence diagnostics
## #' (Gelman & Rubin 1992).
## #'
## #' @return Returns a plot of file type \code{.pdf}.
## #'
## #' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
## #' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
## #' Carlo." Working Paper. Available at
## #' \url{http://imai.princeton.edu/research/files/redist.pdf}.
## #'
## #' Gelman, Andrew and Donald Rubin. (1992) "Inference from iterative simulations
## #' using multiple sequences (with discussion)." Statistical Science.
## #'
## #' Plummer, Martin, Nicky Best, Kate Cowles and Karen Vines. (2006) "CODA:
## #' Convergence Diagnosis and Output Analysis for MCMC." R News.
## #'
## #' @examples
## #' \dontrun{
## #' data(algdat.pfull)
## #'
## #' ## Get an initial partition
## #' set.seed(1)
## #' initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]
## #'
## #' ## 25 precinct, three districts - no pop constraint ##
## #' alg_253 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
## #' popvec = algdat.pfull$precinct.data$pop,
## #' initcds = initcds,nsims = 10000)
## #'
## #' ## Get Republican Dissimilarity Index from simulations
## #' rep_dmi_253 <- redist.segcalc(alg_253,
## #' algdat.pfull$precinct.data$repvote,
## #' algdat.pfull$precinct.data$pop)
## #'
## #' ## Generate diagnostic plots
## #' redist.diagplot(rep_dmi_253, plot = "trace")
## #' redist.diagplot(rep_dmi_253, plot = "autocorr")
## #' redist.diagplot(rep_dmi_253, plot = "densplot")
## #' redist.diagplot(rep_dmi_253, plot = "mean")
## #' }
## #' @export
## redist.diagplot <- function(sumstat,
##                             plot = c("trace", "autocorr", "densplot",
##                                 "mean", "gelmanrubin"),
##                             logit = FALSE, savename = NULL
##                             ){

##     ##############
##     ## Warnings ##
##     ##############
##     if(missing(sumstat)){
##         stop("Please provide a vector or list of summary statistics to the function")
##     }
##     if(!(class(sumstat) %in% c("numeric", "list", "mcmc", "mcmc.list"))){
##         stop("Please provide either a numeric vector, list, or mcmc object")
##     }
##     if(!(plot %in% c("trace", "autocorr", "densplot",
##                      "mean", "gelmanrubin"))){
##         stop("Sorry. We don't currently support that MCMC diagnostic.")
##     }
##     if(plot == "gelmanrubin" & !(class(sumstat) %in% c("list", "mcmc.list"))){
##         stop("If generating a Gelman-Rubin plot, please provide an object of class list or mcmc.list")
##     }
    
##     ########################
##     ## Create mcmc object ##
##     ########################
##     if(class(sumstat) == "numeric"){
##         segout <- mcmc(sumstat)
##     }else if(class(sumstat) == "list"){
##         for(i in 1:length(sumstat)){
##             sumstat[[i]] <- mcmc(sumstat[[i]])
##         }       
##         segout <- mcmc.list(sumstat)
##     }else if(class(sumstat) %in% c("mcmc", "mcmc.list")){
##         segout <- sumstat
##     }
    
##     ## Logit transform
##     if(logit){
##         if(class(segout) == "mcmc"){
##             segout <- log(segout / (1 - segout))
##         }else if(class(segout) == "mcmc.list"){
##             for(i in 1:length(segout)){
##                 segout[[i]] <- log(segout[[i]] / (1 - segout[[i]]))
##             }
##         }
##     }

##     ##################
##     ## Create plots ##
##     ##################
##     if(plot == "trace"){
##         if(!is.null(savename)){
##             pdf(file = paste(savename, ".pdf", sep = ""))
##         }
##         traceplot(segout)
##         if(!is.null(savename)){
##             dev.off()
##         }
##     }
##     if(plot == "autocorr"){
##         if(!is.null(savename)){
##             pdf(file = paste(savename, ".pdf", sep = ""))
##         }
##         autocorr.plot(segout, lag.max = 50)
##         if(!is.null(savename)){
##             dev.off()
##         }
##     }
##     if(plot == "densplot"){
##         if(!is.null(savename)){
##             pdf(file = paste(savename, ".pdf", sep = ""))
##         }
##         densplot(segout)
##         if(!is.null(savename)){
##             dev.off()
##         }
##     }
##     if(plot == "mean"){
##         if(!is.null(savename)){
##             pdf(file = paste(savename, ".pdf", sep = ""))
##         }
##         cumuplot(segout, probs = .5, type = "l", lty = 1)
##         if(!is.null(savename)){
##             dev.off()
##         }
##     }
##     if(plot == "gelmanrubin" & class(segout) == "mcmc.list"){
##         if(!is.null(savename)){
##             pdf(file = paste(savename, ".pdf", sep = ""))
##         }
##         gelman.plot(segout, transform = FALSE)
##         if(!is.null(savename)){
##             dev.off()
##         }
##     }

## }

#' Segregation index calculation for MCMC redistricting.
#'
#' \code{redist.segcalc} calculates the dissimilarity index of segregation (see
#' Massey \& Denton 1987 for more details) for a specified subgroup under any
#' redistricting plan.
#'
#' @usage redist.segcalc(algout, grouppop, fullpop)
#'
#' @param algout A matrix of congressional district assignments or a
#' redist object.
#' @param grouppop A vector of populations for some subgroup of interest.
#' @param fullpop A vector containign the populations of each geographic unit.
#'
#' @return \code{redist.segcalc} returns a vector where each entry is the
#' dissimilarity index of segregation (Massey & Denton 1987) for each
#' redistricting plan in \code{algout}.
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain
#' Monte Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' Massey, Douglas and Nancy Denton. (1987) "The Dimensions of Social
#' Segregation". Social Forces.
#'
#' @examples
#' \dontrun{
#' data(algdat.pfull)
#'
#' ## Code to run the simulations in Figure 4 of Fifield, Higgins,
#' ## Imai and Tarr (2015)
#'
#' ## Get an initial partition
#' set.seed(1)
#' initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]
#'
#' ## Run simulations
#' alg_253 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
#' popvec = algdat.pfull$precinct.data$pop,
#' initcds = initcds, nsims = 10000)
#'
#' ## Get Republican Dissimilarity Index from simulations
#' rep_dmi_253 <- redist.segcalc(alg_253,
#' algdat.pfull$precinct.data$repvote,
#' algdat.pfull$precinct.data$pop)
#' }
#' @export
redist.segcalc <- function(algout,
                           grouppop,
                           fullpop)
{

    ## Warnings
    if(missing(algout) | !(class(algout) %in% c("data.frame", "matrix", "redist"))){
        stop("Please provide either a redist object or a proper matrix of congessional districts")
    }
    if(missing(grouppop)){
        stop("Please provide a vector of sub-group populations to calculate
the segregation index")
    }
    if(missing(fullpop)){
        stop("Please provide a vector of populations for each geographic unit")
    }

    ## If redist object, get the partitions entry
    if(class(algout) == "redist"){
        algout <- algout$partitions
    }
    
    if(!((nrow(algout) == length(grouppop)) &
             (length(grouppop) == length(fullpop)) &
                 (length(fullpop) == nrow(algout)))){
        stop("Please make sure there is a population entry for each geographic unit")
    }

    ## Calculate dissimilarity index
    seg.out <- segregationcalc(algout,
                               grouppop,
                               fullpop)

    ## Return
    return(seg.out)
    
}
    
