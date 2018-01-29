## Deprecated code from the geiger r package
area.between.curves <- function(x, f1, f2){
    xrange = c(min(x), max(x))
    a<-0.0;
    for(i in 1:length(x)) {
        if(x[i]>=xrange[1] & x[i]<=xrange[2]) {
            if(i==1) {
                lhs<-0
            } else if(x[i-1]<xrange[1]) {
                lhs<-xrange[1]
            } else lhs<-x[i-1];
            if(i==length(x)) {
                rhs<-x[i]
            } else if(x[i+1]>xrange[2]) {
                rhs<-xrange[2];
            } else rhs<-x[i+1];
            a<-a+(f2[i]-f1[i])*(rhs-lhs)/2;
        } else if(i!=1)
            if(x[i-1]>=xrange[1] & x[i-1]<=xrange[2]) {
                y1<-f1[i-1]+(f1[i]-f1[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
                y2<-f2[i-1]+(f2[i]-f2[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
                a<-a+(y2-y1)*(xrange[2]-x[i-1])/2;
            } else if(i!=length(x))
                if(x[i+1]>=xrange[1] & x[i+1]<=xrange[2]) {
                    y1<-f1[i]+(f1[i+1]-f1[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
                    y2<-f2[i]+(f2[i+1]-f2[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
                    a<-a+(y2-y1)*(x[i+1]-xrange[1])/2;
                }
    }
    return(a)
}

#' Calculate standard redistricting diagnostic statistics
#'
#' \code{calc.stats} can calculate various statistics used to
#' summarize ensembles of simulated redistricting plans.
#'
#' @usage calc.stats(algout, group1vote, group2vote, margin, swing, npoints, stats)
#'
#' @param algout An object of class "redist" or a matrix of congressional district
#' assignments where columns are different redistricting plans.
#' @param group1vote Vote counts for group 1 (for instance, Democrats)
#' @param group2vote Vote counts for group 2 (for instance, Republicans)
#' @param margin_mc The margin (50% - margin, 50% + margin) for counting "close" plans
#' for the "mc" stat.
#' @param margin_ps The margin (50% - swing, 50% + swing) over which to calculate the
#' partisan symmetry statistic.
#' @param npoints The number of points along which to estimate the integral for the
#' partisan symmetry statistic.
#' @param stats The stats to calculate. "ec" refers to electoral competitiveness,
#' "mc" is a count of plans within a particular vote margin, "ps_simple" is a
#' simple measure of partisan symmetry that counts the expected number of seats
#' given a party vote share minus the expected number of seats that party would receive
#' under the other party's vote share, "eg" is the efficiency gap, and "ps_full" is a
#' full implementation of the partisan symmetry diagnostic that measures the area
#' between the empirical seats-votes curve and a null, balanced seats-votes curve.
#' @param nc Number of cores to parallelize calculation over. Default is detectCores() - 1.
#'
#' @return \code{calc.stats} returns a list of estimated statistics for each plan.
#'
#' @export
#' @importFrom parallel mclapply detectCores
calc.stats <- function(algout, group1vote, group2vote,
                       margin_mc = .05, margin_ps = .05, npoints = 10,
                       stats = c("ec", "mc", "ps_simple", "eg", "ps_full"),
                       nc = NULL){

    ## Check
    if(sum(stats %in% c("ec", "mc", "ps_simple", "eg", "ps_full")) < length(stats)){
        badstat <- stats[!(stats %in% c("ec", "mc", "ps_simple", "eg", "ps_full"))]
        stop("The specified statistics", badstat, "are not supported.")
    }
    if(inherits(algout, "redist")){
        algout <- algout$partitions
    }
    if(length(group1vote) != nrow(algout)){
        stop("The number of geographic units in algout is not equal to the length of group1vote.")
    }
    if(length(group2vote) != nrow(algout)){
        stop("The number of geographic units in algout is not equal to the length of group2vote.")
    }
    if(is.null(nc)){
        nc <- detectCores()-1
    }

    ## Container
    stat_out <- vector(mode = "list", length = length(stats))
    ind <- 1

    ## Electoral competitiveness
    if("ec" %in% stats){
        cat("Calculating electoral competitiveness.\n")
        t_p <- mclapply(
            1:ncol(algout), function(i){
                x <- algout[,i]
                unq <- length(unique(x))
                return(sum(abs(tapply(group1vote, x, sum) /
                               (tapply(group1vote, x, sum) +
                                tapply(group2vote, x, sum)) - .5)) / unq)
            }, mc.cores = nc
        )
        t_e <- mclapply(
            1:ncol(algout), function(i){
                x <- algout[,i]
                unq <- length(unique(x))
                n_d <- tapply(group1vote, x, sum)
                n_r <- tapply(group2vote, x, sum)
                b_r <- sum(n_d > n_r)
                return(abs(b_r / unq - .5))
            }, mc.cores = nc
        )
        alpha <- 1
        beta <- 4/3
        f <- unlist(t_p) * (1 + alpha * unlist(t_e)) * beta

        stat_out[[ind]] <- f
        ind <- ind + 1
        
    }

    ## Close seats
    if("mc" %in% stats){
        cat("Counting number of seats within specified margin.\n")
        closeseats <- mclapply(
            1:ncol(algout), function(i){
                x <- algout[,i]
                demvote <- tapply(group1vote, x, sum)
                repvote <- tapply(group2vote, x, sum)
                return(sum(demvote / (demvote + repvote) > .5 - margin_mc &
                           demvote / (demvote + repvote) < .5 + margin_mc))
            }, mc.cores = nc
        )

        stat_out[[ind]] <- unlist(closeseats)
        ind <- ind + 1

    }

    ## Partisan Symmetry - simple
    if("ps_simple" %in% stats){
        cat("Calculating partisan symmetry by seats-votes inversion.\n")
        
        ## Set up swing
        statebase <- sum(group1vote) /
            (sum(group1vote) + sum(group2vote))
        swingshare <- statebase - (1 - statebase)
        swing <- (group1vote + group2vote) * abs(swingshare)

        ## Induce swing
        group1vote_swing <- group1vote - swing
        group2vote_swing <- group2vote + swing

        ps <- mclapply(
            1:ncol(algout), function(i){
                x <- algout[,i]
                demvote <- tapply(group1vote, x, sum)
                repvote <- tapply(group2vote, x, sum)
                f_v0 <- sum(demvote > repvote) / length(demvote)

                demvote_swing <- tapply(group1vote_swing, x, sum)
                repvote_swing <- tapply(group2vote_swing, x, sum)
                f_v0p <- sum(demvote_swing > repvote_swing) / length(demvote)

                return(abs(f_v0 - (1 - f_v0p)))
            }, mc.cores = nc
        )

        stat_out[[ind]] <- unlist(ps)
        ind <- ind + 1
        
    }

    ## Efficiency gap
    if("eg" %in% stats){
        cat("Calculating efficiency gap.\n")
        egap <- mclapply(
            1:ncol(algout), function(i){
                x <- algout[,i]
                demvote <- tapply(group1vote, x, sum)
                repvote <- tapply(group2vote, x, sum)
                seatmarg <- sum(demvote > repvote) - .5
                votemarg <- sum(demvote) / (sum(demvote) + sum(repvote)) - .5
                return(seatmarg - 2 * votemarg)
            }, mc.cores = nc
        )

        stat_out[[ind]] <- unlist(egap)
        ind <- ind + 1
        
    }

    ## Partisan symmetry - full distribution
    if("ps_full" %in% stats){
        cat("Calculating partisan symmetry with full seats-votes curve. May take a while.\n")
        ## Set up swing
        statebase <- sum(group1vote) /
            (sum(group1vote) + sum(group2vote))
        equal <- .5 - statebase

        v_s <- margin_ps
        inc <- seq(0, v_s, length.out = npoints)

        bias <- unique(c(rev(equal - inc), equal + inc))

        null_df <- data.frame(voteshare = c(0, .5, 1),
                              seatshare = c(0, length(unique(algout[,1]))/2,
                                            length(unique(algout[,1]))))
        lm_null <- lm(seatshare ~ voteshare, null_df)

        voteshare <- statebase + bias
        pred_seatshare <- predict(lm_null, data.frame(voteshare))

        ## Calculate
        pbias <- mclapply(
            1:ncol(algout), function(i){
                x <- algout[,i]

                ## Get seats
                nseats <- rep(NA, length(bias))
                for(j in 1:length(bias)){
                    swing <- (group1vote + group2vote) * abs(bias[j])
                    if(bias[j] < 0){
                        demvote_swing <- group1vote - swing
                        repvote_swing <- group2vote + swing
                    }else{
                        demvote_swing <- group1vote + swing
                        repvote_swing <- group2vote - swing
                    }
                    nseats[j] <- sum(tapply(demvote_swing, x, sum) >
                                     tapply(repvote_swing, x, sum))
                }

                ## Calc area
                gt0_area <- area.between.curves(
                    voteshare[voteshare >= .5],
                    nseats[voteshare >= .5],
                    pred_seatshare[voteshare >= .5]
                )
                lt0_area <- area.between.curves(
                    voteshare[voteshare <= .5],
                    nseats[voteshare <= .5],
                    pred_seatshare[voteshare <= .5]
                )

                return(gt0_area - lt0_area)
            }, mc.cores = nc
        )

        stat_out[[ind]] <- unlist(pbias)

    }

    stats <- stats[match(c("ec", "mc", "ps_simple", "eg", "ps_full"), stats)]
    stats <- stats[!is.na(stats)]
    
    names(stat_out) <- stats
    return(stat_out)
    
}
