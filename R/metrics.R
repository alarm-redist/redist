########################################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/07/15
## Date Modified: 2021/03/11
## Purpose: R function to compute gerrymandering metrics
########################################################

#' Calculate gerrymandering metrics for a set of plans
#'
#' \code{redist.metrics} is used to compute different gerrymandering metrics for a
#' set of maps.
#'
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param measure A vector with a string for each measure desired from list "DSeats", "DVS", "EffGap",
#' "EffGapEqPop", "TauGap", "MeanMedian", "Bias", "BiasV", "Declination",
#' "Responsiveness", "LopsidedWins", "RankedMarginal", and "SmoothedSeat". Use "all" to get all metrics.
#' "DSeats" and "DVS" are always computed, so it is recommended to always return those values.
#' @param rvote A numeric vector with the Republican vote for each precinct.
#' @param dvote A numeric vector with the Democratic vote for each precinct.
#' @param draw A numeric to specify draw number. Defaults to 1 if only one map provided
#' and the column number if multiple maps given. Can also take a factor input, which will become the
#' draw column in the output if its length matches the number of entries in plans. If the `plans` input
#' is a `redist_plans` object, it extracts the `draw` identifier.
#' @param tau A non-negative number for calculating Tau Gap. Only used with option "TauGap". Defaults to 1.
#' @param biasV A value between 0 and 1 to compute bias at. Only used with option "BiasV". Defaults to 0.5.
#' @param respV A value between 0 and 1 to compute responsiveness at. Only used with option "Responsiveness". Defaults to 0.5.
#' @param bandwidth A value between 0 and 1 for computing responsiveness. Only used with option "Responsiveness." Defaults to 0.01.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#'
#'
#' @details This function computes specified compactness scores for a map.  If
#' there is more than one precinct specified for a map, it aggregates to the district level
#' and computes one score.
#'
#' - DSeats is computed as the expected number of Democratic seats with no change in votes.
#' - DVS is the Democratic Vote Share, which is the two party vote share with Democratic votes as the numerator.
#' - EffGap is the Efficiency Gap, calculated with votes directly.
#' - EffGapEqPop is the Efficiency Gap under an Equal Population assumption, calculated with the DVS.
#' - TauGap is the Tau Gap, computed with the Equal Population assumption.
#' - MeanMedian is the Mean Median difference.
#' - Bias is the Partisan Bias computed at 0.5.
#' - BiasV is the Partisan Bias computed at value V.
#' - Declination is the value of declination at 0.5.
#' - Responsiveness is the responsiveness at the user-supplied value with the user-supplied bandwidth.
#' - LopsidedWins computed the Lopsided Outcomes value, but does not produce a test statistic.
#' - RankedMarginal computes the Ranked Marginal Deviation (0-1, smaller is better). This is also known
#' as the "Gerrymandering Index" and is sometimes presented as this value divided by 10000.
#' - SmoothedSeat computes the Smoothed Seat Count Deviation (0-1, smaller is R Bias, bigger is D Bias).
#'
#' @return A tibble with  a column for each specified measure and
#' a column that specifies the map number.
#'
#' @importFrom dplyr select %>% tibble
#'
#' @examples
#' data(fl25)
#' data(fl25_enum)
#' plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#' redist.metrics(plans_05, measure = 'all', rvote = fl25$mccain, dvote = fl25$obama)
#'
#' @references
#' Jonathan N. Katz, Gary King, and Elizabeth Rosenblatt. 2020.
#' Theoretical Foundations and Empirical Evaluations of Partisan Fairness in District-Based Democracies.
#' American Political Science Review, 114, 1, Pp. 164-178.
#'
#' Gregory S. Warrington. 2018. "Quantifying Gerrymandering Using the Vote Distribution."
#' Election Law Journal: Rules, Politics, and Policy. Pp. 39-57.http://doi.org/10.1089/elj.2017.0447
#'
#' Samuel S.-H. Wang. 2016. "Three Tests for Practical Evaluation of Partisan Gerrymandering."
#' Stanford Law Review, 68, Pp. 1263 - 1321.
#'
#' Gregory Herschlag, Han Sung Kang, Justin Luo, Christy Vaughn Graves, Sachet Bangia,
#' Robert Ravier & Jonathan C. Mattingly (2020) Quantifying Gerrymandering in North Carolina,
#' Statistics and Public Policy, 7:1, 30-38, DOI: 10.1080/2330443X.2020.1796400
#'
#' @md
#' @concept analyze
#' @export
redist.metrics <- function(plans, measure = "DSeats", rvote, dvote,
                           tau = 1, biasV = 0.5, respV = 0.5, bandwidth = 0.01,
                           draw = 1, ncores = 1){

  # All measures available:
  all_measures <- c("DSeats", "DVS", "EffGap", "EffGapEqPop", "TauGap",
                    "MeanMedian", "Bias", "BiasV", "Declination", "Responsiveness",
                    "LopsidedWins", "RankedMarginal", "SmoothedSeat")

  # Check Inputs
  if("all" %in% measure){
    measure <-  all_measures
  }
  match.arg(arg = measure,several.ok = TRUE, choices = all_measures)

  if(any(class(plans) %in% 'redist')){
    plans <- plans$plans
  }

  if (inherits(plans, 'redist_plans')) {
    draw <- plans$draw
    plans <- get_plans_matrix(plans)
  }

  if(!is.numeric(plans)){
    stop('Please provide "plans" as a numeric vector or matrix.')
  }
  if(!is.matrix(plans)){
    plans <- as.matrix(plans)
  }
  if(any(is.na(plans))){
    stop('NA value in argument to plans.')
  }

  if(any(is.na(rvote))){
    stop('NA value in argument to rvote.')
  }
  if(any(is.na(dvote))){
    stop('NA value in argument to dvote.')
  }
  if(!is.numeric(rvote)){
    stop('Please provide rvote as a numeric or integer vector.')
  }
  if(!is.numeric(dvote)){
    stop('Please provide rvote as a numeric or integer vector.')
  }

  rvote <- as.integer(rvote)
  dvote <- as.integer(dvote)
  if(length(rvote) != nrow(plans)){
    stop('rvote length and plans row dimension are not equal.')
  }
  if(length(dvote) != nrow(plans)){
    stop('dvote length and plans row dimension are not equal.')
  }

  if(!is.numeric(draw) & !is.factor(draw)){
    stop('Please provide "draw" as a numeric.')
  }

  if(!is.numeric(ncores)){
    stop('Please provide "ncores" as a numeric.')
  }


  # Precompute a few useful variables
  nd <- length(unique(plans[, 1]))
  totvote <- sum(rvote) + sum(dvote)
  nmap <- ncol(plans)
  dists <- sort(unique(plans[, 1]))

  # Aggregate to Precinct and get baseline DVS + Seats - compute here to avoid multiple computation
  rcounts <- agg_p2d(vote = rvote, dm = plans, nd = nd)
  dcounts <- agg_p2d(vote = dvote, dm = plans, nd = nd)
  dseat_vec <- dseats(dm = plans, rcounts = rcounts, dcounts = dcounts, nd = nd)
  dvs <- DVS(dcounts = dcounts, rcounts = rcounts)



  # Create return tibble:
  if (!(is.factor(draw) && length(draw) == nd * nmap)) {
    if (nmap != 1) {
      draw <- rep(draw + (1:ncol(plans)) - 1, each = nd)
    } else {
      draw <- rep(draw, nd)
    }
  }

  metrics <- tibble(district = rep(x = dists, nmap),
                 DSeats = rep(NA_real_, nd*nmap),
                 DVS = rep(NA_real_, nd*nmap),
                 EffGap = rep(NA_real_, nd*nmap),
                 EffGapEqPop = rep(NA_real_, nd*nmap),
                 TauGap = rep(NA_real_, nd*nmap),
                 MeanMedian = rep(NA_real_, nd*nmap),
                 Bias = rep(NA_real_, nd*nmap),
                 BiasV = rep(NA_real_, nd*nmap),
                 Declination = rep(NA_real_, nd*nmap),
                 Responsiveness = rep(NA_real_, nd*nmap),
                 LopsidedWins = rep(NA_real_, nd*nmap),
                 RankedMarginal = rep(NA_real_, nd*nmap),
                 SmoothedSeat = rep(NA_real_, nd*nmap),
                 draw = draw) %>%
    dplyr::select(all_of(c("district", measure)), draw)

  # Compute Metrics if desired:
  if("DSeats" %in% measure){
    metrics[['DSeats']] <- rep(dseat_vec, each = nd)
  }
  if('DVS' %in% measure){
    metrics[['DVS']] <- c(dvs)
  }
  if("EffGap" %in% measure){
    eg <- effgap(dcounts = dcounts, rcounts = rcounts, totvote = totvote)
    metrics[['EffGap']] <-  rep(eg, each = nd)
  }
  if("EffGapEqPop" %in% measure){
    egep <- effgapEP(dvs = dvs, dseat_vec = dseat_vec, nd = nd)
    metrics[["EffGapEqPop" ]] <- rep(egep, each = nd)
  }
  if("TauGap" %in% measure){
    tg <- taugap(tau = tau, dvs = dvs, dseat_vec = dseat_vec, nd = nd)
    metrics[["TauGap"]] <- rep(tg, each = nd)
  }
  if("MeanMedian" %in% measure){
    mm <- meanmedian(dvs = dvs)
    metrics[["MeanMedian"]]<- rep(mm, each = nd)
  }
  if("Bias" %in% measure){
    b <- bias(dvs = dvs, nd = nd)
    metrics[["Bias"]] <- rep(b, each = nd)
  }
  if("BiasV" %in% measure){
    bv <- biasatv(dvs = dvs, v = biasV, nd = nd)
    metrics[["BiasV"]] <- rep(bv, each = nd)
  }
  if("Declination" %in% measure){
     dec <- declination(dvs = dvs, dseat_vec = dseat_vec, nd = nd)
    metrics[["Declination"]] <- rep(dec, each = nd)
  }
  if("Responsiveness" %in% measure){
    resp <- responsiveness(dvs = dvs, v = respV, nd = nd, bandwidth = bandwidth)
    metrics[["Responsiveness"]] <- rep(resp, each = nd)
  }
  if("LopsidedWins" %in% measure){
    lw <- lopsidedwins(dvs = dvs, dseat_vec = dseat_vec, nd = nd)
    metrics[["LopsidedWins"]] <- rep(lw, each = nd)
  }
  if('RankedMarginal' %in% measure){
    RMDev <- RankedMarginalDev(dvs = dvs)
    metrics[["RankedMarginal"]] <- rep(RMDev, each = nd)
  }
  if('SmoothedSeat' %in% measure){
    SSCD <- smoothseat(dvs = dvs, nd = nd)
    metrics[['SmoothedSeat']] <- rep(SSCD, each = nd)
  }

  # Return computed results
  return(metrics)
}
