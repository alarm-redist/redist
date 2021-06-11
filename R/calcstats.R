#' Segregation index calculation for MCMC redistricting.
#'
#' \code{redist.segcalc} calculates the dissimilarity index of segregation (see
#' Massey \& Denton 1987 for more details) for a specified subgroup under any
#' redistricting plan.
#'
#' @param plans A matrix of congressional district assignments or a
#' redist object.
#' @param group_pop A vector of populations for some subgroup of interest.
#' @param total_pop A vector containing the populations of each geographic unit.
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
#' \donttest{
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#'
#' ## 25 precinct, three districts - no pop constraint ##
#' alg_253 <- redist.flip(adj = fl25_adj, total_pop = fl25$pop,
#'                        init_plan = init_plan, nsims = 10000)
#'
#' ## Get Republican Dissimilarity Index from simulations
#' rep_dmi_253 <- redist.segcalc(alg_253, fl25$mccain, fl25$pop)
#' }
#' @concept analyze
#' @export
redist.segcalc <- function(plans, group_pop, total_pop){
    ## Warnings
    if(missing(plans) | !any(class(plans) %in% c("data.frame", "matrix", "redist"))){
        stop("Please provide either a redist object or a proper matrix of congessional districts")
    }
    if(missing(group_pop)){
        stop("Please provide a vector of sub-group populations to ",
             "calculate the segregation index")
    }
    if(missing(total_pop)){
        stop("Please provide a vector of populations for each geographic unit")
    }

    ## If redist object, get the partitions entry
    if(all(class(plans) == 'redist')){
        plans <- plans$plans
    }

    if(!((nrow(plans) == length(group_pop)) &
             (length(group_pop) == length(total_pop)) &
                 (length(total_pop) == nrow(plans)))){
        stop("Please make sure there is a population entry for each geographic unit")
    }

    ## Calculate dissimilarity index
    seg.out <- segregationcalc(plans,
                               group_pop,
                               total_pop)

    ## Return
    return(seg.out)

}
