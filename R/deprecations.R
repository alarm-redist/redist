#' Segregation index calculation for MCMC redistricting.
#'
#' \code{redist.segcalc} calculates the dissimilarity index of segregation (see
#' Massey & Denton 1987 for more details) for a specified subgroup under any
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
#' fl25$init_plan <- init_plan
#'
#' ## 25 precinct, three districts - no pop constraint ##
#' fl_map <- redist_map(fl25, existing_plan = 'init_plan', adj = fl25_adj)
#' alg_253 <- redist_flip(fl_map, nsims = 10000)
#'
#'
#' ## Get Republican Dissimilarity Index from simulations
#' # old: rep_dmi_253 <- redist.segcalc(alg_253, fl25$mccain, fl25$pop)
#' rep_dmi_253 <- seg_dissim(alg_253, fl25, mccain, pop)  |>
#'     redistmetrics::by_plan(ndists = 3)
#' }
#' @concept analyze
#' @export
redist.segcalc <- function(plans, group_pop, total_pop) {
    .Deprecated("seg_dissim")

    ## If redist object, get the partitions entry
    if (all(class(plans) == "redist")) {
        plans <- plans$plans
    }

    if (!((nrow(plans) == length(group_pop)) &
        (length(group_pop) == length(total_pop)) &
        (length(total_pop) == nrow(plans)))) {
        cli_abort("Please make sure there is a population entry for each geographic unit")
    }

    nd <- dplyr::n_distinct((plans[, 1]))
    out <- redistmetrics::seg_dissim(plans,
        shp = data.frame(), group_pop = group_pop,
        total_pop = total_pop
    )
    out[seq(1, length(out), by = nd)]
}
#' Compute Competitiveness
#'
#' Currently only implements the competitiveness function in equation (5)
#' of Cho & Liu 2016.
#'
#'
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param rvote A numeric vector with the Republican vote for each precinct.
#' @param dvote A numeric vector with the Democratic vote for each precinct.
#' @param alpha A numeric value for the alpha parameter for the talisman metric
#' @param beta A numeric value for the beta parameter for the talisman metric
#' @return Numeric vector with competitiveness scores
#' @export
#'
#' @concept analyze
#' @examples
#' data(fl25)
#' data(fl25_enum)
#'
#' plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#' # old: comp <- redist.competitiveness(plans_05, fl25$mccain, fl25$obama)
#' comp <- compet_talisman(plans_05, fl25, mccain, obama)
#'
redist.competitiveness <- function(plans, rvote, dvote, alpha = 1, beta = 1) {
    .Deprecated("compet_talisman")
    nd <- length(unique(plans[, 1]))
    redistmetrics::compet_talisman(plans = plans, shp = data.frame(),
        rvote = rvote, dvote = dvote,
        alpha = alpha, beta = beta)[seq(1, nd*ncol(plans), by = nd)]
}
##############################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/01/20
## Date Modified: 2022/01/13
## Purpose: R function to compute compactness
##############################################


#' Calculate compactness measures for a set of plans
#'
#' \code{redist.compactness} is used to compute different compactness statistics for a
#' shapefile. It currently computes the Polsby-Popper, Schwartzberg score, Length-Width Ratio,
#' Convex Hull score, Reock score, Boyce Clark Index, Fryer Holden score, Edges Removed number,
#' and the log of the Spanning Trees.
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required unless "EdgesRemoved"
#' and "logSpanningTree" with adjacency provided.
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param measure A vector with a string for each measure desired. "PolsbyPopper",
#' "Schwartzberg", "LengthWidth", "ConvexHull", "Reock", "BoyceClark", "FryerHolden",
#' "EdgesRemoved", "FracKept", and "logSpanningTree" are implemented. Defaults to "PolsbyPopper". Use "all" to
#' return all implemented measures.
#' @param total_pop A numeric vector with the population for every observation. Is
#' only necessary when "FryerHolden" is used for measure. Defaults to NULL.
#' @param adj A zero-indexed adjacency list. Only used for "PolsbyPopper",
#' EdgesRemoved" and "logSpanningTree". Created with \code{redist.adjacency} if not
#' supplied and needed. Default is NULL.
#' @param draw A numeric to specify draw number. Defaults to 1 if only one map provided
#' and the column number if multiple maps given. Can also take a factor input, which will become the
#' draw column in the output if its length matches the number of entries in plans. If the `plans` input
#' is a `redist_plans` object, it extracts the `draw` identifier.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' @param counties A numeric vector from 1:ncounties corresponding to counties. Required for "logSpanningTree".
#' @param planarize a number, indicating the CRS to project the shapefile to if
#' it is latitude-longitude based. Set to FALSE to avoid planarizing.
#' @param ppRcpp Boolean, whether to run Polsby Popper and Schwartzberg using Rcpp.
#' It has a higher upfront cost, but quickly becomes faster.
#' Becomes TRUE if ncol(district_membership > 8) and not manually set.
#' @param perim_path it checks for an Rds, if no rds exists at the path,
#' it creates an rds with borders and saves it.
#' This can be created in advance with [redistmetrics::prep_perims()].
#' @param perim_df A dataframe output from [redistmetrics::prep_perims()].
#'
#' @details This function computes specified compactness scores for a map.  If
#' there is more than one shape specified for a single district, it combines
#' them, if necessary, and computes one score for each district.
#'
#'
#' Polsby-Popper is computed as \deqn{\frac{4*\pi*A(d)}{P(d)^2}} where A is the area
#' function, the district is d, and P is the perimeter function. All  values are between
#' 0 and 1, where larger values are more compact.
#'
#' Schwartzberg is computed as \deqn{\frac{P(d)}{2*\pi*\sqrt{\frac{A(d)}{\pi}}}}
#' where A is the area function, the district is d, and P is the perimeter function.
#' All  values are between 0 and 1, where larger values are more compact.
#'
#' The Length Width ratio is computed as \deqn{\frac{length}{width}} where length
#' is the shorter of the maximum x distance and the maximum y distance. Width is
#' the longer of the two values. All  values are between 0 and 1, where larger
#' values are more compact.
#'
#' The Convex Hull score is computed as \deqn{\frac{A(d)}{A(CVH)}} where A is the area
#' function, d is the district, and CVH is the convex hull of the district. All
#' values are between 0 and 1, where larger values are more compact.
#'
#' The Reock score is computed as \deqn{\frac{A(d)}{A(MBC)}} where A is the area
#' function, d is the district, and MBC is the minimum bounding circle of the
#' district. All values are between 0 and 1, where larger values are more compact.
#'
#' The Boyce Clark Index is computed as \deqn{1 - \sum_{1}^{16}\{\frac{|\frac{r_i}{\sum_ir_i}*100-6.25 |\}}{200}}.
#' The \eqn{r_i} are the distances of the 16 radii computed from the geometric
#' centroid of the shape to the most outward point of the shape that intersects
#' the radii, if the centroid is contained within the shape.  If the centroid
#' lies outside of the shape, a point on the surface is used, which will naturally
#' incur a penalty to the score. All  values are between 0 and 1,
#' where larger values are more compact.
#'
#' The Fryer Holden score for each district is computed with \deqn{Pop\odot D(precinct)^2},
#' where \eqn{Pop} is the population product matrix.  Each element is the
#' product of the i-th and j-th precinct's populations.  D represents the distance,
#' where the matrix is the distance between each precinct.  To fully compute this
#' index, for any map, the sum of these values should be used as the numerator.
#' The denominator can be calculated from the full enumeration of districts as the
#' smallest calculated numerator. This produces very large numbers, where smaller
#' values are more compact.
#'
#' The log spanning tree measure is the logarithm of the product of the
#' number of spanning trees which can be drawn on each district.
#'
#' The edges removed measure is number of edges removed from the underlying adjacency graph.
#' A smaller number of edges removed is more compact.
#'
#' The fraction kept measure is the fraction of edges that were not removed from the
#' underlying adjacency graph. This takes values 0 - 1, where 1 is more compact.
#'
#' @return A tibble with a column that specifies the district, a column for
#' each specified measure, and a column that specifies the map number.
#'
#' @references Boyce, R., & Clark, W. 1964. The Concept of Shape in Geography.
#' Geographical Review, 54(4), 561-572.
#'
#' Cox, E. 1927. A Method of Assigning Numerical and Percentage Values to the
#' Degree of Roundness of Sand Grains. Journal of Paleontology, 1(3), 179-183.
#'
#' Fryer R, Holden R. 2011. Measuring the Compactness of Political Districting Plans.
#' Journal of Law and Economics.
#'
#' Harris, Curtis C. 1964. “A scientific method of districting”.
#' Behavioral Science 3(9), 219–225.
#'
#' Maceachren, A. 1985. Compactness of Geographic Shape: Comparison and
#' Evaluation of Measures. Geografiska Annaler. Series B, Human Geography, 67(1),
#' 53-67.
#'
#' Polsby, Daniel D., and Robert D. Popper. 1991. “The Third Criterion:
#' Compactness as a procedural safeguard against partisan gerrymandering.”
#' Yale Law & Policy Review 9 (2): 301–353.
#'
#' Reock, E. 1961. A Note: Measuring Compactness as a Requirement of Legislative
#' Apportionment. Midwest Journal of Political Science, 5(1), 70-74.
#'
#' Schwartzberg, Joseph E. 1966. Reapportionment, Gerrymanders, and the Notion
#' of Compactness. Minnesota Law Review. 1701.
#'
#' @importFrom dplyr tibble %>%
#' @importFrom sf st_cast st_bbox st_centroid st_within st_point_on_surface st_coordinates
#' @importFrom sf st_linestring st_intersection st_area st_crs st_is_longlat st_length
#' @importFrom sf st_convex_hull st_crs<- st_geometry st_distance st_union st_touches st_is_valid
#' @importFrom sf st_is_longlat
#' @importFrom dplyr select all_of arrange bind_rows rename summarize
#' @importFrom stats dist
#'
#' @concept analyze
#' @examples
#' data(fl25)
#' data(fl25_enum)
#'
#' plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#'
#' # old redist.compactness(
#' #     shp = fl25, plans = plans_05[, 1:3],
#' #     measure = c("PolsbyPopper", "EdgesRemoved")
#' # )
#' comp_polsby(plans_05[, 1:3], fl25)
#' comp_edges_rem(plans_05[, 1:3], fl25, fl25$adj)
#' @export redist.compactness
redist.compactness <- function(shp = NULL,
                               plans,
                               measure = c("PolsbyPopper"),
                               total_pop = NULL, adj = NULL, draw = 1,
                               ncores = 1, counties = NULL, planarize = 3857,
                               ppRcpp, perim_path, perim_df) {
    repl = c(
        PolsbyPopper = "comp_polsby",
        Schwartzberg = "comp_schwartz",
        LengthWidth = "comp_lw",
        ConvexHull = "comp_ch",
        Reock = "comp_reock",
        BoyceClark = "comp_bc",
        FryerHolden = "comp_fh",
        EdgesRemoved = "comp_edges_rem",
        FracKept = "comp_frac_kept",
        logSpanningTree = "comp_log_st"
    )[measure]
    .Deprecated(repl)

    # Check Inputs
    if (is.null(shp) & is.null(adj)) {
        cli_abort("Please provide a {.arg shp} or {.arg adj} argument.")
    }

    if (!is.null(shp)) {
        if ("SpatialPolygonsDataFrame" %in% class(shp)) {
            shp <- sf::st_as_sf(shp)
        } else if (!inherits(shp, "sf")) {
            cli_abort("Please provide {.arg shp} as a SpatialPolygonsDataFrame or sf object.")
        }


        if (isTRUE(st_is_longlat(st_geometry(shp)))) {
            if (!is.null(st_crs(shp)) & !is.null(planarize) && !isFALSE(planarize)) {
                shp <- sf::st_transform(shp, planarize)
            }
        }
    }

    if (inherits(shp, "redist_map") & missing(adj)) {
        adj <- get_adj(shp)
    }

    if (inherits(plans, "redist_plans")) {
        draw <- plans$draw
        plans <- get_plans_matrix(plans)
    }

    if (!is.numeric(plans)) {
        cli_abort("Please provide {.arg plans} as a numeric vector, matrix, or {.cls redist_plans}.")
    }


    possible_measures <- c(
        "PolsbyPopper", "Schwartzberg", "LengthWidth", "ConvexHull",
        "Reock", "BoyceClark", "FryerHolden", "EdgesRemoved", "FracKept",
        "logSpanningTree"
    )
    if ("all" %in% measure) {
        measure <- possible_measures
    }

    match.arg(
        arg = measure, several.ok = TRUE,
        choices = possible_measures
    )

    if ("FryerHolden" %in% measure & is.null(total_pop)) {
        cli_abort("Please provide a {.arg total_pop} argument when FryerHolden is specified.")
    }

    if ("FryerHolden" %in% measure) {
        if (!any(class(total_pop) %in% c("numeric", "integer"))) {
            cli_abort("Please provide {.arg total_pop} as a numeric or integer.")
        }
    }

    if (!is.numeric(draw) & !is.factor(draw)) {
        cli_abort("Please provide {.arg draw} as a numeric or factor.")
    }

    if (!is.numeric(ncores)) {
        cli_abort('Please provide "ncores" as a numeric.')
    }
    if (("logSpanningTree" %in% measure) & is.null(counties)) {
        cli_abort("Please provide {.arg counties}.")
    }


    # Compute compactness scores
    dists <- sort(unique(c(plans)))
    nd <- length(dists)

    if (!is.matrix(plans)) {
        plans <- as.matrix(plans)
    }
    V <- nrow(plans)

    if (missing(ppRcpp)) {
        if (ncol(plans) > 8 || !missing(perim_path) || !missing(perim_df)) {
            ppRcpp <- TRUE
        } else {
            ppRcpp <- FALSE
        }
    }

    nmap <- ncol(plans)
    if (!(is.factor(draw) && length(draw) == nd*nmap)) {
        if (nmap != 1) {
            draw <- rep(draw + (1:ncol(plans)) - 1, each = nd)
        } else {
            draw <- rep(draw, nd)
        }
    }

    # Initialize object
    comp <- tibble(
        district = rep(x = dists, nmap),
        PolsbyPopper = rep(NA_real_, nd*nmap),
        Schwartzberg = rep(NA_real_, nd*nmap),
        LengthWidth = rep(NA_real_, nd*nmap),
        ConvexHull = rep(NA_real_, nd*nmap),
        Reock = rep(NA_real_, nd*nmap),
        BoyceClark = rep(NA_real_, nd*nmap),
        FryerHolden = rep(NA_real_, nd*nmap),
        EdgesRemoved = rep(NA_real_, nd*nmap),
        FracKept = rep(NA_real_, nd*nmap),
        logSpanningTree = rep(NA_real_, nd*nmap),
        draw = draw
    ) %>%
        dplyr::select(all_of(c("district", measure)), all_of(measure), draw)

    # Compute Specified Scores for provided districts
    if ("PolsbyPopper" %in% measure) {
        comp$PolsbyPopper <- redistmetrics::comp_polsby(plans,
            shp = shp,
            use_Rcpp = ppRcpp,
            perim_path = perim_path,
            perim_df = perim_df,
            epsg = planarize,
            ncores = ncores
        )
    }

    if ("Schwartzberg" %in% measure) {
        comp$Schwartzberg <- redistmetrics::comp_schwartz(
            plans = plans,
            shp = shp,
            use_Rcpp = ppRcpp,
            perim_path = perim_path,
            perim_df = perim_df,
            epsg = planarize,
            ncores = ncores
        )
    }

    if ("ConvexHull" %in% measure) {
        comp$ConvexHull <- redistmetrics::comp_ch(
            plans = plans,
            shp = shp,
            epsg = planarize,
            ncores = ncores
        )
    }

    if ("Reock" %in% measure) {
        comp$Reock <- redistmetrics::comp_reock(
            plans = plans,
            shp = shp,
            epsg = planarize,
            ncores = ncores
        )
    }

    if ("BoyceClark" %in% measure) {
        comp$BoyceClark <- redistmetrics::comp_bc(
            plans = plans,
            shp = shp,
            epsg = planarize,
            ncores = ncores
        )
    }

    if ("FryerHolden" %in% measure) {
        comp$FryerHolden <- redistmetrics::comp_fh(
            plans = plans,
            shp = shp,
            total_pop = total_pop,
            epsg = planarize,
            ncores = ncores
        )
    }

    if ("LengthWidth" %in% measure) {
        comp$LengthWidth <- redistmetrics::comp_lw(
            plans = plans,
            shp = shp,
            epsg = planarize,
            ncores = ncores
        )
    }


    if (any(measure %in% c("EdgesRemoved", "logSpanningTree", "FracKept")) & is.null(adj)) {
        adj <- redist.adjacency(shp)
    }

    if ("logSpanningTree" %in% measure) {
        comp$logSpanningTree <- redistmetrics::comp_log_st(plans, shp,
            counties = counties,
            adj = adj
        )
    }

    if ("EdgesRemoved" %in% measure) {
        comp$EdgesRemoved <- redistmetrics::comp_edges_rem(plans, shp, adj = adj)
    }

    if ("FracKept" %in% measure) {
        comp$FracKept <- redistmetrics::comp_frac_kept(plans, shp, adj = adj)
    }

    # Return results
    comp
}


#' Calculate Group Proportion by District
#'
#' \code{redist.group.percent} computes the proportion that a group makes up in
#' each district across a matrix of maps.
#'
#' @param plans A matrix with one row
#' for each precinct and one column for each map. Required.
#' @param group_pop A numeric vector with the population of the group for every precinct.
#' @param total_pop A numeric vector with the population for every precinct.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#'
#' @return matrix with percent for each district
#'
#' @export
#' @concept analyze
#'
#' @examples
#' data(fl25)
#' data(fl25_enum)
#'
#' cd <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#' fl25_map = redist_map(fl25, ndists=3, pop_tol=0.1)
#' fl25_plans = redist_plans(cd, fl25_map, algorithm="enumpart")
#'
#' group_frac(fl25_map, BlackPop, TotPop, fl25_plans)
redist.group.percent <- function(plans, group_pop, total_pop, ncores = 1) {
    .Deprecated("group_frac()")
    if (!is.numeric(group_pop) || !is.numeric(total_pop))
        cli_abort("{.arg group_pop} and {.arg total_pop} must be numeric vectors.")
    if (!is.matrix(plans))
        cli_abort("{.arg plans} must be a matrix.")

    if (!is.matrix(plans)) {
        plans <- as.matrix(plans)
    }

    if (length(total_pop) != nrow(plans))
        cli_abort("{.arg plans} and {.arg total_pop} must have the same number of precincts.")
    if (length(group_pop) != nrow(plans))
        cli_abort("{.arg plans} and {.arg group_pop} must have the same number of precincts.")

    ndists <- max(plans[, 1])
    if (ndists ==  length(unique(plans[, 1])) - 1) {
        plans <- plans + 1
        ndists <- ndists + 1
    }
    group_pct(plans, group_pop, total_pop, ndists)
}




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
#' # old: redist.metrics(plans_05, measure = "DSeats", rvote = fl25$mccain, dvote = fl25$obama)
#' part_dseats(plans_05, fl25, mccain, obama)
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
                           draw = 1) {
    repl = c(
        DSeats = "part_dseats",
        DVS = "part_dvs",
        EffGap = "part_egap",
        EffGapEqPop = "part_egap_ep",
        TauGap = "part_tau_gap",
        MeanMedian = "part_mean_median",
        Bias = "part_bias",
        BiasV = "part_bias",
        Declination = "part_decl",
        Responsiveness = "part_resp",
        LopsidedWins = "part_lop_wins",
        RankedMarginal = "part_rmd",
        SmoothedSeat = "part_sscd"
    )[measure]
    .Deprecated(repl)

    # All measures available:
    all_measures <- c("DSeats", "DVS", "EffGap", "EffGapEqPop", "TauGap",
        "MeanMedian", "Bias", "BiasV", "Declination", "Responsiveness",
        "LopsidedWins", "RankedMarginal", "SmoothedSeat")

    # Check Inputs
    if ("all" %in% measure) {
        measure <-  all_measures
    }
    match.arg(arg = measure, several.ok = TRUE, choices = all_measures)

    if (inherits(plans, "redist_plans")) {
        draw <- plans$draw
        plans <- get_plans_matrix(plans)
    }

    if (!is.numeric(plans)) {
        cli_abort("Please provide {.arg plans} as a numeric vector, matrix, or {.cls redist_plans}.")
    }
    if (!is.matrix(plans)) {
        plans <- as.matrix(plans)
    }
    if (any(is.na(plans))) {
        cli_abort("{.val NA} in argument to {.arg plans}.")
    }

    if (any(is.na(rvote))) {
        cli_abort("{.val NA} value in argument to {.arg rvote}.")
    }
    if (any(is.na(dvote))) {
        cli_abort("{.val NA} value in argument to {.arg dvote}.")
    }
    if (!is.numeric(rvote)) {
        cli_abort("Please provide {.arg rvote} as a numeric or integer vector.")
    }
    if (!is.numeric(dvote)) {
        cli_abort("Please provide {.arg dvote} as a numeric or integer vector.")
    }

    rvote <- as.integer(rvote)
    dvote <- as.integer(dvote)
    if (length(rvote) != nrow(plans)) {
        cli_abort("{.arg rvote} length and {.arg plans} row dimension are not equal.")
    }
    if (length(dvote) != nrow(plans)) {
        cli_abort("{.arg dvote} length and {.arg plans} row dimension are not equal.")
    }

    if (!is.numeric(draw) & !is.factor(draw)) {
        cli_abort('Please provide "draw" as a numeric or factor.')
    }

    # Precompute a few useful variables
    nd <- length(unique(plans[, 1]))
    nmap <- ncol(plans)
    dists <- sort(unique(plans[, 1]))

    # Create return tibble:
    if (!(is.factor(draw) && length(draw) == nd*nmap)) {
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
    if ("DSeats" %in% measure) {
        metrics[["DSeats"]] <- redistmetrics::part_dseats(plans = plans, shp = data.frame(),
            rvote = rvote, dvote = dvote)
    }
    if ("DVS" %in% measure) {
        metrics[["DVS"]] <- redistmetrics::part_dvs(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("EffGap" %in% measure) {
        metrics[["EffGap"]] <- redistmetrics::part_egap(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("EffGapEqPop" %in% measure) {
        metrics[["EffGapEqPop"]] <- redistmetrics::part_egap_ep(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("TauGap" %in% measure) {
        metrics[["TauGap"]] <- redistmetrics::part_tau_gap(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("MeanMedian" %in% measure) {
        metrics[["MeanMedian"]] <- redistmetrics::part_mean_median(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("Bias" %in% measure) {
        metrics[["Bias"]] <- redistmetrics::part_bias(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote, v = 0.5)
    }
    if ("BiasV" %in% measure) {
        metrics[["BiasV"]] <- redistmetrics::part_bias(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote, v = biasV)
    }
    if ("Declination" %in% measure) {
        metrics[["Declination"]] <- redistmetrics::part_decl(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("Responsiveness" %in% measure) {
        metrics[["Responsiveness"]] <- redistmetrics::part_resp(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote,
            v = respV, bandwidth = bandwidth)
    }
    if ("LopsidedWins" %in% measure) {
        metrics[["LopsidedWins"]] <- redistmetrics::part_lop_wins(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("RankedMarginal" %in% measure) {
        metrics[["RankedMarginal"]] <- redistmetrics::part_rmd(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }
    if ("SmoothedSeat" %in% measure) {
        metrics[["SmoothedSeat"]] <- redistmetrics::part_sscd(plans = plans, shp = data.frame(),
            dvote = dvote, rvote = rvote)
    }

    # Return computed results
    metrics
}

#' Count County Splits
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param counties A vector of county names or county ids.
#'
#' @return integer vector with one number for each map
#'
#' @concept analyze
#' @export
redist.splits <- function(plans, counties) {
    .Deprecated("splits_admin")
    if (missing(plans)) {
        cli_abort("Please provide an argument to {.arg plans}.")
    }
    if (inherits(plans, "redist_plans")) {
        plans <- get_plans_matrix(plans)
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol = 1)
    }
    if (any(class(plans) %in% c("numeric", "matrix"))) {
        plans <- redist.reorder(plans)
    }

    if (missing(counties)) {
        cli_abort("Please provide an argument to {.arg counties}.")
    }


    redistmetrics::splits_admin(plans = plans, shp = data.frame(),
        admin = counties)
}

#' Counts the Number of Counties within a District
#'
#' Counts the total number of counties that are found within a district.
#' This does not subtract out the number of counties that are found completely
#' within a district.
#'
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param counties A vector of county names or county ids.
#'
#' @return integer matrix where each district is a
#'
#' @concept analyze
#' @export
#'
#' @examples
#' data(iowa)
#' ia <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop, pop_tol = 0.01)
#' plans <- redist_smc(ia, 50, silent = TRUE)
#' #old redist.district.splits(plans, ia$region)
#' splits_count(plans, ia, region)
redist.district.splits <- function(plans, counties) {
    .Deprecated("splits_count")
    if (missing(plans)) {
        stop("Please provide an argument to plans.")
    }
    if (inherits(plans, "redist_plans")) {
        plans <- get_plans_matrix(plans)
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol = 1)
    }
    if (!any(class(plans) %in% c("numeric", "matrix"))) {
        stop('Please provide "plans" as a matrix.')
    }

    if (missing(counties)) {
        stop("Please provide an argument to counties.")
    }
    if (class(counties) %in% c("character", "numeric", "integer")) {
        county_id <- redist.county.id(counties)
    } else {
        stop('Please provide "counties" as a character, numeric, or integer vector.')
    }


    dist_cty_splits(plans - 1, community = county_id - 1, length(unique(plans[, 1])))
}


#' Counts the Number of Counties Split Between 3 or More Districts
#'
#' Counts the total number of counties that are split across more than 2 districts.
#'
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param counties A vector of county names or county ids.
#'
#' @return integer matrix where each district is a
#'
#' @concept analyze
#' @export
#'
#' @examples
#' data(iowa)
#' ia <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop, pop_tol = 0.01)
#' plans <- redist_smc(ia, 50, silent = TRUE)
#' #old redist.multisplits(plans, ia$region)
#' splits_multi(plans, ia, region)
redist.multisplits <- function(plans, counties) {
    .Deprecated("splits_multi")
    if (missing(plans)) {
        cli_abort("Please provide an argument to {.arg plans}.")
    }
    if (inherits(plans, "redist_plans")) {
        plans <- get_plans_matrix(plans)
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol = 1)
    }
    if (any(class(plans) %in% c("numeric", "matrix"))) {
        plans <- redist.reorder(plans)
    }

    if (missing(counties)) {
        cli_abort("Please provide an argument to {.arg counties}.")
    }

    redistmetrics::splits_multi(plans = plans, shp = data.frame(),
        admin = counties)
}


#' Counts the Number of Municipalities Split Between Districts
#'
#' Counts the total number of municpalities that are split.
#' Municipalities in this interpretation do not need to cover the entire state, which
#' differs from counties.
#'
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param munis A vector of municipality names or ids.
#'
#' @return integer vector of length ndist by ncol(plans)
#'
#' @concept analyze
#' @export
#'
#' @examples
#' data(iowa)
#' ia <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop, pop_tol = 0.01)
#' plans <- redist_smc(ia, 50, silent = TRUE)
#' ia$region[1:10] <- NA
#' #old redist.muni.splits(plans, ia$region)
#' splits_sub_admin(plans, ia, region)
redist.muni.splits <- function(plans, munis) {
    .Deprecated("splits_sub_admin")
    redistmetrics::splits_sub_admin(plans = plans, shp = data.frame(),
        sub_admin = munis)
}

