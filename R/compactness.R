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
#' This can be created in advance with \code{redist.prep.polsbypopper}.
#' @param perim_df A dataframe output from \code{redist.prep.polsbypopper}
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
#' redist.compactness(
#'     shp = fl25, plans = plans_05[, 1:3],
#'     measure = c("PolsbyPopper", "EdgesRemoved")
#' )
#' @export redist.compactness
redist.compactness <- function(shp = NULL,
                               plans,
                               measure = c("PolsbyPopper"),
                               total_pop = NULL, adj = NULL, draw = 1,
                               ncores = 1, counties = NULL, planarize = 3857,
                               ppRcpp, perim_path, perim_df) {
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
                shp <- st_transform(shp, planarize)
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
        cli_inform("{.pkg redist} 4.0.0 fixes a bug where Schwartzberg was reporting 1/Schwartzberg. Current values are correct.",
            .frequency = "once", .frequency_id = "schwartz"
        )
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



#' Prep Polsby Popper Perimeter Dataframe
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required unless "EdgesRemoved"
#' and "logSpanningTree" with adjacency provided.
#' @param planarize a number, indicating the CRS to project the shapefile to if
#' it is latitude-longitude based. Set to FALSE to avoid planarizing.
#' @param perim_path A path to save an Rds
#' @param ncores the number of cores to parallelize over
#'
#' @return A perimeter dataframe
#' @export
#'
#' @importFrom sf st_buffer st_is_valid st_geometry<- st_touches st_transform
#' @examples
#' data(fl25)
#' perim_df <- redistmetrics::prep_perims(shp = fl25)
redist.prep.polsbypopper <- function(shp, planarize = 3857, perim_path, ncores = 1) {
    .Deprecated(new = "redistmetrics::prep_perims()")
    if (missing(shp)) {
        cli_abort("Please provide an argument to {.arg shp}.")
    }
    redistmetrics::prep_perims(
        shp = shp, epsg = planarize, perim_path = perim_path,
        ncores = ncores
    )
}
