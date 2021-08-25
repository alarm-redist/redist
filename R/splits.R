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
    if (missing(plans)) {
        stop('Please provide an argument to plans.')
    }
    if (inherits(plans, 'redist_plans')) {
        plans <- get_plans_matrix(plans)
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol = 1)
    }
    if (!any(class(plans) %in% c('numeric', 'matrix'))) {
        stop('Please provide "plans" as a matrix.')
    }

    if (missing(counties)) {
        stop('Please provide an argument to counties.')
    }
    if (class(counties) %in% c('character', 'numeric', 'integer')) {
        uc <- unique(sort(counties))
        county_id <- rep(0, nrow(plans))
        for (i in 1:nrow(plans)) {
            county_id[i] <- which(uc == counties[i])
        }
    } else{
        stop('Please provide "counties" as a character, numeric, or integer vector.')
    }


    splits(plans, community = county_id, length(unique(plans[, 1])), 1)
}

#' Identify which counties are split by a plan
#'
#' @param plan A vector of precinct/unit assignments
#' @param counties A vector of county names or county ids.
#'
#' @return A logical vector which is \code{TRUE} for precincts belonging to
#'   counties which are split
#'
#' @concept analyze
#' @export
is_county_split = function(plan, counties) {
    counties = as.integer(as.factor(counties))
    (tapply(plan, counties, FUN=function(y) length(unique(y))) > 1)[counties]
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
#' splits <- redist.district.splits(plans, ia$region)
redist.district.splits <- function(plans, counties) {
    if (missing(plans)) {
        stop('Please provide an argument to plans.')
    }
    if (inherits(plans, 'redist_plans')) {
        plans <- get_plans_matrix(plans)
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol = 1)
    }
    if (!any(class(plans) %in% c('numeric', 'matrix'))) {
        stop('Please provide "plans" as a matrix.')
    }

    if (missing(counties)) {
        stop('Please provide an argument to counties.')
    }
    if (class(counties) %in% c('character', 'numeric', 'integer')) {
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
#' splits <- redist.multisplits(plans, ia$region)
redist.multisplits <- function(plans, counties) {
    if (missing(plans)) {
        stop('Please provide an argument to plans.')
    }
    if (inherits(plans, 'redist_plans')) {
        plans <- get_plans_matrix(plans)
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol = 1)
    }
    if (!any(class(plans) %in% c('numeric', 'matrix'))) {
        stop('Please provide "plans" as a matrix.')
    }

    if (missing(counties)) {
        stop('Please provide an argument to counties.')
    }
    if (class(counties) %in% c('character', 'numeric', 'integer')) {
        county_id <- redist.county.id(counties)
    } else {
        stop('Please provide "counties" as a character, numeric, or integer vector.')
    }

    splits(plans - 1, community = county_id - 1, length(unique(plans[, 1])), 2)
}
