##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/11/07
## Purpose: redist constraints for a tidy workflow
##############################################


#######################
# constructors and reconstructors

# internal constructor
new_redist_constr <- function(constr = list(), data = tibble()) {
    constr <- reconstruct.redist_constr(constr, NULL)

    if (length(constr) > 0) {
        if (is.null(names(constr))) cli_abort("Null names.")
        if (any(names(constr) == "")) cli_abort("Empty names.")
        for (el in constr) {
            if (!is.list(el)) cli_abort("Not a nested list")
            classes <- vapply(el, class, character(1))
            if (length(classes) == 0 || any(classes != "list"))
                cli_abort("Not a nested list")
        }
    }

    stopifnot(is.data.frame(data))
    attr(constr, "data") <- data

    constr
}

validate_redist_constr <- function(constr) {
    if (!is.list(constr)) cli_abort("Not a list")
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")

    constr
}

#' Set up constraints for sampling
#'
#' `redist_constr` objects are used to specify constraints when sampling
#' redistricting plans with [redist_smc()] and [redist_mergesplit()]. Each
#' constraint is specified as a function which scores a given plan. Higher
#' scores are penalized and sampled less frequently.
#'
#' The `redist_constr` object keeps track of sampling constraints in a nested list.
#' You can view the exact structure of this list by calling [str()].
#' Constraints may be added by using one of the following functions:
#'
#' `r paste0("* [", setdiff(ls("package:redist")[grep("add_constr_", ls("package:redist"))], "add_constr_qps"), "()]", collapse="\n")`
#'
#' More information about each constraint can be found on the relevant constraint page.
#'
#' @param map a [redist_map()] object; the map that will be used in sampling
#'
#' @returns a `redist_constr` object, which is just a list with a certain nested structure.
#'
#' @examples
#' data(iowa)
#' map_ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
#' constr <- redist_constr(map_ia)
#' constr <- add_constr_splits(constr, strength = 1.5, admin = region)
#' print(constr)
#'
#' @md
#' @concept simulate
#' @export
redist_constr <- function(map = tibble()) {
    new_redist_constr(data = map)
}

# properly nest
reconstruct.redist_constr <- function(constr, old) {
    # handle deprecated constraint names
    if ("vra" %in% names(constr)) {
        .Deprecated("grp_pow", old = "vra")
        constr <- add_to_constr(constr, "grp_pow", constr$vra)
        constr$vra <- NULL
    }
    if ("hinge" %in% names(constr)) {
        .Deprecated("grp_hinge", old = "hinge")
        constr <- add_to_constr(constr, "grp_hinge", constr$hinge)
        constr$hinge <- NULL
    }

    for (i in seq_along(constr)) {
        classes <- vapply(constr[[i]], class, character(1))
        if (any(classes != "list")) {
            constr[[i]] <- list(constr[[i]])
        }
    }

    class(constr) <- c("redist_constr", "list")
    constr
}

#######################
# constraint helpers

# helper
add_to_constr <- function(constr, name, new_constr) {
    cur_constr <- constr[[name]]

    if (!is.null(cur_constr)) {
        classes <- vapply(cur_constr, class, character(1))
        if (any(classes != "list")) {
            constr[[name]] <- list(cur_constr, new_constr)
        } else {
            constr[[name]] <- c(cur_constr, list(new_constr))
        }
    } else {
        constr[[name]] <- list(new_constr)
    }

    constr
}

#' Sampling constraints
#'
#' The [redist_smc()] and [redist_mergesplit()] algorithms in this package allow
#' for additional constraints on the redistricting process to be encoded in the
#' target distribution for sampling. These functions are provided to specify
#' these constraints. All arguments are quoted and evaluated in the context of
#' the data frame provided to [redist_constr()].
#'
#' All constraints are fed into a Gibbs measure, with coefficients on each
#' constraint set by the corresponding `strength` parameter.
#' The strength can be any real number, with zero corresponding to no constraint.
#' Higher and higher `strength` values will eventually cause the algorithm's
#' accuracy and efficiency to suffer. Whenever you use constraints, be sure to
#' check all sampling diagnostics.
#'
#' The `status_quo` constraint adds a term measuring the variation of
#' information distance between the plan and the reference, rescaled to \[0, 1\].
#'
#' The `grp_hinge` constraint takes a list of target group percentages. It
#' matches each district to its nearest target percentage, and then applies a
#' penalty of the form \eqn{\sqrt{max(0, tgt - grouppct)}}, summing across
#' districts. This penalizes districts which are below their target percentage.
#' Use [plot.redist_constr()] to visualize the effect of this constraint and
#' calibrate `strength` appropriately.
#'
#' The `grp_inv_hinge` constraint takes a list of target group percentages. It
#' matches each district to its nearest target percentage, and then applies a
#' penalty of the form \eqn{\sqrt{max(0, grouppct - tgt)}}, summing across
#' districts. This penalizes districts which are above their target percentage.
#' Use [plot.redist_constr()] to visualize the effect of this constraint and
#' calibrate `strength` appropriately.
#'
#' The `grp_pow` constraint (for expert use) adds a term of the form
#' \eqn{(|tgtgroup-grouppct||tgtother-grouppct|)^{pow})}, which
#' encourages districts to have group shares near either `tgt_group`
#' or `tgt_other`.  Values of `strength` depend heavily on the values of these
#' parameters and especially the `pow` parameter.
#' Use [plot.redist_constr()] to visualize the effect of this constraint and
#' calibrate `strength` appropriately.
#'
#' The `compet` constraint encourages competitiveness by applying the `grp_pow`
#' constraint with target percentages set to 50%. For convenience, it is
#' specified with Democratic and Republican vote shares.
#'
#' The `incumbency` constraint adds a term counting the number of districts
#' containing paired-up incumbents.
#' Values of `strength` should generally be small, given that the underlying values are counts.
#'
#' The `splits` constraint adds a term counting the number of
#' counties which are split once or more.
#' Values of `strength` should generally be small, given that the underlying values are counts.
#'
#' The `multisplits` constraint adds a term counting the number of
#' counties which are split twice or more.
#' Values of `strength` should generally be small, given that the underlying values are counts.
#'
#' The `total_splits` constraint adds a term counting the total number of times
#' each county is split, summed across counties (i.e., counting the number of
#' excess district-county pairs). Values of `strength` should generally be
#' small, given that the underlying values are counts.
#'
#' The `edges_rem` constraint adds a term counting the number of edges removed from the
#' adjacency graph. This is only usable with `redist_flip()`, as other algorithms
#' implicitly use this via the `compactness` parameter. Values of `strength` should
#' generally be small, given that the underlying values are counts.
#'
#' The `log_st` constraint constraint adds a term counting the log number of spanning
#' trees. This is only usable with `redist_flip()`, as other algorithms
#' implicitly use this via the `compactness` parameter.
#'
#' The `polsby` constraint adds a term encouraging compactness as defined by the
#' Polsby Popper metric. Values of `strength` may be of moderate size.
#'
#' The `fry_hold` constraint adds a term encouraging compactness as defined by the
#' Fryer Holden metric. Values of `strength` should be extremely small, as the
#' underlying values are massive when the true minimum Fryer Holden denominator is not known.
#'
#' The `segregation` constraint adds a term encouraging segregation among minority groups,
#' as measured by the dissimilarity index.
#'
#' The `pop_dev` constraint adds a term encouraging plans to have smaller population deviations
#' from the target population.
#'
#' The `custom` constraint allows the user to specify their own constraint using
#' a function which evaluates districts one at a time. The provided function
#' `fn` should take two arguments: a vector describing the current plan
#' assignment for each unit as its first argument, and an integer describing the
#' district which to evaluate in the second argument. `which([plans == distr])`
#' would give the indices of the units that are assigned to a district `distr`
#' in any iteration. The function must return a single scalar for each plan -
#' district combination, where a value of 0 indicates no penalty is applied. If
#' users want to penalize an entire plan, they can have the penalty function
#' return a scalar that does not depend on the district. It is important that
#' `fn` not use information from precincts not included in `distr`, since in the
#' case of SMC these precincts may not be assigned any district at all (`plan`
#' will take the value of 0 for these precincts). The flexibility of this
#' constraint comes with an additional computational cost, since the other
#' constraints are written in C++ and so are more performant.
#'
#' @param constr A [redist_constr()] object
#' @param strength The strength of the constraint. Higher values mean a more restrictive constraint.
#' @param score_districts_only Whether or not to apply the constraints to
#' districts only. If constraints are only applied to districts then it will
#' likely cause a drop in efficiency in the final round. If splitting plans all
#' the way this does not affect the final target distribution.
#'
#' @examples
#' data(iowa)
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
#' constr <- redist_constr(iowa_map)
#' constr <- add_constr_splits(constr, strength = 1.5, admin = name)
#' constr <- add_constr_grp_hinge(constr, strength = 100,
#'     dem_08, tot_08, tgts_group = c(0.5, 0.6))
#' # encourage districts to have the same number of counties
#' constr <- add_constr_custom(constr, strength = 1000, fn = function(plan, distr) {
#'     # notice that we only use information on precincts in `distr`
#'     abs(sum(plan == distr) - 99/4)
#' })
#' print(constr)
#'
#' @md
#' @concept simulate
#' @name constraints
NULL

#' @param current The reference map for the status quo constraint.
#' @rdname constraints
#' @export
add_constr_status_quo <- function(constr, strength, current,
                                  score_districts_only = TRUE ) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    is.logical(score_districts_only)
    data <- attr(constr, "data")
    if (missing(current)) current <- get_existing(data)

    new_constr <- list(strength = strength,
        current = eval_tidy(enquo(current), data),
        score_districts_only = score_districts_only)
    if (is.null(current) || length(new_constr$current) != nrow(data))
        cli_abort("{.arg current} must be provided, and must have as many
                  precincts as the {.cls redist_map}")
    new_constr$n_current <- max(new_constr$current)

    add_to_constr(constr, "status_quo", new_constr)
}

#' @param group_pop A vector of group population
#' @param total_pop A vector of total population. Defaults to the population vector used for sampling.
#' @param tgt_group,tgt_other Target group shares for the power-type constraint.
#' @param pow The exponent for the power-type constraint.
#'
#' @rdname constraints
#' @export
add_constr_grp_pow <- function(constr, strength, group_pop, total_pop = NULL,
                               tgt_group = 0.5, tgt_other = 0.5, pow = 1.0,
                               score_districts_only = FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    is.logical(score_districts_only)
    data <- attr(constr, "data")

    new_constr <- list(strength = strength, score_districts_only=score_districts_only,
        group_pop = eval_tidy(enquo(group_pop), data),
        total_pop = eval_tidy(enquo(total_pop), data),
        tgt_group = tgt_group, tgt_other = tgt_other, pow = pow)
    if (is.null(new_constr$total_pop)) {
        if (!is.null(attr(data, "pop_col"))) {
            new_constr$total_pop <- data[[attr(data, "pop_col")]]
        } else {
            cli_abort("{.arg total_pop} missing.")
        }
    }

    stopifnot(length(new_constr$group_pop) == nrow(data))
    stopifnot(length(new_constr$total_pop) == nrow(data))

    add_to_constr(constr, "grp_pow", new_constr)
}

#' @param tgts_group A vector of target group shares for the hinge-type constraint.
#' @rdname constraints
#' @export
add_constr_grp_hinge <- function(constr, strength, group_pop, total_pop = NULL,
                                 tgts_group = c(0.55),
                                 score_districts_only = FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    # if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    is.logical(score_districts_only)
    data <- attr(constr, "data")

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        group_pop = eval_tidy(enquo(group_pop), data),
        total_pop = eval_tidy(enquo(total_pop), data),
        tgts_group = tgts_group
        )
    if (is.null(new_constr$total_pop)) {
        if (!is.null(attr(data, "pop_col"))) {
            new_constr$total_pop <- data[[attr(data, "pop_col")]]
        } else {
            cli_abort("{.arg total_pop} missing.")
        }
    }

    stopifnot(length(new_constr$group_pop) == nrow(data))
    stopifnot(length(new_constr$total_pop) == nrow(data))

    add_to_constr(constr, "grp_hinge", new_constr)
}


#' @param tgts_group A vector of target group shares for the hinge-type constraint.
#' @rdname constraints
#' @export
add_constr_grp_inv_hinge <- function(constr, strength, group_pop, total_pop = NULL,
                                     tgts_group = c(0.55),
                                     score_districts_only = FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    # if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data <- attr(constr, "data")

    new_constr <- list(strength = strength,
                       score_districts_only = score_districts_only,
        group_pop = eval_tidy(enquo(total_pop), data) - eval_tidy(enquo(group_pop), data),
        total_pop = eval_tidy(enquo(total_pop), data),
        tgts_group = tgts_group)
    if (is.null(new_constr$total_pop)) {
        if (!is.null(attr(data, "pop_col"))) {
            new_constr$total_pop <- data[[attr(data, "pop_col")]]
        } else {
            cli_abort("{.arg total_pop} missing.")
        }
    }

    stopifnot(length(new_constr$group_pop) == nrow(data))
    stopifnot(length(new_constr$total_pop) == nrow(data))

    add_to_constr(constr, "grp_inv_hinge", new_constr)
}
#' @param dvote,rvote A vector of Democratic or Republican vote counts
#' @rdname constraints
#' @export
add_constr_compet <- function(constr, strength, dvote, rvote, pow = 0.5,
                              score_districts_only = FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data <- attr(constr, "data")

    new_constr <- list(strength = strength,
                       score_districts_only = score_districts_only,
        dvote = eval_tidy(enquo(dvote), data),
        rvote = eval_tidy(enquo(rvote), data),
        pow = pow)
    stopifnot(length(new_constr$dvote) == nrow(data))
    stopifnot(length(new_constr$rvote) == nrow(data))

    add_to_constr(constr, "compet", new_constr)
}

#' @param incumbents A vector of unit indices for incumbents. For example, if
#' three incumbents live in the precincts that correspond to rows 1, 2, and
#' 100 of your [redist_map], entering incumbents = c(1, 2, 100) would avoid
#' having two or more incumbents be in the same district.
#' @rdname constraints
#' @export
add_constr_incumbency <- function(constr, strength, incumbents,
                                  score_districts_only = TRUE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data <- attr(constr, "data")

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        incumbents = eval_tidy(enquo(incumbents), data))

    add_to_constr(constr, "incumbency", new_constr)
}

#' @param admin A vector indicating administrative unit membership
#' @rdname constraints
#' @export
add_constr_splits <- function(constr, strength, admin, score_districts_only=TRUE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data <- attr(constr, "data")

    admin <- eval_tidy(enquo(admin), data)
    if (is.null(admin)) {
        cli_abort("{.arg admin} may not be {.val NULL}.")
    }
    if (any(is.na(admin))) {
        cli_abort("{.arg admin} many not contain {.val NA}s.")
    }
    admin <- vctrs::vec_group_id(admin)

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        admin = admin,
        n = length(unique(admin)))

    add_to_constr(constr, "splits", new_constr)
}

#' @rdname constraints
#' @export
add_constr_multisplits <- function(constr, strength, admin,
                                   score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data <- attr(constr, "data")

    admin <- eval_tidy(enquo(admin), data)
    if (is.null(admin)) {
        cli_abort("{.arg admin} may not be {.val NULL}.")
    }
    if (any(is.na(admin))) {
        cli_abort("{.arg admin} many not contain {.val NA}s.")
    }

    admin <- vctrs::vec_group_id(admin)

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        admin = admin,
        n = length(unique(admin)))
    add_to_constr(constr, "multisplits", new_constr)
}

#' @rdname constraints
#' @export
add_constr_total_splits <- function(constr, strength, admin,
                                    score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data <- attr(constr, "data")

    admin <- eval_tidy(enquo(admin), data)
    if (is.null(admin)) {
        cli_abort("{.arg admin} may not be {.val NULL}.")
    }
    if (any(is.na(admin))) {
        cli_abort("{.arg admin} many not contain {.val NA}s.")
    }

    admin <- vctrs::vec_group_id(admin)

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        admin = admin,
        n = length(unique(admin)))
    add_to_constr(constr, "total_splits", new_constr)
}

#' @rdname constraints
#' @export
add_constr_pop_dev <- function(constr, strength, score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results.")
    data <- attr(constr, "data")

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only)
    add_to_constr(constr, "pop_dev", new_constr)
}

#' @rdname constraints
#' @export
add_constr_segregation <- function(constr, strength, group_pop, total_pop = NULL,
                                   score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results.")
    data <- attr(constr, "data")

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        group_pop = eval_tidy(enquo(group_pop), data),
        total_pop = eval_tidy(enquo(total_pop), data))
    if (is.null(new_constr$total_pop)) {
        if (!is.null(attr(data, "pop_col"))) {
            new_constr$total_pop <- data[[attr(data, "pop_col")]]
        } else {
            cli_abort("{.arg total_pop} missing.")
        }
    }
    if (is.null(new_constr$group_pop)) {
        cli_abort("{.arg group_pop} missing.")
    }

    stopifnot(length(new_constr$group_pop) == nrow(data))
    stopifnot(length(new_constr$total_pop) == nrow(data))
    add_to_constr(constr, "segregation", new_constr)
}

#' @param perim_df A dataframe output from `redistmetrics::prep_perims`
#' @rdname constraints
#' @export
add_constr_polsby <- function(constr, strength, perim_df = NULL,
                              score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results.")
    data <- attr(constr, "data")

    if (!inherits(data, "sf")) {
        cli_abort("Input to {.fun redist_constr} must be a {.cls sf} object.")
    }

    areas <- sf::st_area(data)

    if (is.null(perim_df)) {
        perim_df <- redistmetrics::prep_perims(data)
    }

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        from = perim_df$origin,
        to = perim_df$touching,
        area = areas,
        perimeter = perim_df$edge)

    add_to_constr(constr, "polsby", new_constr)
}

#' @rdname constraints
#' @param ssdmat Squared distance matrix for Fryer Holden constraint
#' @param denominator Fryer Holden minimum value to normalize by. Default is 1 (no normalization).
#' @export
add_constr_fry_hold <- function(constr, strength, total_pop = NULL, ssdmat = NULL, denominator = 1,
                                score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results.")
    data <- attr(constr, "data")

    total_pop <- eval_tidy(enquo(total_pop), data)
    if (is.null(total_pop)) {
        if (!is.null(attr(data, "pop_col"))) {
            total_pop <- data[[attr(data, "pop_col")]]
        } else {
            cli_abort("{.arg total_pop} missing.")
        }
    }
    if (is.null(ssdmat)) {
        ssdmat <- calcPWDh(sf::st_coordinates(sf::st_centroid(data)))
    }

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        total_pop = total_pop,
        ssdmat = ssdmat,
        denominator = denominator)


    add_to_constr(constr, "fry_hold", new_constr)
}

#' @rdname constraints
#' @export
add_constr_log_st <- function(constr, strength, admin = NULL,
                              score_districts_only = FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results.")
    data <- attr(constr, "data")

    admin <- eval_tidy(enquo(admin), data)
    if (is.null(admin)) {
        admin <- rep(1, nrow(data))
    }
    if (any(is.na(admin))) {
        cli_abort("{.arg admin} many not contain {.val NA}s.")
    }

    admin <- vctrs::vec_group_id(admin)

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        admin = admin)

    add_to_constr(constr, "log_st", new_constr)
}

#' @rdname constraints
#' @export
add_constr_edges_rem <- function(constr, strength, score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results.")
    data <- attr(constr, "data")

    new_constr <- list(strength = strength, score_districts_only=FALSE)

    add_to_constr(constr, "edges_removed", new_constr)
}

#' @param cities A vector containing zero entries for non-cities and non-zero entries for each city for `qps`.
#' @noRd
add_constr_qps <- function(constr, strength, cities, total_pop = NULL,
                           score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results.")
    data <- attr(constr, "data")

    new_constr <- list(strength = strength,
                       score_districts_only=score_districts_only,
        cities = eval_tidy(enquo(cities), data))
    new_constr$n_cty <- max(new_constr$cities) + 1

    cli::cli_inform("The QPS constraint is not officially supported and may disappear.",
        .frequency = "once")
    add_to_constr(constr, "qps", new_constr)
}

# utilty functions for parsing ASTs
find_env <- function(name, env = rlang::caller_env()) {
    if (identical(env, rlang::empty_env())) {
        NULL
    } else if (rlang::env_has(env, name)) {
        env
    } else {
        find_env(name, rlang::env_parent(env))
    }
}
extract_vars = function(expr) {
    if (rlang::is_syntactic_literal(expr)) {
        NULL
    } else if (rlang::is_symbol(expr)) {
        rlang::as_string(expr)
    } else if (rlang::is_call(expr)) {
        c(rlang::as_string(expr[[1]]), unlist(sapply(expr[-1], extract_vars)))
    }
}


#' @param fn A function
#' @rdname constraints
#' @export
add_constr_custom <- function(constr, strength, fn, score_districts_only=FALSE) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")

    args <- rlang::fn_fmls(fn)
    if (length(args) != 2) cli_abort("Function must take exactly two arguments.")

    constr_env = rlang::fn_env(fn)
    constr_env <- rlang::env(constr_env)
    # every symbol used in the function (except the 2 arguments)
    var_names = setdiff(
        all.names(rlang::fn_body(fn)),
        names(args)
    )

    for (nm in var_names) {
        found = find_env(nm, constr_env)
        if (!is.null(found) &&
                !identical(found, rlang::base_env()) &&
                !identical(found, constr_env) &&
                !identical(found, rlang::pkg_env("redist"))) {
            constr_env[[nm]] = get(nm, envir=found)
        }
    }

    if (!is.null(plan <- get_existing(attr(constr, "data")))) {
        out <- tryCatch(fn(plan, min(plan)), error = function(e) {
            cli_abort(c("Ran into an error testing custom constraint
                        on the existing plan:",
                "x" = e$message))
        })
        if (!is.numeric(out) || length(out) != 1 || !is.finite(out))
            cli_abort(c("Evaluting custom constraint on the existing plan failed.",
                "*" = "The constraint function should return a single scalar value.",
                "*" = "Make sure that your constraint function tests all edge cases
                             and never returns {.val {NA}} or {.val {Inf}}."))
    }

    rlang::fn_env(fn) <- constr_env

    new_constr <- list(strength = strength, score_districts_only=score_districts_only,
                       fn = fn)
    add_to_constr(constr, "custom", new_constr)
}

#######################
# generics

#' Generic to print redist_constr
#' @param x redist_constr
#' @param header if FALSE, then suppress introduction / header line
#' @param details if FALSE, then suppress the details of each constraint
#' @param \dots additional arguments
#' @method print redist_constr
#' @return Prints to console and returns input redist_constr
#' @export
print.redist_constr <- function(x, header = TRUE, details = TRUE, ...) {
    if (header)
        cli_text("A {.cls redist_constr} with {length(x)} constraint{?s}")

    print_constr <- function(x) {
        if (details) {
            idx_strength <- which(names(x) == "strength")
            str(x[-idx_strength], no.list = T, comp.str = "   ", give.attr = FALSE)
        }
    }

    x <- unlist(x, recursive = FALSE)
    for (nm in names(x)) {
        if("score_districts_only" %in% x[[nm]] && x[[nm]]$score_districts_only){
            score_str <- "districts only"
        }else{
            score_str <- "all regions"
        }
        if (startsWith(nm, "status_quo")) {
            cli::cli_bullets(c("*" = "A status quo constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "grp_pow")) {
            cli::cli_bullets(c("*" = "A (power-type) group share constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "grp_hinge")) {
            cli::cli_bullets(c("*" = "A (hinge-type) group share constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "grp_inv_hinge")) {
            cli::cli_bullets(c("*" = "An (inverse-hinge-type) group share constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "compet")) {
            cli::cli_bullets(c("*" = "A competitiveness constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "incumbency")) {
            cli::cli_bullets(c("*" = "An incumbency constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "splits")) {
            cli::cli_bullets(c("*" = "A splits constraint of strength {x[[nm]]$strength} applied to {score_str}"))
        } else if (startsWith(nm, "multisplits")) {
            cli::cli_bullets(c("*" = "A multisplits constraint of strength {x[[nm]]$strength} applied to {score_str}"))
        } else if (startsWith(nm, "total_splits")) {
            cli::cli_bullets(c("*" = "A total splits constraint of strength {x[[nm]]$strength} applied to {score_str}"))
        } else if (startsWith(nm, "custom")) {
            cli::cli_bullets(c("*" = "A custom constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "edges_rem")) {
            cli::cli_bullets(c("*" = "An (edges-removed-type) compactness constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "log_st")) {
            cli::cli_bullets(c("*" = "A (log-spanning-tree-type) compactness constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "polsby")) {
            cli::cli_bullets(c("*" = "A (Polsby-Popper-type) compactness constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "fry_hold")) {
            cli::cli_bullets(c("*" = "A (Fryer-Holden-type) compactness constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "pop_dev")) {
            cli::cli_bullets(c("*" = "A population deviation constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "segregation")) {
            cli::cli_bullets(c("*" = "A dissimilarity segregation constraint of strength {x[[nm]]$strength} applied to {score_str}"))
            print_constr(x[[nm]])
        } else {
            cli::cli_bullets(c("*" = "An unknown constraint {.var {nm}} applied to {score_str}"))
            print_constr(x[[nm]])
        }
    }
}


#' Visualize constraints
#'
#' Plots the constraint strength versus some running variable. Currently
#' supports visualizing the `grp_hinge`, `grp_inv_hinge`, and `grp_pow`
#' constraints.
#'
#' @param x A [redist_constr] object.
#' @param y Ignored.
#' @param type What type of constraint to visualize. Currently supports only
#' `"group"`, for visualizing constraint strength by group share.
#' @param xlim Range of group shares to visualize.
#' @param \dots additional arguments (ignored)
#'
#' @method plot redist_constr
#' @return A ggplot object
#'
#' @examples
#' data(iowa)
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
#' constr <- redist_constr(iowa_map)
#' constr <- add_constr_grp_hinge(constr, strength = 30,
#'                                dem_08, tot_08, tgts_group = 0.5)
#' constr <- add_constr_grp_hinge(constr, strength = -20,
#'                                dem_08, tot_08, tgts_group = 0.3)
#' plot(constr)
#'
#' @concept prepare
#' @export
plot.redist_constr <- function(x, y, type="group", xlim=c(0, 1), ...) {
    if (type != "group") cli_abort("Only {.arg type = \"group\"} is currently supported.")

    out <- tibble(share = seq(xlim[1], xlim[2], by = .001),
                  penalty = 0)

    if ("grp_pow" %in% names(x)) {
        for (obj in x$grp_pow) {
            out$penalty = out$penalty + obj$strength * (
                abs(out$share - obj$tgt_group) * abs(out$share - obj$tgt_other)
                )^obj$pow
        }
    }

    warn_multiple = FALSE
    if ("grp_hinge" %in% names(x)) {
        for (obj in x$grp_hinge) {
            if (length(obj$tgts_group) > 1) warn_multiple = TRUE
            out$penalty = out$penalty + obj$strength * sqrt(pmax(0.0, obj$tgts_group[1] - out$share))
        }
    }

    if ("grp_inv_hinge" %in% names(x)) {
        for (obj in x$grp_inv_hinge) {
            if (length(obj$tgts_group) > 1) warn_multiple = TRUE
            out$penalty = out$penalty + obj$strength * sqrt(pmax(0.0, out$share - obj$tgts_group[1]))
        }
    }

    if (warn_multiple) {
        cli_warn("Multiple group-share targets found; only plotting first.")
    }

    ggplot(out, aes(x=.data$share, y=.data$penalty)) +
        geom_path() +
        ggplot2::scale_x_continuous("Group share of district population",
                                    labels=function(x) paste0(round(100*x), "%")) +
        labs(y = "Penalty")
}

#' @method str redist_constr
#' @export
str.redist_constr <- function(object, give.attr = FALSE, ...) {
    NextMethod("str", object, give.attr = give.attr, ...)
}

#' @method as.list redist_constr
#' @export
as.list.redist_constr <- function(x, ...) {
    class(x) <- "list"
    attr(x, "data") <- NULL
    x
}
