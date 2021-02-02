##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################

#######################
# constructors and reconstructors

#' A set of redistricting plans
#'
#' A \code{redist_plans} object is essentially a data frame of summary
#' information on each district and each plan, along with the matrix of district
#' assignments and information about the simulation process used to generate the
#' plans.
#'
#' The first two columns of the data frame will be \code{draw}, a factor indexing
#' the simulation draw, and \code{district}, an integer indexing the districts
#' within a plan. The data frame will therefore have \code{n_sims*n_distr} rows.
#' The data frame will at all times be grouped at least by \code{draw}.
#' As a data frame, the usual \code{\link{dplyr}} methods will work.
#'
#' Other useful methods for \code{redist_plans} objects:
#' * \code{\link{add_reference}}
#' * \code{\link{pullback}}
#' * \code{\link{get_plan_matrix}}
#' * \code{\link{get_plan_weights}}
#' * \code{\link{get_sampling_info}}
#' * \code{\link{as.matrix.redist_plans}}
#' * \code{\link{plot.redist_plans}}
#'
#' @name redist_plans
#' @md
NULL

# plans has n_precinct columns and n_sims rows
# map is a redist_map
# algorithm is one of "smc" or "mcmc"
# wgt is the weights before any resampling or truncation
# ... will depend on the algorithm
new_redist_plans = function(plans, map, algorithm, wgt, resampled=TRUE, ...) {
    n_sims = ncol(plans)
    stopifnot(n_sims >= 1)

    n_prec = nrow(plans)
    n_distr = attr(map, "n_distr")

    prec_pop = map[[attr(map, "pop_col")]]
    distr_pop = pop_tally(plans, prec_pop, n_distr)

    attr_names = c("redist_attr", "plans", "algorithm", "wgt", "resampled",
                   "merge_idx", "pop", names(list(...)))

    x = structure(
        dplyr::group_by(
            tibble::tibble(draw = rep(as.factor(1:n_sims), each=n_distr),
                           district = rep(1:n_distr, n_sims),
                           pop = as.numeric(distr_pop)),
            .data$draw),
        plans=plans, algorithm=algorithm, wgt=wgt,
        resampled=resampled, merge_idx=attr(map, "merge_idx"),
        pop=prec_pop, redist_attr=attr_names, ...,
        class=c("redist_plans", "grouped_df", "tbl_df", "tbl", "data.frame")
    )

    if (!is.null(attr(map, "existing_col")))
        add_reference(x, get_existing(map))
    else
        x
}

validate_redist_plans = function(x) {
    stopifnot(names(x)[1] == "draw")
    stopifnot(is.factor(x$draw))
    stopifnot(names(x)[2] == "district")
    stopifnot(is.integer(x$district))

    x
}

reconstruct.redist_plans = function(data, old) {
    if (!missing(old)) {
        for (name in attr(old, "redist_attr")) {
            if (is.null(attr(data, name)))
                attr(data, name) = attr(old, name)
        }

        if (is.null(attr(data, "groups")))
            attr(data, "groups") = attr(old, "groups")
    }

    class(data) = c("redist_plans", "grouped_df", "tbl_df", "tbl", "data.frame")

    data
}


#################
# Getters

#' Extract the matrix of district assignments from a redistricting simulation
#'
#' @param x the \code{redist_plans} object
#' @param ... ignored
#'
#' @export
get_plan_matrix = function(x) {
    stopifnot(inherits(x, "redist_plans"))
    attr(x, "plans")
}
# internal -- no check performed!
set_plan_matrix = function(x, mat) {
    attr(x, "plans") = mat
    x
}

#' Extract the sampling weights from a redistricting simulation
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns the weights, with an additional attribute \code{resampled}
#'   indicating whether the plans have been resampled according to these
#'   weights.
#'
#' @export
get_plan_weights = function(plans) {
    stopifnot(inherits(plans, "redist_plans"))
    wgt = attr(plans, "wgt")
    attr(wgt, "resampled") = attr(plans, "resampled")
    wgt
}

#' Extract the sampling information from a redistricting simulation
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns a list of parameters and information about the sampling problem.
#'
#' @export
get_sampling_info = function(plans) {
    stopifnot(inherits(plans, "redist_plans"))
    all_attr = attributes(plans)

    all_attr$names = NULL
    all_attr$row.names = NULL
    all_attr$class = NULL
    all_attr$plans = NULL
    all_attr$redist_attr = NULL

    all_attr
}

#######################
# Other useful functions

#' Add a reference plan to a set of plans
#'
#' This function facilitates comparing an existing (i.e., non-simulated)
#' redistricting plan to a set of simulated plans.
#'
#' @param plans a \code{redist_plans} object
#' @param ref_plan an integer vector containing the reference plan. It will be
#' renumbered to 1..\code{n_distr}.
#' @param name a human-readable name for the referece plan.
#'
#' @returns a modified \code{redist_plans} object containing the reference plan
add_reference = function(plans, ref_plan, name="<ref>") {
    stopifnot(inherits(plans, "redist_plans"))
    stopifnot(is.character(name))

    plan_m = get_plan_matrix(plans)
    stopifnot(is.numeric(ref_plan))
    stopifnot(length(ref_plan) == nrow(plan_m))

    ref_plan = as.integer(as.factor(ref_plan))
    n_distr = max(ref_plan)
    stopifnot(n_distr == max(plan_m[,1]))

    # first the matrix
    plan_m = cbind(ref_plan, plan_m)
    colnames(plan_m)[1] = name

    # then the dataframe
    distr_pop = pop_tally(matrix(ref_plan, ncol=1), attr(plans, "pop"), n_distr)
    fct_levels = c(name, levels(plans$draw))
    new_draw = rep(factor(fct_levels, levels=fct_levels), each=n_distr)
    x = dplyr::bind_rows(
            tibble::tibble(district = 1:n_distr,
                           pop = as.numeric(distr_pop)),
            plans[,-1] # 1 is 'draw' by defn
        ) %>%
        dplyr::mutate(draw = new_draw, .before="district")

    reconstruct.redist_plans(x, set_plan_matrix(plans, plan_m))
}

#' Pull back plans to unmerged units
#'
#' Merging map units through \code{\link{merge_by}} or \code{\link{summarize}}
#' changes the indexing of each unit.  Use this function to take a set of
#' redistricting plans from a \code{redist} algorithm and re-index them to
#' be compatible with the original set of units.
#'
#' @param plans a \code{redist_plans} object
#'
#' @returns a new, re-indexed, \code{redist_plans} object
#'
#' @export
pullback = function(plans) {
    stopifnot(inherits(plans, "redist_plans"))

    merge_idx = attr(plans, "merge_idx")
    if (is.null(merge_idx)) {
        warning("No merged indexing found.")
        return(plans)
    }

    attr(plans, "merge_idx") = NULL
    set_plan_matrix(plans, get_plan_matrix(plans)[merge_idx,])
}


#######################
# generics

#' @rdname get_plan_matrix
#' @method as.matrix redist_plans
#' @export
as.matrix.redist_plans = function(x, ...) get_plan_matrix(x)

#' @method dplyr_row_slice redist_plans
#' @export
dplyr_row_slice.redist_plans = function(data, i, ...) {
    if (is.logical(i)) i = which(i)

    draws_left = unique(as.integer(data$draw[i]))
    y = vctrs::vec_slice(data, i)
    plans_m = get_plan_matrix(data)

    # if we don't have every district present in every row
    # this check is necessary but not sufficient for what we want
    distrs = table(y$district)
    n_draws = length(draws_left)
    n_distr = max(plans_m[,1])
    if (!all.equal(range(distrs), rep(n_draws, 2)) || length(distrs) != n_distr)
        warning("Some districts may have been dropped. ",
                "This will prevent summary statistics from working correctly.\n",
                "To avoid this message, coerce using `as_tibble`.")

    set_plan_matrix(y, plans_m[,draws_left])
}

# 'template' is the old df
#' @method dplyr_reconstruct redist_plans
#' @export
dplyr_reconstruct.redist_plans = function(data, template) {
    reconstruct.redist_plans(data, template)
}

#' @method print redist_plans
#' @importFrom utils str
#' @export
print.redist_plans = function(x, ...) {
    plans_m = get_plan_matrix(x)
    n_samp = ncol(plans_m)
    if (is.null(colnames(plans_m))) {
        n_ref = 0
    } else {
        n_ref = sum(nchar(colnames(plans_m))>0)
        n_samp = n_samp - n_ref
    }

    cat(cli::pluralize("{n_samp} sampled plan{?s} "))
    if (n_ref > 0)
        cat(cli::pluralize("and {n_ref} reference plan{?s} "))
    cli::cat_line("with ", max(plans_m[,1]), " districts from a ",
        nrow(plans_m), "-unit map,\n  drawn using ",
        c(mcmc="Markov chain Monte Carlo",
          smc="Sequential Monte Carlo")[attr(x, "algorithm")])

    merge_idx = attr(x, "merge_idx")
    if (!is.null(merge_idx))
        cli::cat_line("Merged from another map with reindexing:",
                      utils::capture.output(str(merge_idx, vec.len=2)))

    if (attr(x, "resampled"))
        cli::cat_line("With plans resampled from weights")
    else
        cli::cat_line("With plans not resampled from weights")

    cli::cat_line("Plans matrix:", utils::capture.output(str(plans_m, give.attr=F)))

    utils::getS3method("print", "tbl")(x)

    invisible(x)
}

#' Summary plots for \code{redist_plans}
#'
#' @param x the \code{redist_plans} object.
#' @param y ignored
#' @param ... ignored
#'
#' @export
plot.redist_plans = function(x, y, ...) {
    ggplot(NULL, aes(x=x$pop)) +
        ggplot2::geom_histogram(bins=1+ceiling(ncol(get_plan_matrix(x)) / 10)) +
        ggplot2::labs(x="District population", y=NULL)
}
