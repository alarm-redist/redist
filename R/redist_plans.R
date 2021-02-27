##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################


# constructors and reconstructors -----------------------------------------

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
#' As a data frame, the usual \code{\link{dplyr}} methods will work.
#'
#' Other useful methods for \code{redist_plans} objects:
#' * \code{\link{add_reference}}
#' * \code{\link{pullback}}
#' * \code{\link{number_by}}
#' * \code{\link{match_numbers}}
#' * \code{\link{prec_assignment}}
#' * \code{\link{get_plan_matrix}}
#' * \code{\link{get_plan_weights}}
#' * \code{\link{get_sampling_info}}
#' * \code{\link{subset_sampled}}
#' * \code{\link{subset_ref}}
#' * \code{\link{as.matrix.redist_plans}}
#' * \code{\link{plot.redist_plans}}
#'
#' @name redist_plans
#' @concept analyze
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

    structure(tibble::tibble(draw = rep(as.factor(1:n_sims), each=n_distr),
                             district = rep(1:n_distr, n_sims),
                             pop = as.numeric(distr_pop)),
              plans=plans, algorithm=algorithm, wgt=wgt,
              resampled=resampled, merge_idx=attr(map, "merge_idx"),
              pop=prec_pop, redist_attr=attr_names, ...,
              class=c("redist_plans", "tbl_df", "tbl", "data.frame"))
}

validate_redist_plans = function(x) {
    stopifnot(names(x)[1] == "draw")
    stopifnot(is.factor(x$draw))

    x
}

reconstruct.redist_plans = function(data, old) {
    if (colnames(data)[1] != "draw")
        return(data)

    if (!missing(old)) {
        for (name in attr(old, "redist_attr")) {
            if (is.null(attr(data, name)))
                attr(data, name) = attr(old, name)
        }
    }

    classes = c("tbl_df", "tbl", "data.frame")
    if (inherits(data, "grouped_df"))
        classes = c("grouped_df", classes)

    class(data) = c("redist_plans", classes)

    data
}




# getters / setters ------------------------------------------------------------


#' Extract the matrix of district assignments from a redistricting simulation
#'
#' @param x the \code{redist_plans} object
#' @param ... ignored
#'
#' @concept analyze
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
#' @concept analyze
#' @export
get_plan_weights = function(plans) {
    stopifnot(inherits(plans, "redist_plans"))
    wgt = attr(plans, "wgt")
    attr(wgt, "resampled") = attr(plans, "resampled")
    wgt
}

get_n_ref = function(x) {
    stopifnot(inherits(x, "redist_plans"))
    plans_m = get_plan_matrix(x)
    if (is.null(colnames(plans_m))) 0 else sum(nchar(colnames(plans_m))>0)
}

#' Extract the sampling information from a redistricting simulation
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns a list of parameters and information about the sampling problem.
#'
#' @concept simulate
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

#' Subset to sampled or reference draws
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns a \code{redist_plans} object, with only rows corresponding to
#' simulated (or reference) draws remaining.
#'
#' @concept analyze
#' @export
subset_sampled = function(plans) {
    n_ref = get_n_ref(plans)
    dplyr::filter(plans, as.integer(.data$draw) > n_ref)
}

#' @rdname subset_sampled
#' @export
subset_ref = function(plans) {
    n_ref = get_n_ref(plans)
    dplyr::filter(plans, as.integer(.data$draw) <= n_ref)
}



# generics ----------------------------------------------------------------


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
    if ("district" %in% colnames(y)) {
        distrs = table(as.integer(y$district))
        n_draws = length(draws_left)
        n_distr = max(plans_m[,1])
        if (!all.equal(range(distrs), rep(n_draws, 2)) || length(distrs) != n_distr)
            warning("Some districts may have been dropped. ",
                    "This will prevent summary statistics from working correctly.\n",
                    "To avoid this message, coerce using `as_tibble`.")
    }

    if (length(draws_left) != ncol(plans_m)) {
        attr(y, "wgt") = attr(y, "wgt")[draws_left]
        y = set_plan_matrix(y, plans_m[,draws_left, drop=F])
    }

    y
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
    n_ref = get_n_ref(x)
    n_samp = ncol(plans_m) - n_ref

    cat(cli::pluralize("{n_samp} sampled plan{?s} "))
    if (n_ref > 0)
        cat(cli::pluralize("and {n_ref} reference plan{?s} "))
    if (ncol(plans_m) == 0) return(invisible(x))

    cli::cat_line("with ", max(plans_m[,1]), " districts from a ",
        nrow(plans_m), "-unit map,\n  drawn using ",
        c(mcmc="Markov chain Monte Carlo",
          smc="Sequential Monte Carlo",
          mergesplit="Merge-split Markov chain Monte Carlo")[attr(x, "algorithm")])

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
#' @param ... passed on to the underlying function
#' @param type the name of the plotting function to use. Will have \code{plot_},
#' prepended to it; e.g., use \code{type="plan"} to call \code{plot_plan()}.
#'
#' @concept visualize
#' @export
plot.redist_plans = function(x, ..., type="hist") {
    get(paste0("plot_", type))(x, ...)
}

#' Plot a histogram of a summary statistic
#'
#' Plots a histogram of a statistic of a \code{\link{redist_plans}} object,
#' with a reference line for each reference plan, if applicable.
#'
#' @param x the \code{redist_plans} object.
#' @param qty \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the statistic.
#' @param bins the number of bins to use in the histogram
#' @param ... passed on to \code{\link[ggplot2]{geom_histogram}}
#'
#' @concept visualize
#' @export
plot_hist = function(x, qty, bins=ceiling(sqrt(nrow(x))*1.5), ...) {
    p = ggplot(subset_sampled(x), aes({{ qty }})) +
        ggplot2::geom_histogram(..., bins=bins) +
        labs(y="Number of plans", color="Plan")
    if (get_n_ref(x) > 0)
        p = p + ggplot2::geom_vline(aes(xintercept={{ qty }}, color=.data$draw),
                                    data=subset_ref(x))
    p
}

#' @rdname plot_hist
#' @export
hist.redist_plans = function(x, qty, ...) {
    qty = rlang::enquo(qty)
    plot_hist(x, !!qty, ...)
}

#' Plot box and whisker plots for each district
#'
#' Plots a boxplot of a quantity of interest across districts, with districts
#' optionally sorted by this quantity. Adds reference points for each reference
#' plan, if applicable.
#'
#' @param x the \code{redist_plans} object.
#' @param qty \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the
#'   quantity of interest.
#' @param sort set to \code{"asc"} to sort districts in ascending order of
#'   \code{qty} (the default), \code{"desc"} for descending order, or
#'   \code{FALSE} or \code{"none"} for no sorting.
#' @param ... passed on to \code{\link[ggplot2]{geom_boxplot}}
#'
#' @concept visualize
#' @export
plot_distr_qtys = function(x, qty, sort="asc", ...) {
    if (isFALSE(sort) || sort == "none") {
        x = dplyr::group_by(x, .data$draw) %>%
            dplyr::mutate(.distr_no = as.factor(.data$district))
    } else {
        ord = if (sort == "asc") 1  else if (sort == "desc") -1 else
            stop("`sort` not recognized: ", sort)
        x = dplyr::group_by(x, .data$draw) %>%
            dplyr::mutate(.distr_no = as.factor(rank(ord * {{ qty }})))
    }

    p = ggplot(subset_sampled(x), aes(.data$.distr_no, {{ qty }})) +
        ggplot2::geom_boxplot(..., outlier.size=1)

    if (isFALSE(sort) || sort == "none")
        p = p + labs(x="District", color="Plan", shape="Plan")
    else
        p = p + labs(x="Ordered district", color="Plan", shape="Plan")

    if (get_n_ref(x) > 0)
        p = p + ggplot2::geom_point(aes(color=.data$draw, shape=.data$draw), size=2, data=subset_ref(x))

    p
}

#' Plot a district assignment
#'
#' @param x a \code{redist_plans} object.
#' @param draw the plan to plot. Will match the \code{draw} column of \code{x}.
#' @param geom the geometry to use
#'
#' @concept visualize
#' @export
plot_plan = function(x, draw, geom) {
    stopifnot(inherits(x, "redist_plans"))
    stopifnot()

    draw_idx = match(as.character(draw), levels(x$draw))
    distr_assign = get_plan_matrix(x)[,draw_idx]
    if (inherits(geom, "redist_map")) {
        distr_colors = as.factor(color_graph(get_graph(geom), distr_assign))
    } else {
        distr_colors = as.factor(distr_assign)
    }

    PAL = c("#6D9537", "#364B6F", "#E59A20", "#9A9BB9", "#2A4E45")
    sf::st_sf(geom) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(District = distr_colors) %>%
    ggplot(aes(fill=.data$District)) +
        ggplot2::geom_sf(size=0) +
        ggplot2::guides(fill=FALSE) +
        ggplot2::scale_fill_manual(values=PAL) +
        theme_void()
}
