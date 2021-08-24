##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################


# helper: TRUE if `x is constant within `grps`
is_const_num = function(x, grps) {
    all(tapply(x, grps, FUN=function(y) diff(range(y)) == 0))
}

#' Summary plots for \code{\\link{redist_plans}}
#'
#' If no arguments are passed, defaults to plotting the sampling weights for
#' the \code{\link{redist_plans}} object. If no weights exist, plots district
#' populations.
#'
#' @param x the \code{redist_plans} object.
#' @param ... passed on to the underlying function
#' @param type the name of the plotting function to use. Will have
#'   \code{redist.plot.}, prepended to it; e.g., use \code{type="plans"} to call
#'   \code{\link{redist.plot.plans}}.
#'
#' @concept plot
#' @export
plot.redist_plans = function(x, ..., type="distr_qtys") {
    if (rlang::dots_n(...) == 0) {
        wgts = get_plans_weights(subset_sampled(x))
        if (is.null(wgts))
            return(redist.plot.distr_qtys(x, total_pop, size=0.1))
        n = length(wgts)
        iqr = IQR(wgts)
        bins = min(max(round(diff(range(wgts)) / (2 * iqr / n^(1/3))), 3), 100)
        if (iqr == 0) bins = 3

        ggplot(NULL, aes(x=wgts)) +
            geom_histogram(bins=bins) +
            ggplot2::scale_x_continuous(name="Weights", trans="log10") +
            ggplot2::labs(y=NULL, title="Plan weights")
    } else {
        get(paste0("redist.plot.", type))(x, ...)
    }
}

#' Plot a histogram of a summary statistic
#'
#' Plots a histogram of a statistic of a \code{\link{redist_plans}} object,
#' with a reference line for each reference plan, if applicable.
#'
#' @param plans the \code{redist_plans} object.
#' @param qty \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the statistic.
#' @param bins the number of bins to use in the histogram. Defaults to Freedman-Diaconis rule.
#' @param ... passed on to \code{\link[ggplot2]{geom_histogram}}
#'
#' @returns A ggplot
#'
#' @examples
#' library(dplyr)
#' data(iowa)
#'
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05)
#' plans = redist_smc(iowa, nsims=100, silent=TRUE)
#' group_by(plans, draw) %>%
#'     summarize(pop_dev = max(abs(total_pop / mean(total_pop) - 1))) %>%
#'     redist.plot.hist(pop_dev)
#'
#' @concept plot
#' @export
redist.plot.hist = function(plans, qty, bins=NULL, ...) {
    stopifnot(inherits(plans, "redist_plans"))
    if (missing(qty))
        stop("Must provide a quantity to make the histogram from.", call.=FALSE)

    val = rlang::eval_tidy(rlang::enquo(qty), plans)
    rg = diff(range(val, na.rm=T))
    is_int = isTRUE(all.equal(as.integer(val), val)) && rg <= 100
    if (is.null(bins)) {
        if (is_int) {
            bins = 3*rg + 1
        } else { # Freedman-Diaconis
            n = length(val)
            if (is_const_num(val, plans$draw)) {
                n = n / nrow(get_plans_matrix(plans))
            }
            iqr = IQR(val, na.rm=T)
            if (iqr > 0)
                bins = max(round(rg / (2 * iqr / n^(1/3))), 3)
            else
                bins = 3
        }
    }

    percent = function(x) sprintf("%1.0f%%", 100*x)
    p = ggplot(subset_sampled(plans), aes({{ qty }})) +
        ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density*width)), ...,
                                boundary=0.5*is_int, bins=bins) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)),
                                    labels = percent) +
        labs(y="Fraction of plans")
    if (get_n_ref(plans) > 0)
        p = p + labs(color="Plan") +
            ggplot2::geom_vline(aes(xintercept={{ qty }}, color=.data$draw),
                                size=1.15, data=subset_ref(plans))
    p
}

#' @rdname redist.plot.hist
#' @param x \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the statistic.
#' @export
hist.redist_plans = function(x, qty, ...) {
    if (missing(qty))
        stop("Must provide a quantity to make the histogram from.", call.=FALSE)
    qty = rlang::enquo(qty)
    redist.plot.hist(x, !!qty, ...)
}

#' Scatter plot of plan summary statistics
#'
#' Makes a scatterplot of two quantities of interest across districts or plans.
#'
#' @param plans the \code{redist_plans} object.
#' @param x \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the
#'   quantity to plot on the horizontal axis.
#' @param y \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the
#'   quantity to plot on the vertical axis.
#' @param ... passed on to \code{\link[ggplot2:geom_point]{geom_point}}.
#' @param bigger if TRUE, make the point corresponding to the reference plan larger.
#'
#' @returns A ggplot
#'
#' @examples
#' library(dplyr)
#' data(iowa)
#'
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05, total_pop = pop)
#' plans = redist_smc(iowa, nsims=100, silent=TRUE)
#' plans %>%
#'     mutate(comp = distr_compactness(iowa)) %>%
#'     group_by(draw) %>%
#'     summarize(pop_dev = max(abs(total_pop / mean(total_pop) - 1)),
#'               comp = comp[1]) %>%
#'     redist.plot.scatter(pop_dev, comp)
#'
#' @concept plot
#' @export
redist.plot.scatter = function(plans, x, y, ..., bigger=TRUE) {
    stopifnot(inherits(plans, "redist_plans"))

    p = ggplot(subset_sampled(plans), aes(x={{ x }}, y={{ y }})) +
        ggplot2::geom_point(...)
    if (get_n_ref(plans) > 0) {
        p = p + labs(color="Plan")
        if (bigger) {
            p = p + ggplot2::geom_point(aes(color=.data$draw), shape=15,
                                        size=3, data=subset_ref(plans))
        } else {
            p = p + ggplot2::geom_point(aes(color=.data$draw), shape=15,
                                        data=subset_ref(plans))
        }
    }

    p
}

#' Plot quantities by district
#'
#' Plots a boxplot of a quantity of interest across districts, with districts
#' optionally sorted by this quantity. Adds reference points for each reference
#' plan, if applicable.
#'
#' @param plans the \code{redist_plans} object.
#' @param qty \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the
#'   quantity of interest.
#' @param sort set to \code{"asc"} to sort districts in ascending order of
#'   \code{qty} (the default), \code{"desc"} for descending order, or
#'   \code{FALSE} or \code{"none"} for no sorting.
#' @param geom the geom to use in plotting the simulated districts: either
#'   \code{"jitter"} or \code{"boxplot"}
#' @param color_thresh if a number, the threshold to use in coloring the points.
#'   Plans with quantities of interest above the threshold will be colored
#'   differently than plans below the threshold.
#' @param ... passed on to \code{\link[ggplot2]{geom_boxplot}}
#'
#' @returns A ggplot
#'
#' @examples
#' library(dplyr)
#' data(iowa)
#'
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05, total_pop = pop)
#' plans = redist_smc(iowa, nsims=100, silent=TRUE)
#' plans %>%
#'     mutate(pct_dem = group_frac(iowa, dem_08, tot_08)) %>%
#'     redist.plot.distr_qtys(pct_dem)
#'
#' @concept plot
#' @export
redist.plot.distr_qtys = function(plans, qty, sort="asc", geom="jitter",
                                  color_thresh=NULL, ...) {
    stopifnot(inherits(plans, "redist_plans"))

    if (isFALSE(sort) || sort == "none") {
        plans = dplyr::group_by(plans, .data$draw) %>%
            dplyr::mutate(.distr_no = as.factor(.data$district))
    } else {
        ord = if (sort == "asc") 1  else if (sort == "desc") -1 else
            stop("`sort` not recognized: ", sort)
        plans = dplyr::group_by(plans, .data$draw) %>%
            dplyr::mutate(.distr_no = as.factor(rank(ord * {{ qty }},
                                                     ties.method="first")))
    }

    val = eval_tidy(enquo(qty), plans)
    if (is_const_num(val, plans$draw)) {
        warning("`", rlang::as_label(enquo(qty)),
                "` is constant across districts. ",
                "Consider using `hist()` instead.", call.=FALSE)
    }

    if (is.null(color_thresh)) {
        p = ggplot(subset_sampled(plans), aes(.data$.distr_no, {{ qty }}))
    } else {
        stopifnot(is.numeric(color_thresh))
        p = ggplot(subset_sampled(plans), aes(.data$.distr_no, {{ qty }},
                                              color = {{ qty }} >= color_thresh)) +
            ggplot2::guides(color="none")
    }

    if (geom == "jitter") {
        p = p + ggplot2::geom_jitter(...)
    } else if (geom == "boxplot") {
        p = p + ggplot2::geom_boxplot(..., outlier.size=1)
    } else {
        stop('`geom` must be either "jitter" or "boxplot"')
    }

    if (isFALSE(sort) || sort == "none")
        p = p + labs(x="District")
    else
        p = p + labs(x="Ordered district")

    if (get_n_ref(plans) > 0) {
        if (is.null(color_thresh)) {
            p = p + labs(color="Plan", shape="Plan")
            if (geom == "jitter") {
                p = p + ggplot2::geom_segment(aes(as.integer(.data$.distr_no)-0.5,
                                                  xend=as.integer(.data$.distr_no)+0.5,
                                                  yend={{ qty }},
                                                  color=.data$draw),
                             data=subset_ref(plans), size=1.2)
            } else {
                p = p + ggplot2::geom_point(aes(color=.data$draw), shape=15,
                                            size=2, data=subset_ref(plans))
            }
        } else {
            if (geom == "jitter") {
                p = p + labs(lty="Plan") +
                    ggplot2::geom_segment(aes(as.integer(.data$.distr_no)-0.5,
                                              xend=as.integer(.data$.distr_no)+0.5,
                                              yend={{ qty }},
                                              lty=.data$draw),
                                          data=subset_ref(plans),
                                          size=1.2, color="black")
            } else {
                p = p + labs(shape="Plan") +
                    ggplot2::geom_point(aes(shape=.data$draw), color="black",
                                        size=2, data=subset_ref(plans))
            }
        }
    }

    p
}

#' Plot a district assignment
#'
#' @param plans a \code{redist_plans} object.
#' @param draws the plan(s) to plot. Will match the \code{draw} column of \code{x}.
#' @param geom the \code{redist_map} geometry to use
#' @param qty the quantity to plot. Defaults to the district assignment.
#' @param interactive if \code{TRUE}, show an interactive map in the viewer
#'   rather than a static map. Only uses the first element of \code{draws}
#' @param ... additional arguments passed to the plotting functions.
#'
#' @returns A ggplot
#'
#' @examples
#' library(dplyr)
#' data(iowa)
#'
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05, total_pop = pop)
#' plans = redist_smc(iowa, nsims=100, silent=TRUE)
#' redist.plot.plans(plans, c(1, 2, 3, 4), iowa)
#'
#' @concept plot
#' @export
redist.plot.plans = function(plans, draws, geom, qty=NULL, interactive=FALSE, ...) {
    stopifnot(inherits(plans, "redist_plans"))
    m = get_plans_matrix(plans)
    stopifnot(nrow(geom) == nrow(m))

    if (interactive) {
        if (length(draws) > 1)
            warning("Only first of `draws` plotted when `interactive=TRUE`")

        draw = draws[1]
        draw_idx = match(as.character(draw), levels(plans$draw))
        qty = eval_tidy(enquo(qty), plans[plans$draw == as.character(draw), ])
        if (is.null(qty)) {
            qty = as.factor(m[, draw_idx])
        } else {
            qty = qty[m[, draw_idx]]
        }
        return(redist.plot.interactive(geom, fill=qty, ...))
    }

    plot_single = function(draw) {
        draw_idx = match(as.character(draw), levels(plans$draw))
        lab = rlang::quo_text(enquo(qty))
        title = if (suppressWarnings(is.na(as.numeric(draw)))) draw else paste0("Plan #", draw)

        qty = eval_tidy(enquo(qty), plans[plans$draw == as.character(draw), ])
        if (is.null(qty)) {
            qty = as.factor(m[, draw_idx])
        } else {
            qty = qty[m[, draw_idx]]
        }

        redist.plot.map(geom, fill=qty, fill_label=lab, ...) +
            ggplot2::labs(title=title)
    }

    if (length(draws) == 1) {
        plot_single(draws)
    } else {
        plots = lapply(draws, plot_single)
        patchwork::wrap_plots(plots)
    }
}

utils::globalVariables(c("density", "width", "total_pop"))
