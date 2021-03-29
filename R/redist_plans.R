##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################


# constructors and reconstructors -----------------------------------------


# plans has n_precinct columns and n_sims rows
# map is a redist_map
# algorithm is one of "smc" or "mcmc"
# wgt is the weights before any resampling or truncation
# ... will depend on the algorithm
new_redist_plans = function(plans, map, algorithm, wgt, resampled=TRUE, ...) {
    n_sims = ncol(plans)
    stopifnot(n_sims >= 1)

    n_prec = nrow(plans)
    ndists = attr(map, "ndists")

    prec_pop = map[[attr(map, "pop_col")]]
    distr_pop = pop_tally(plans, prec_pop, ndists)

    attr_names = c("redist_attr", "plans", "algorithm", "wgt", "resampled",
                   "merge_idx", "prec_pop", names(list(...)))

    structure(tibble(draw = rep(as.factor(1:n_sims), each=ndists),
                             district = rep(1:ndists, n_sims),
                             total_pop = as.numeric(distr_pop)),
              plans=plans, algorithm=algorithm, wgt=wgt,
              resampled=resampled, merge_idx=attr(map, "merge_idx"),
              prec_pop=prec_pop, redist_attr=attr_names, ...,
              class=c("redist_plans", "tbl_df", "tbl", "data.frame"))
}

validate_redist_plans = function(x) {
    stopifnot(names(x)[1] == "draw")
    stopifnot(is.factor(x$draw))

    plan_m = attr(x, "plans")
    stopifnot(!is.null(plan_m))

    min_distr = colmin(plan_m)
    max_distr = colmax(plan_m)
    stopifnot(all(min_distr == 1))
    stopifnot(all(diff(max_distr) == 0))

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

#' A set of redistricting plans
#'
#' A \code{redist_plans} object is essentially a data frame of summary
#' information on each district and each plan, along with the matrix of district
#' assignments and information about the simulation process used to generate the
#' plans.
#'
#' The first two columns of the data frame will be \code{draw}, a factor indexing
#' the simulation draw, and \code{district}, an integer indexing the districts
#' within a plan. The data frame will therefore have \code{n_sims*ndists} rows.
#' As a data frame, the usual \code{\link{dplyr}} methods will work.
#'
#' Other useful methods for \code{redist_plans} objects:
#' * \code{\link{add_reference}}
#' * \code{\link{subset_sampled}}
#' * \code{\link{subset_ref}}
#' * \code{\link{pullback}}
#' * \code{\link{number_by}}
#' * \code{\link{match_numbers}}
#' * \code{\link{prec_assignment}}
#' * \code{\link{get_plan_matrix}}
#' * \code{\link{get_plan_weights}}
#' * \code{\link{get_sampling_info}}
#' * \code{\link{as.matrix.redist_plans}}
#' * \code{\link{plot.redist_plans}}
#'
#' @param plans a matrix with \code{n_precinct} columns and \code{n_sims} rows,
#'   or a single vector of precinct assignments.
#' @param map a \code{\link{redist_map}} object
#' @param algorithm the algorithm used to generate the plans (usually "smc" or "mcmc")
#' @param wgt the weights to use, if any.
#' @param ... Other named attributes to set
#'
#' @returns a new \code{redist_plans} object.
#'
#' @examples
#' data(iowa)
#'
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05)
#' rsg_plan = redist.rsg(iowa$adj, iowa$pop, ndists=4, pop_tol=0.05)$plan
#' redist_plans(rsg_plan, iowa, "rsg")
#'
#' @md
#' @concept analyze
#' @export
redist_plans = function(plans, map, algorithm, wgt=NULL, ...) {
    if (is.numeric(plans) && length(plans) == nrow(map)) {
        plans = matrix(as.integer(plans), ncol=1)
    }
    stopifnot(is.matrix(plans))
    stopifnot(nrow(plans) == nrow(map))
    stopifnot(inherits(map, "redist_map"))

    if (min(plans) == 0L) plans = plans + 1L

    obj = new_redist_plans(plans, map, algorithm, wgt=wgt,
                           resampled=FALSE, ...)
    validate_redist_plans(obj)
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
#' @rdname get_plan_matrix
#' @method as.matrix redist_plans
#' @export
as.matrix.redist_plans = function(x, ...) get_plan_matrix(x)

# internal -- no check performed!
set_plan_matrix = function(x, mat) {
    attr(x, "plans") = mat
    x
}

#' Extract the sampling weights from a redistricting simulation.
#'
#' May be \code{NULL} if no weights exist (MCMC or optimization methods).
#'
#' @param plans,object the \code{redist_plans} object
#'
#' @returns A numeric vector of weights, with an additional attribute
#'   \code{resampled} indicating whether the plans have been resampled according
#'   to these weights.
#'
#' @concept analyze
#' @export
get_plan_weights = function(plans) {
    stopifnot(inherits(plans, "redist_plans"))
    wgt = attr(plans, "wgt")
    if (!is.null(wgt))
        attr(wgt, "resampled") = attr(plans, "resampled")
    wgt
}

#' @rdname get_plan_weights
#' @param ... Ignored.
#' @importFrom stats weights
#' @method weights redist_plans
#' @export
weights.redist_plans = function(object, ...) {
    get_plan_weights(object)
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
        ndists = max(plans_m[,1])
        if (!all.equal(range(distrs), rep(n_draws, 2)) || length(distrs) != ndists)
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

    if (n_samp == 1) {
        cat("1 sampled plan ")
    } else {
        cat(n_samp, "sampled plans ")
    }

    if (n_ref == 1) {
        cat("and 1 reference plan ")
    } else if (n_ref > 1) {
        cat("and", n_ref, "reference plans ")
    }
    if (ncol(plans_m) == 0) return(invisible(x))

    alg_name = c(mcmc="Markov chain Monte Carlo",
                 smc="Sequential Monte Carlo",
                 mergesplit="Merge-split Markov chain Monte Carlo",
                 rsg="random seed-and-grow",
                 crsg="compact random seed-and-grow",
                 enumpart="Enumpart",
                 shortburst="short bursts")[attr(x, "algorithm")]
    if (is.na(alg_name)) alg_name = "an unknown algorithm"

    cat("with ", max(plans_m[,1]), " districts from a ",
        nrow(plans_m), "-unit map,\n  drawn using ", alg_name, "\n", sep="")

    merge_idx = attr(x, "merge_idx")
    if (!is.null(merge_idx))
        cat("Merged from another map with reindexing:",
            utils::capture.output(str(merge_idx, vec.len=2)), "\n", sep="")

    if (!is.null(attr(x, "wgt"))) {
        if (attr(x, "resampled"))
            cat("With plans resampled from weights\n")
        else
            cat("With plans not resampled from weights\n")
    }

    cat("Plans matrix:", utils::capture.output(str(plans_m, give.attr=F)),
        "\n", sep="")

    utils::getS3method("print", "tbl")(x)

    invisible(x)
}

#' Summary plots for \code{redist_plans}
#'
#' @param x the \code{redist_plans} object.
#' @param ... passed on to the underlying function
#' @param type the name of the plotting function to use. Will have
#'   \code{redist.plot.}, prepended to it; e.g., use \code{type="plans"} to call
#'   \code{\link{redist.plot.plans}}.
#'
#' @concept plot
#' @export
plot.redist_plans = function(x, ..., type="hist") {
    if (rlang::dots_n(...) == 0) {
        ggplot(tibble(wgt=get_plan_weights(x)), aes(x=.data$wgt)) +
            geom_histogram() +
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

    if (is.null(bins)) {
        val = rlang::eval_tidy(rlang::enquo(qty), plans)
        n = length(val)
        iqr = IQR(val)
        if (iqr > 0)
            bins = max(round(diff(range(val)) / (2 * iqr / n^(1/3))), 3)
        else
            bins = 3
    }

    p = ggplot(subset_sampled(plans), aes({{ qty }})) +
        ggplot2::geom_histogram(..., bins=bins) +
        labs(y="Number of plans")
    if (get_n_ref(plans) > 0)
        p = p + labs(color="Plan") +
            ggplot2::geom_vline(aes(xintercept={{ qty }}, color=.data$draw),
                                data=subset_ref(plans))
    p
}

#' @rdname redist.plot.hist
#' @param x \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the statistic.
#' @export
hist.redist_plans = function(x, qty, ...) {
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
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05)
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
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05)
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
            dplyr::mutate(.distr_no = as.factor(rank(ord * {{ qty }})))
    }

    if (is.null(color_thresh)) {
        p = ggplot(subset_sampled(plans), aes(.data$.distr_no, {{ qty }}))
    } else {
        stopifnot(is.numeric(color_thresh))
        p = ggplot(subset_sampled(plans), aes(.data$.distr_no, {{ qty }},
                                              color = {{ qty }} >= color_thresh)) +
            ggplot2::guides(color=FALSE)
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
#' @param geom the geometry to use
#'
#' @returns A ggplot
#'
#' @examples
#' library(dplyr)
#' data(iowa)
#'
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05)
#' plans = redist_smc(iowa, nsims=100, silent=TRUE)
#' redist.plot.plans(plans, c(1, 2, 3, 4), iowa)
#'
#' @concept plot
#' @export
redist.plot.plans = function(plans, draws, geom) {
    stopifnot(inherits(plans, "redist_plans"))

    plot_single = function(draw) {
        draw_idx = match(as.character(draw), levels(plans$draw))
        distr_assign = get_plan_matrix(plans)[, draw_idx]
        if (inherits(geom, "redist_map")) {
            distr_colors = as.factor(color_graph(get_adj(geom), distr_assign))
        } else {
            distr_colors = as.factor(distr_assign)
        }

        PAL = c("#6D9537", "#364B6F", "#E59A20", "#9A9BB9", "#2A4E45")
        title = if (suppressWarnings(is.na(as.numeric(draw)))) draw else paste0("Plan #", draw)
        sf::st_sf(geom) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(District = distr_colors) %>%
        ggplot(aes(fill=.data$District)) +
            ggplot2::geom_sf(size=0) +
            ggplot2::guides(fill=FALSE) +
            ggplot2::scale_fill_manual(values=PAL) +
            ggplot2::labs(title=title) +
            theme_void()
    }

    if (length(draws) == 1) {
        plot_single(draws)
    } else {
        plots = lapply(draws, plot_single)
        patchwork::wrap_plots(plots)
    }
}
