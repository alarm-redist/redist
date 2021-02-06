##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################

# redist_map functions ----


#' Check that a \code{redist_map} object is contiguous
#'
#' @param x the object
#'
#' @return \code{TRUE} if contiguous.
#' @concept prepare
#'
#' @export
is_contiguous = function(x) {
    stopifnot(inherits(x, "redist_map"))
    all(contiguity(get_graph(x), rep(1, nrow(x))) == 1)
}


## merge helpers
# checks if is a proportion/pct
is_prop = function(x) is.double(x) && min(x, na.rm=T) >= 0 && max(x, na.rm=T) <= 1
# checks if is not a  proportion/pct
is_nonprop = function(x) is.numeric(x) && (min(x, na.rm=T) < 0 || max(x, na.rm=T) > 1)
# checks if x is constant within levels of rel
is_const_rel = function(rel) {
    function(x) {
        !is.numeric(x) && all(tapply(x, rel, FUN=function(y) length(unique(y))) == 1)
    }
}

#' Merge map units
#'
#' In performing a county-level or cores-based analysis it is often necessary to
#' merge several units together into a larger unit.  This function performs this
#' operation, modifying the adjacency graph as needed and attempting to properly
#' aggregate other data columns.
#'
#' @param .data a \code{\link{redist_map}} object
#' @param key \code{\link[tidyr:tidyr_tidy_select]{<tidy-select>}} the column to merge by
#' @param by_existing if an existing assignment is present, whether to also group by it
#' @param drop_geom whether to drop the geometry column. Recommended, as
#'   otherwise a costly geometric merge is required.
#'
#' @returns A merged \code{\link{redist_map}} object
#'
#' @concept prepare
#' @export
merge_by = function(.data, key, by_existing=TRUE, drop_geom=TRUE) {
    .data = as_redist_map(.data)

    key_val = rlang::eval_tidy(rlang::enquo(key), .data)
    pop_col = attr(.data, "pop_col")

    if (drop_geom)
        .data = sf::st_drop_geometry(.data)

    col = attr(.data, "existing_col")
    unique_chr = function(x) paste(unique(x), collapse="~")
    if (!is.null(col) && by_existing) {
        dplyr::group_by(.data, dplyr::across(dplyr::all_of(col)), {{ key }}) %>%
            dplyr::summarize(dplyr::across(where(is_prop),
                                           ~ weighted.mean(., w=.data[[pop_col]], na.rm=T)),
                             dplyr::across(where(is.numeric), sum, na.rm=T),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is.character), unique_chr))
    } else {
        dplyr::group_by(.data, {{ key }}) %>%
            dplyr::summarize(dplyr::across(where(is_prop),
                                           ~ weighted.mean(., w=.data[[pop_col]], na.rm=T)),
                             dplyr::across(where(is.numeric), sum, na.rm=T),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is.character), unique_chr)) %>%
            `attr<-`("existing_col", NULL)
    }
}

#' @rdname redist.identify.cores
#'
#' @param .data a \code{\link{redist_map}} object
#' @param within the core is defined to be at least this number of steps within
#'   district boundaries
#'
#' @concept prepare
#' @export
make_cores = function(.data=get0(".", parent.frame()), within=1, focus=NULL) {
    if (!inherits(.data, "redist_map"))
        stop("Must provide `.data` if not called within a pipe")

    redist.identify.cores(get_graph(.data), get_existing(.data), within, focus, simplify=TRUE)
}



# redist_plans functions ----


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
#' @concept analyze
#' @export
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

    attr(plans, "wgt") = c(0, get_plan_weights(plans))

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
#' @concept analyze
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

#' @rdname redist.compactness
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param ... passed on to \code{redist.compactness}
#'
#' @concept analyze
#' @export
distr_compactness = function(map, measure="EdgesRemoved", .data=get0(".", parent.frame()), ...) {
    if (!inherits(.data, "redist_plans"))
        stop("Must provide `.data` if not called within a pipe")

    redist.compactness(map, get_plan_matrix(.data), measure,
                       population=map[[attr(map, "pop_col")]],
                       adjacency=get_graph(map), ...)[[measure]]
}

#' @rdname redist.group.percent
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
group_frac = function(map, group_pop, full_pop=map[[attr(map, "pop_col")]],
                          .data=get0(".", parent.frame())) {
    if (!inherits(.data, "redist_plans"))
        stop("Must provide `.data` if not called within a pipe")

    group_pop = rlang::eval_tidy(rlang::enquo(group_pop), map)
    full_pop = rlang::eval_tidy(rlang::enquo(full_pop), map)
    as.numeric(redist.group.percent(get_plan_matrix(.data), group_pop, full_pop))
}


#' Confidence Intervals for Importance Sampling Estimates
#'
#' Builds a confidence interval for the mean of a vector of interest,
#' given importance sampling weights.
#'
#' @param x \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the vector to
#'   build importance sampling confidence intervals for.
#' @param conf The confidence level for the intervals.
#' @param .data a \code{\link{redist_plans}} object
#'
#' @returns A tibble with three columns: \code{X}, \code{X_lower}, and
#'   \code{X_upper}, where \code{X} is the name of the vector of interest,
#'   containing the mean and confidence interval. When used inside
#'   \code{\link[dplyr:summarise]{summarize()}} this will create three columns in the
#'   output data.
#'
#' @concept analyze
#' @export
imp_confint = function(x, conf=0.95, .data=get0(".", parent.frame())) {
    if (!inherits(.data, "redist_plans"))
        stop("Must provide `.data` if not called within a pipe")

    y = rlang::eval_tidy(rlang::enquo(x), .data)
    ci = redist.smc_is_ci(y, get_plan_weights(.data), conf)

    tibble::tibble("{{ x }}" := mean(y),
                   "{{ x }}_lower" := ci[1],
                   "{{ x }}_upper" := ci[2])
}
