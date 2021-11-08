##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################

# tidy accessor/helper functons ----

#' Helper function to get current plans/map objects
#' Traverses call stack to find plans object passed to dplyr verbs
#' @noRd
get_cur_df = function(dplyr_funcs) {
    calls = sys.calls()
    frames = sys.frames()
    for (i in rev(seq_along(calls))) {
        call = calls[[i]]
        frame = frames[[i]]
        if (is.null(rlang::call_name(call))) next
        if (any(vapply(dplyr_funcs,
                       function(x) identical(x, rlang::call_fn(call, frame)),
                       logical(1)))) {
            return(rlang::env_get(frame, ".data"))
        }
    }
    return(NULL)
}

#' Helper function to get current map object
#' @noRd
cur_map = function(verbs=c("mutate", "summarize", "merge_by",
                           "filter", "arrange", "transmute")) {
    get_cur_df(list(mutate=mutate.redist_map,
                    transmute=transmute.redist_map,
                    summarize=summarise.redist_map,
                    merge_by=merge_by,
                    filter=filter.redist_map,
                    arrange=arrange.redist_map)[verbs])
}

#' Helper function to get current plans object
#' @noRd
cur_plans = function(verbs=c("mutate", "summarize", "filter",
                             "arrange", "transmute")) {
    get_cur_df(list(mutate=mutate.redist_plans,
                    transmute=transmute.redist_plans,
                    summarize=summarise.redist_plans,
                    filter=filter.redist_plans,
                    arrange=arrange.redist_plans)[verbs])
}
#

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
    if (!inherits(x, "redist_map")) cli_abort("{.arg x} must be a {.cls redist_map}")
    all(contiguity(get_adj(x), rep(1, nrow(x))) == 1)
}


## merge helpers
# checks if is a proportion/pct
is_prop = function(x) is.double(x) && min(x, na.rm=TRUE) >= 0 && max(x, na.rm=TRUE) <= 1
# checks if is not a  proportion/pct
is_nonprop = function(x) is.numeric(x) && (min(x, na.rm=TRUE) < 0 || max(x, na.rm=TRUE) > 1)
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
#' @param ... \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} the column(s) to merge by
#' @param by_existing if an existing assignment is present, whether to also group by it
#' @param drop_geom whether to drop the geometry column. Recommended, as
#'   otherwise a costly geometric merge is required.
#' @param collapse_chr if \code{TRUE}, preserve character columns by collapsing
#'   their values. For example, a county name column in Iowa might be merged and
#'   have entries such as "Cedar~Clinton~Des Moines". Set to \code{FALSE} to
#'   drop character columns instead.
#'
#' @returns A merged \code{\link{redist_map}} object
#'
#' @concept prepare
#' @export
merge_by = function(.data, ..., by_existing=TRUE, drop_geom=TRUE, collapse_chr=TRUE) {
    .data = as_redist_map(.data)

    dots = rlang::enquos(...)
    key_val = rlang::eval_tidy(dots[[1]], .data)
    pop_col = attr(.data, "pop_col")

    if (drop_geom && is(.data, "sf"))
        .data = sf::st_drop_geometry(.data)

    col = attr(.data, "existing_col")
    unique_chr = function(x) paste(unique(x), collapse="~")
    is_col_chr = if (collapse_chr) is.character else (function(x) FALSE)
    if (!is.null(col) && by_existing) {
        dplyr::group_by(.data, dplyr::across(dplyr::all_of(col)), !!!dots) %>%
            dplyr::summarize(dplyr::across(where(is_prop),
                                           ~ weighted.mean(., w=.data[[pop_col]], na.rm=TRUE)),
                             dplyr::across(where(is.numeric), sum, na.rm=TRUE),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is_col_chr), unique_chr),
                             .groups="drop")
    } else {
        dplyr::group_by(.data, !!!dots) %>%
            dplyr::summarize(dplyr::across(where(is_prop),
                                           ~ weighted.mean(., w=.data[[pop_col]], na.rm=TRUE)),
                             dplyr::across(where(is.numeric), sum, na.rm=TRUE),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is_col_chr), unique_chr),
                             .groups="drop") %>%
            `attr<-`("existing_col", NULL)
    }
}

#' @rdname redist.identify.cores
#' @order 1
#'
#' @param .data a \code{\link{redist_map}} object
#' @param boundary the core is defined to be at least this number of steps within
#'   district boundaries
#'
#' @concept prepare
#' @export
make_cores = function(.data=cur_map(), boundary=1, focus=NULL) {
    if (is.null(.data))
        cli_abort("Must provide {.arg .data} if not called within a {.pkg dplyr} verb")
    if (!inherits(.data, "redist_map")) cli_abort("{.arg .data} must be a {.cls redist_map}")

    existing = get_existing(.data)
    if (is.null(existing))
        cli_abort(c("No existing plan found from which to compute cores.",
                    ">"="Add one using the {.arg existing_plan} argument to {.fun redist_map}"))

    redist.identify.cores(adj=get_adj(.data),
                          plan=as.integer(as.factor(existing)),
                          boundary=boundary, focus=focus, simplify=TRUE)
}


# redist_plans functions ----


#' Add a reference plan to a set of plans
#'
#' This function facilitates comparing an existing (i.e., non-simulated)
#' redistricting plan to a set of simulated plans.
#'
#' @param plans a \code{redist_plans} object
#' @param ref_plan an integer vector containing the reference plan. It will be
#' renumbered to 1..\code{ndists}.
#' @param name a human-readable name for the reference plan. Defaults to the
#' name of \code{ref_plan}.
#'
#' @returns a modified \code{redist_plans} object containing the reference plan
#' @concept analyze
#' @export
add_reference = function(plans, ref_plan, name=NULL) {
    if (!inherits(plans, "redist_plans")) cli_abort("{.arg plans} must be a {.cls redist_plans}")
    if (isTRUE(attr(plans, "partial")))
        cli_abort("Reference plans not supported for partial plans objects.")

    plan_m = get_plans_matrix(plans)
    if (!is.numeric(ref_plan)) cli_abort("{.arg ref_plan} must be numeric")
    if (length(ref_plan) != nrow(plan_m))
        cli_abort("{.arg ref_plan} must have the same number of precincts as {.arg plans}")

    if (is.null(name)) {
        ref_str = deparse(substitute(ref_plan))
        if (stringr::str_detect(ref_str, stringr::fixed("$")))
            name = strsplit(ref_str, "$", fixed=TRUE)[[1]][2]
        else
            name = ref_str
    } else {
        if (!is.character(name)) cli_abort("{.arg ref_plan} must be a {.cls chr}")
    }

    ref_plan = as.integer(as.factor(ref_plan))
    ndists = max(ref_plan)
    if (ndists != max(plan_m[, 1]))
        cli_abort("{.arg ref_plan} has a different number of districts than {.arg plans}")

    # first the matrix
    plan_m = cbind(ref_plan, plan_m)
    colnames(plan_m)[1] = name

    # then the dataframe
    prec_pop = attr(plans, "prec_pop")
    if (!is.null(prec_pop))
        distr_pop = pop_tally(matrix(ref_plan, ncol=1), prec_pop, ndists)
    else
        distr_pop = rep(NA_real_, ndists)

    if (name %in% levels(plans$draw)) cli_abort("Reference plan name already exists")
    fct_levels = c(name, levels(plans$draw))
    new_draw = rep(factor(fct_levels, levels=fct_levels), each=ndists)
    x = dplyr::bind_rows(
            tibble(district = 1:ndists,
                           total_pop = as.numeric(distr_pop)),
            plans[, -match("draw", names(plans))]
        ) %>%
        dplyr::mutate(draw = new_draw, .before="district")

    exist_wgts = get_plans_weights(plans)
    if (!is.null(exist_wgts))
        attr(plans, "wgt") = c(0, exist_wgts)

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
#' @param map optionally, a \code{redist_map} object, which will be used to set the new population vector
#'
#' @returns a new, re-indexed, \code{redist_plans} object
#'
#' @concept analyze
#' @export
pullback = function(plans, map=NULL) {
    if (!inherits(plans, "redist_plans")) cli_abort("{.arg plans} must be a {.cls redist_plans}")

    merge_idx = attr(plans, "merge_idx")
    if (is.null(merge_idx)) {
        cli_warn("No merged indexing found.")
        return(plans)
    }

    attr(plans, "merge_idx") = NULL
    if (inherits(NULL, "redist_map")) {
        attr(plans, "prec_pop") = map[[attr(map, "pop")]]
    } else {
        attr(plans, "prec_pop") = NULL
    }

    set_plan_matrix(plans, get_plans_matrix(plans)[merge_idx,])
}


# helper function for match_numbers
find_numbering = function(plan, ref, pop) {
    joint = plan_joint(ref, plan, pop)
    tot_pop = sum(pop)

    renumb = solve_hungarian(1 - joint / tot_pop)[, 2]

    list(renumb = renumb,
         shared = sum(diag(joint[, renumb])) / tot_pop)
}

#' Renumber districts to match an existing plan
#'
#' District numbers in simulated plans are by and large random.  This
#' function attempts to renumber the districts across all simulated plans to
#' match the numbers in a provided plan, using the Hungarian algorithm.
#'
#' @param data a \code{redist_plans} object
#' @param plan a character vector giving the name of the plan to match to (e.g.,
#'   for a reference plan), or an integer vector containing the plan itself.
#' @param col the name of a new column to store the vector of population overlap
#'   with the reference plan: the fraction of the total population who are in
#'   the same district under each plan and the reference plan. Set to
#'   \code{NULL} if no column should be created.
#'   renumbering options in any plan.
#'
#' @returns a modified \code{redist_plans} object. New district numbers will be
#' stored as an ordered factor variable in the \code{district} column. The
#' district numbers in the plan matrix will match the levels of this factor.
#'
#' @examples
#' data(iowa)
#'
#' iowa_map = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05)
#' plans = redist_smc(iowa_map, 100, silent=TRUE)
#' match_numbers(plans, "cd_2010")
#'
#' @concept analyze
#' @export
match_numbers = function(data, plan, col="pop_overlap") {
    if (!inherits(data, "redist_plans")) cli_abort("{.arg data} must be a {.cls redist_plans}")
    if (!"district" %in% colnames(data)) cli_abort("Missing {.field district} colun in {.arg data}")

    plan_mat = get_plans_matrix(data)
    if (is.character(plan)) plan = plan_mat[,plan]
    plan = factor(plan, ordered=TRUE)
    ndists = length(levels(plan))
    pop = attr(data, "prec_pop")

    if (is.null(pop)) cli_abort("{.field prec_pop} attribute in {.arg data} required.")
    if (max(plan_mat[,1]) != ndists)
        cli_abort("Can't match numbers on a subset of a {.cls redist_plans}")

    # compute renumbering and extract info
    best_renumb = apply(plan_mat, 2, find_numbering, as.integer(plan), pop)
    renumb = as.integer(vapply(best_renumb, function(x) x$renumb, integer(ndists)))

    if (!is.null(col))
        data[[col]] = as.numeric(vapply(best_renumb, function(x) rep(x$shared, ndists),
                                        numeric(ndists)))

    renumb_mat = renumber_matrix(plan_mat, renumb)
    colnames(renumb_mat) = colnames(plan_mat)
    data = set_plan_matrix(data, renumb_mat)
    data$district = factor(levels(plan)[renumb], levels(plan), ordered=TRUE)

    orig_groups = dplyr::group_vars(data)
    dplyr::group_by(data, .data$draw) %>%
        dplyr::arrange(.data$district, .by_group=TRUE) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(orig_groups)))
}

#' Renumber districts to match a quantity of interest
#'
#' District numbers in simulated plans are by and large random.  This
#' function will renumber the districts across all simulated plans in order
#' of a provided quantity of interest.
#'
#' @param data a \code{redist_plans} object
#' @param x \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the quantity of interest.
#' @param desc \code{TRUE} if district should be sorted in descending order.
#'
#' @returns a modified \code{redist_plans} object. New district numbers will be
#' stored as an ordered factor variable in the \code{district} column. The
#' district numbers in the plan matrix will match the levels of this factor.
#'
#' @concept analyze
#' @export
number_by = function(data, x, desc=FALSE) {
    if (!inherits(data, "redist_plans")) cli_abort("{.arg data} must be a {.cls redist_plans}")
    if (!"district" %in% colnames(data)) cli_abort("Missing {.field district} colun in {.arg data}")

    ord = 1 - 2*desc
    m = get_plans_matrix(data)
    orig_groups = dplyr::group_vars(data)
    dplyr::group_by(data, .data$draw) %>%
        dplyr::mutate(district = rank(ord * {{ x }})) %>%
        set_plan_matrix(`colnames<-`(renumber_matrix(m, .$district), colnames(m))) %>%
        dplyr::arrange(district, .by_group=TRUE) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(orig_groups)))
}



#' Helper function to check types for tidy wrappers
#' @noRd
check_tidy_types = function(map, .data) {
    if (!is.null(map) && !inherits(map, "redist_map"))
        cli_abort("{.arg map} must be a {.cls redist_map}")
    if (is.null(.data))
        cli_abort("Must provide {.arg .data} if not called within a {.pkg dplyr} verb")
    if (!inherits(.data, "redist_plans"))
        cli_abort("{.arg data} must be a {.cls redist_plans}")
}


#' @rdname redist.parity
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param ... passed on to \code{redist.parity}
#'
#' @concept analyze
#' @export
plan_parity <- function(map, .data = cur_plans(), ...) {
    check_tidy_types(map, .data)
    idxs = unique(as.integer(.data$draw))
    ndists = attr(map, "ndists")
    total_pop = map[[attr(map, "pop_col")]]
    if (is.null(total_pop)) cli_abort("Population vector missing from {.arg map}")

    rep(max_dev(get_plans_matrix(.data)[, idxs, drop=FALSE], total_pop, ndists),
        each = ndists)
}

#' @rdname redist.compactness
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param ... passed on to \code{redist.compactness}
#'
#' @concept analyze
#' @export
distr_compactness = function(map, measure="FracKept", .data=cur_plans(), ...) {
    check_tidy_types(map, .data)

    # districts not in ascending order
    if (length(unique(diff(as.integer(.data$district)))) > 2)
        cli_warn("Districts not sorted in ascending order; output may be incorrect.")

    idxs = unique(as.integer(.data$draw))
    redist.compactness(shp=map, plans=get_plans_matrix(.data)[, idxs, drop=FALSE],
                       measure=measure, total_pop=map[[attr(map, "pop_col")]],
                       adj=get_adj(map), ...)[[measure]]
}

#' @rdname redist.group.percent
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
group_frac = function(map, group_pop, total_pop=map[[attr(map, "pop_col")]],
                          .data=cur_plans()) {
    check_tidy_types(map, .data)
    # districts not in ascending order
    if (length(unique(diff(as.integer(.data$district)))) > 2)
        cli_warn("Districts not sorted in ascending order; output may be incorrect.")

    idxs = unique(as.integer(.data$draw))
    group_pop = rlang::eval_tidy(rlang::enquo(group_pop), map)
    total_pop = rlang::eval_tidy(rlang::enquo(total_pop), map)
    as.numeric(redist.group.percent(plans=get_plans_matrix(.data)[, idxs, drop=FALSE],
                                    group_pop=group_pop, total_pop=total_pop))
}

#' @rdname redist.segcalc
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
segregation_index = function(map, group_pop, total_pop=map[[attr(map, "pop_col")]],
                          .data=cur_plans()) {
    check_tidy_types(map, .data)
    idxs = unique(as.integer(.data$draw))
    group_pop = rlang::eval_tidy(rlang::enquo(group_pop), map)
    total_pop = rlang::eval_tidy(rlang::enquo(total_pop), map)
    plan_m = get_plans_matrix(.data)[, idxs, drop=FALSE]
    rep(as.numeric(redist.segcalc(plans=plan_m, group_pop=group_pop,
                                  total_pop=total_pop)),
        each=attr(map, "ndists"))
}

#' @rdname redist.metrics
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param ... passed on to \code{redist.metrics}
#'
#' @concept analyze
#' @export
partisan_metrics = function(map, measure, rvote, dvote, ...,
                            .data=cur_plans()) {
    check_tidy_types(map, .data)
    # districts not in ascending order
    if (length(unique(diff(as.integer(.data$district)))) > 2)
        cli_warn("Districts not sorted in ascending order; output may be incorrect.")

    idxs = unique(as.integer(.data$draw))
    rvote = rlang::eval_tidy(rlang::enquo(rvote), map)
    dvote = rlang::eval_tidy(rlang::enquo(dvote), map)
    as.numeric(redist.metrics(plans=get_plans_matrix(.data)[, idxs, drop=FALSE],
                              measure=measure, rvote=rvote, dvote=dvote, ...)[[measure]])
}

#' @rdname redist.competitiveness
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
competitiveness = function(map, rvote, dvote, .data=cur_plans()) {
    check_tidy_types(map, .data)
    idxs = unique(as.integer(.data$draw))
    rvote = rlang::eval_tidy(rlang::enquo(rvote), map)
    dvote = rlang::eval_tidy(rlang::enquo(dvote), map)
    rep(redist.competitiveness(plans=get_plans_matrix(.data)[, idxs, drop=FALSE],
                               rvote=rvote, dvote=dvote),
        each = attr(map, "ndists"))
}

#' @rdname redist.splits
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
county_splits = function(map, counties, .data=cur_plans()) {
    check_tidy_types(map, .data)
    idxs = unique(as.integer(.data$draw))
    counties = rlang::eval_tidy(rlang::enquo(counties), map)
    rep(redist.splits(plans=get_plans_matrix(.data)[, idxs, drop=FALSE], counties=counties),
        each = attr(map, "ndists"))
}


#' Extract the last plan from a set of plans
#'
#' @param plans A \code{\link{redist_plans}} object
#'
#' @returns An integer vector containing the final plan assignment.
#'
#' @concept analyze
#' @export
last_plan = function(plans) {
    plan_m = get_plans_matrix(plans)
    plan_m[, ncol(plan_m)]
}

#' Extract the district assignments for a precinct across all simulated plans
#'
#' @param prec the precinct number
#' @param .data a \code{\link{redist_plans}} object
#'
#' @return integer vector, a row from a plans matrix
#'
#' @concept analyze
#' @export
prec_assignment = function(prec, .data=cur_plans()) {
    check_tidy_types(NULL, .data)

    m = get_plans_matrix(.data)
    if (is.integer(prec)) {
        if (prec <= 0 || prec > nrow(m))
            cli_abort(c("{.arg prec} out of bounds",
                        "i"="There are {nrow(m)} precincts in these plans."))
    } else {
        cli_abort("{.arg prec} must be an integer index")
    }

    idxs = unique(as.integer(.data$draw))
    assignment = m[prec, idxs, drop=FALSE]
    if ("district" %in% colnames(.data) && is.factor(.data$district)) {
        lev = levels(.data$district)
        assignment = factor(lev[assignment], lev, ordered=is.ordered(.data$district))
    }

    assignment
}

#' Compute a matrix of precinct co-occurrences
#'
#' For a map with `n` precincts Returns an `n`-by-`n` matrix, where each
#' entry measures the fraction of the plans in which the row and column
#' precincts were in the same district.
#'
#' @param plans a [redist_plans] object.
#' @param which [`<data-masking>`][dplyr::dplyr_data_masking] which plans to
#'   compute the co-occurrence over.  Defaults to all.
#' @param sampled_only if `TRUE`, do not include reference plans.
#'
#' @return a symmetric matrix the size of the number of precincts.
#'
#' @concept analyze
#' @md
#' @export
prec_cooccurrence = function(plans, which=NULL, sampled_only=TRUE) {
    if (sampled_only)
        plans = subset_sampled(plans)
    which = eval_tidy(enquo(which), plans)
    plan_m = get_plans_matrix(plans)
    if (is.null(which))
        which = seq_len(ncol(plan_m))
    prec_cooccur(plan_m, which)
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
imp_confint = function(x, conf=0.95, .data=cur_plans()) {
    check_tidy_types(NULL, .data)

    idxs = unique(as.integer(.data$draw))
    y = rlang::eval_tidy(rlang::enquo(x), .data)
    ci = redist.smc_is_ci(y, get_plans_weights(.data)[, idxs, drop=FALSE], conf)

    tibble("{{ x }}" := mean(y),
                   "{{ x }}_lower" := ci[1],
                   "{{ x }}_upper" := ci[2])
}
