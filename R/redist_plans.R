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
new_redist_plans = function(plans, map, algorithm, wgt, resampled=TRUE, ndists=attr(map, "ndists"), ...) {
    n_sims = ncol(plans)
    if (n_sims < 1) cli_abort("Need at least one simulation draw.")

    n_prec = nrow(plans)
    map_dists = attr(map, "ndists")
    partial = ndists < map_dists

    prec_pop = map[[attr(map, "pop_col")]]
    if (!partial) {
        distr_range = 1:ndists
        distr_pop = pop_tally(plans, prec_pop, ndists)
    } else {
        distr_range = 1:ndists - 1L
        distr_pop = pop_tally(plans, prec_pop, ndists)
        pl_tmp = plans + 1L
        distr_pop = pop_tally(pl_tmp, prec_pop, ndists)
    }

    attr_names = c("redist_attr", "plans", "ndists", "algorithm", "wgt",
                   "resampled", "ndists", "merge_idx", "prec_pop",
                   names(list(...)))

    structure(tibble(draw = rep(as.factor(1:n_sims), each=ndists),
                             district = rep(distr_range, n_sims),
                             total_pop = as.numeric(distr_pop)),
              plans=plans, ndists=ndists, algorithm=algorithm, wgt=wgt,
              resampled=resampled, merge_idx=attr(map, "merge_idx"),
              prec_pop=prec_pop, redist_attr=attr_names, ...,
              class=c("redist_plans", "tbl_df", "tbl", "data.frame"))
}

validate_redist_plans = function(x) {
    if (names(x)[1] != "draw") cli_abort("First column must be named \"{.field draw}\"")
    if (!is.factor(x$draw)) cli_abort("{.field draw} column must be a factor")

    plan_m = attr(x, "plans")
    if (is.null(plan_m)) cli_abort("Missing plans matrix")

    min_distr = colmin(plan_m)
    max_distr = colmax(plan_m)
    if (any(min_distr != 1) || any(diff(max_distr) != 0))
        cli_abort("District numbers must start at 1 and run sequentially to the number of districts.")

    x
}

reconstruct.redist_plans = function(data, old) {
    if (!("draw" %in% colnames(data)))
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

    if (inherits(data, "rowwise_df"))
        classes = c("rowwise_df", classes)

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
#' * \code{\link{is_county_split}}
#' * \code{\link{prec_assignment}}
#' * \code{\link{plan_distances}}
#' * \code{\link{get_plans_matrix}}
#' * \code{\link{get_plans_weights}}
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
#' iowa = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05, total_pop = pop)
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
    if (!is.matrix(plans)) cli_abort("{.arg plans} must be a matrix.")
    if (nrow(plans) != nrow(map)) cli_abort("{.arg plans} matrix must have as many rows as {.arg map} has precincts.")
    if (!inherits(map, "redist_map")) cli_abort("{.arg map} must be a {.cls redist_map}")

    if (min(plans) == 0L) plans = plans + 1L
    storage.mode(plans) = "integer"

    obj = new_redist_plans(plans, map, algorithm, wgt=wgt,
                           resampled=FALSE, ...)
    validate_redist_plans(obj)
}


# getters / setters ------------------------------------------------------------


#' Extract the matrix of district assignments from a redistricting simulation
#'
#' @param x the \code{redist_plans} object
#' @param ... ignored
#' @return matrix
#' @concept analyze
#' @export
get_plans_matrix = function(x) {
    if (!inherits(x, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    attr(x, "plans")
}
#' @rdname get_plans_matrix
#' @method as.matrix redist_plans
#' @return matrix
#' @export
as.matrix.redist_plans = function(x, ...) get_plans_matrix(x)

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
#'   to these weights. If weights have been resampled, this returns the weights
#'   before resampling (i.e., they do not correspond to the resampled plans).
#'
#' @concept analyze
#' @export
get_plans_weights = function(plans) {
    if (!inherits(plans, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    wgt = attr(plans, "wgt")
    if (!is.null(wgt))
        attr(wgt, "resampled") = attr(plans, "resampled")
    wgt
}

#' @rdname get_plans_weights
#' @param ... Ignored.
#' @importFrom stats weights
#' @method weights redist_plans
#' @return numeric vector
#' @export
weights.redist_plans = function(object, ...) {
    get_plans_weights(object)
}

get_n_ref = function(x) {
    if (!inherits(x, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    plans_m = get_plans_matrix(x)
    if (is.null(colnames(plans_m))) return(0)
    refs = which(nchar(colnames(plans_m)) > 0)
    length(unique(colnames(plans_m)[refs]))
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
    if (!inherits(plans, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
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
    plans_m = get_plans_matrix(plans)
    if (is.null(colnames(plans_m))) return(plans)
    refs = which(nchar(colnames(plans_m)) > 0)
    dplyr::filter(plans, !as.integer(.data$draw) %in% refs)
}

#' @rdname subset_sampled
#' @export
subset_ref = function(plans) {
    plans_m = get_plans_matrix(plans)
    if (is.null(colnames(plans_m))) return(plans)
    refs = which(nchar(colnames(plans_m)) > 0)
    dplyr::filter(plans, as.integer(.data$draw) %in% refs)
}

#' Extract the Metropolis Hastings Acceptance Rate
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns a numeric acceptance rate
#'
#' @concept simulate
#' @export
get_mh_acceptance_rate <- function(plans){
    if (!inherits(plans, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    alg <- attr(plans, "algorithm")

    if (alg %in% c("flip", "mergesplit")){
        attr(plans, "mh_acceptance")
    } else {
        NA_real_
    }
}

# generics ----------------------------------------------------------------


#' @method dplyr_row_slice redist_plans
#' @export
dplyr_row_slice.redist_plans = function(data, i, ...) {
    if (is.logical(i)) i = which(i)

    draws = rle(as.integer(data$draw))
    draws$values = seq_along(draws$values)
    draws_left = unique(inverse.rle(draws)[i])
    y = vctrs::vec_slice(data, i)
    plans_m = get_plans_matrix(data)

    # if we don't have every district present in every row
    # this check is necessary but not sufficient for what we want
    if (length(i) > 0 && "district" %in% colnames(y)) {
        distrs = table(as.integer(y$district))
        ndists = max(plans_m[,1])
        if (any(distrs != distrs[1]) || length(distrs) != ndists)
            cli_warn(c("Some districts may have been dropped. This will prevent summary statistics from working correctly.",
                    ">"="To avoid this message, coerce using {.fun as_tibble}."))
    }

    if (length(draws_left) != ncol(plans_m)) {
        attr(y, "wgt") = attr(y, "wgt")[draws_left]
        y = set_plan_matrix(y, plans_m[, draws_left, drop=FALSE])
    }

    if (is.factor(y$draw)) {
        y$draw <- droplevels(y$draw)
    }

    y
}

# 'template' is the old df
#' @method dplyr_reconstruct redist_plans
#' @export
dplyr_reconstruct.redist_plans = function(data, template) {
    reconstruct.redist_plans(data, template)
}

#' Combine multiple sets of redistricting plans
#'
#' Only works when all the sets are compatible---generated from the same map,
#' with the same number of districts.  Sets of plans will be indexed by the
#' `chain` column.
#'
#' @param ... The [`redist_plans`] objects to combine.  If named arguments are
#' provided, the names will be used in the `chain` column; otherwise, numbers
#' will be used for the `chain` column.
#' @param deparse.level Ignored.
#'
#' @return A new [`redist_plans`] object.
#'
#' @md
#' @concept analyze
#' @export
rbind.redist_plans = function(..., deparse.level=1) {
    objs = rlang::list2(...)
    n_obj = length(objs)
    if (n_obj == 1) return(objs[[1]])

    # check types
    n_prec = nrow(get_plans_matrix(objs[[1]]))
    prec_pop = attr(objs[[1]], "prec_pop")
    ndists = attr(objs[[1]], "ndists")
    constr = attr(objs[[1]], "constraints")
    resamp = attr(objs[[1]], "resampled")
    comp = attr(objs[[1]], "compactness")
    for (i in 2:n_obj) {
        if (nrow(get_plans_matrix(objs[[i]])) != n_prec)
            cli_abort("Number of precincts must match for all sets of plans.")
        if (!identical(attr(objs[[i]], "prec_pop"), prec_pop))
            cli_abort("Precinct populations must match for all sets of plans.")
        if (!identical(attr(objs[[i]], "ndists"), ndists))
            cli_abort("Number of districts must match for all sets of plans.")
        if (attr(objs[[i]], "resampled") != resamp)
            cli_abort("Some sets of plans are resampled while others are not.")
        if (!is.null(comp)) {
            if (attr(objs[[i]], "compactness") != comp) {
                cli_warn("Compactness values differ across sets of plans.")
                comp = NA
            }
        } else {
            if (!is.null(attr(objs[[i]], 'compactness'))) {
                cli_warn('Some compactness values were non-NULL. Set to {.val NA}.')
                comp <- NA
            }
        }
        if (!identical(attr(objs[[i]], "constraints"), constr)) {
            cli_warn("Constraints do not match for all sets of plans.")
            constr = NA
        }
    }

    ret = bind_rows(lapply(objs, dplyr::as_tibble), .id="chain")
    ret = reconstruct.redist_plans(ret, objs[[1]])
    attr(ret, "adapt_k_thresh") = NA
    attr(ret, "seq_alpha") = NA
    attr(ret, "pop_temper") = NA
    attr(ret, "compactness") = comp
    attr(ret, "constraints") = constr
    attr(ret, "ndists") = ndists
    attr(ret, "prec_pop") = prec_pop
    attr(ret, "plans") = do.call(cbind, lapply(objs, function(x) get_plans_matrix(x)))
    attr(ret, "wgt") = do.call(c, lapply(objs, function(x) get_plans_weights(x)))
    attr(ret, "n_eff") = sum(do.call(c, lapply(objs, function(x) attr(x, "n_eff"))))

    ret
}


#' Diagnostic information on sampled plans
#'
#' Currently only supported for [redist_smc()] output.
#'
#' @param object a [redist_plans] object
#' @param \dots additional arguments (ignored)
#'
#' @method summary redist_plans
#' @return A data frame containing diagnostic information, invisibly.
#'
#' @export
summary.redist_plans = function(object, ...) {
    algo = attr(object, "algorithm")
    name = deparse(substitute(object))
    if (algo == "smc") {
        diagn = attr(object, "diagnostics")
        plans_m = get_plans_matrix(object)
        n_ref = get_n_ref(object)
        n_samp = ncol(plans_m) - n_ref
        n_distr = attr(object, "ndists")

        fmt_comma = function(x) format(x, digits=0, big.mark=",")
        cli_text("{.strong SMC:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")

        est_div = plans_diversity(object)
        div_rg = format(quantile(est_div, c(0.1, 0.9)), digits=2)
        div_bad = (mean(est_div) <= 0.35) || (mean(est_div <= 0.05) > 0.2)
        div_warn = if (div_bad) " [WARNING: LOW DIVERSITY]" else ""
        cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}{.strong {div_warn}}")
        cat("\n")

        out = tibble(n_eff = c(diagn$step_n_eff, diagn$n_eff),
                     eff = c(diagn$step_n_eff, diagn$n_eff)/n_samp,
                     accept_rate = c(diagn$accept_rate, NA),
                     sd_log_wgt = diagn$sd_lp,
                     max_unique = diagn$unique_survive,
                     est_k = c(diagn$est_k, NA))

        tbl_print = as.data.frame(out)
        min_n = max(0.05*n_samp, min(0.4*n_samp, 100))
        bottlenecks = with(tbl_print, pmin(max_unique, n_eff) < min_n)
        tbl_print$bottleneck = ifelse(bottlenecks, "     *     ", "")
        tbl_print$n_eff = with(tbl_print,
                str_glue("{fmt_comma(n_eff)} ({sprintf('%0.1f%%', 100*eff)})"))
        tbl_print$eff = NULL
        tbl_print$accept_rate = with(tbl_print, sprintf('%0.1f%%', 100*accept_rate))
        max_pct = with(tbl_print, max_unique/(-n_samp * expm1(-1)))
        tbl_print$max_unique = with(tbl_print,
                str_glue("{fmt_comma(max_unique)} ({sprintf('%3.0f%%', 100*max_pct)})"))

        colnames(tbl_print) = c("Eff. samples (%)", "Accept. rate",
                                "Log weights s.d.", " Max. unique",
                                "Est. k", "Bottleneck?")
        rownames(tbl_print) = c(paste("Split", seq_len(n_distr-1)), "Final resample")

        print(tbl_print, digits=2)
        cat("\n")

        cli::cli_li(cli::col_grey("
            Watch out for low effective samples, very low acceptance rates (less than 1%),
            large std. devs. of the log weights (more than 10 or so),
            and low numbers of unique plans."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Check for potential bottlenecks.
                        Increase the number of samples.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance. If the accpetance rate drops
                        quickly in the final splits, try increasing
                        {.arg pop_temper} by 0.01.")
        }
        if (any(bottlenecks)) {
            cli::cli_li("{.strong Bottlenecks:} Consider weakening or removing
                        constraints, or increasing the population tolerance.
                        If the accpetance rate drops quickly in the final splits,
                        try increasing {.arg pop_temper} by 0.01.
                        To visualize what geographic areas may be causing problems,
                        try running the following code. Highlighted areas are
                        those that may be causing the bottleneck.")
            code = str_glue("plot(<map object>, colMeans(as.matrix({name}) %in% c({paste(which(bottlenecks), collapse=', ')})))")
            cli::cat_line("    ", cli::code_highlight(code, "Material"))
        }
    } else {
        cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
    }

    invisible(out)
}


#' Print method for \code{redist_plans}
#'
#' @param x a [redist_plans] object
#' @param \dots additional arguments (ignored)
#'
#' @method print redist_plans
#' @return The original object, invisibly.
#' @export
print.redist_plans = function(x, ...) {
    plans_m = get_plans_matrix(x)
    n_ref = get_n_ref(x)
    n_samp = ncol(plans_m) - n_ref
    nd = attr(x, "ndists")
    if (is.null(nd)) nd = max(plans_m[,1])

    if (n_ref > 0)
        cli_text("A {.cls redist_plans} containing {n_samp} sampled plan{?s} and
                 {n_ref} reference plan{?s}")
    else
        cli_text("A {.cls redist_plans} containing {n_samp} sampled plan{?s}")

    if (ncol(plans_m) == 0) return(invisible(x))

    alg_name = c(mcmc="Flip Markov chain Monte Carlo",
                 smc="Sequential Monte Carlo",
                 mergesplit="Merge-split Markov chain Monte Carlo",
                 rsg="random seed-and-grow",
                 crsg="compact random seed-and-grow",
                 enumpart="Enumpart",
                 shortburst="short bursts")[attr(x, "algorithm")]
    if (is.na(alg_name)) alg_name = "an unknown algorithm"

    cli_text("Plans have {nd} districts from a {nrow(plans_m)}-unit map,
             and were drawn using {alg_name}.")

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

    cat("Plans matrix:", utils::capture.output(str(plans_m, give.attr=FALSE)),
        "\n", sep="")

    utils::getS3method("print", "tbl")(x)

    invisible(x)
}
