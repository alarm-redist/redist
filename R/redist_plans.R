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
#' @return matrix
#' @concept analyze
#' @export
get_plans_matrix = function(x) {
    stopifnot(inherits(x, "redist_plans"))
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
#'   to these weights.
#'
#' @concept analyze
#' @export
get_plans_weights = function(plans) {
    stopifnot(inherits(plans, "redist_plans"))
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
    stopifnot(inherits(x, "redist_plans"))
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

#' Extract the Metropolis Hastings Acceptance Rate
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns a numeric acceptance rate
#'
#' @concept simulate
#' @export
get_mh_acceptance_rate <- function(plans){
    stopifnot(inherits(plans, "redist_plans"))
    alg <- attr(plans, 'algorithm')

    if( alg %in% c('flip', 'mergesplit')){
        attr(plans, 'mh_acceptance')
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
    if ("district" %in% colnames(y)) {
        distrs = table(as.integer(y$district))
        ndists = max(plans_m[,1])
        if (any(distrs != distrs[1]) || length(distrs) != ndists)
            warning("Some districts may have been dropped. ",
                    "This will prevent summary statistics from working correctly.\n",
                    "To avoid this message, coerce using `as_tibble`.")
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
    constr = attr(objs[[1]], "constraints")
    resamp = attr(objs[[1]], "resampled")
    comp = attr(objs[[1]], "compactness")
    for (i in 2:n_obj) {
        if (nrow(get_plans_matrix(objs[[i]])) != n_prec)
            stop("Number of precincts must match for all sets of plans.")
        if (!identical(attr(objs[[i]], "prec_pop"), prec_pop))
            stop("Precinct populations must match for all sets of plans.")
        if (attr(objs[[i]], "resampled") != resamp)
            stop("Some sets of plans are resampled while others are not.")
        if (attr(objs[[i]], "compactness") != comp) {
            warning("Compactness values differ across sets of plans.")
            comp = NA
        }
        if (!identical(attr(objs[[i]], "constraints"), constr)) {
            warning("Constraints do not match for all sets of plans.")
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
    attr(ret, "plans") = do.call(cbind, lapply(objs, function(x) get_plans_matrix(x)))
    attr(ret, "wgt") = do.call(c, lapply(objs, function(x) get_plans_weights(x)))
    attr(ret, "n_eff") = sum(do.call(c, lapply(objs, function(x) attr(x, "n_eff"))))

    ret
}


#' Print method for redist_plans
#' @param x redist_plans object
#' @param \dots additional arguments
#' @method print redist_plans
#' @importFrom utils str
#' @return prints to console
#' @export
print.redist_plans = function(x, ...) {
    plans_m = get_plans_matrix(x)
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

    alg_name = c(mcmc="Flip Markov chain Monte Carlo",
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

    cat("Plans matrix:", utils::capture.output(str(plans_m, give.attr=FALSE)),
        "\n", sep="")

    utils::getS3method("print", "tbl")(x)

    invisible(x)
}
