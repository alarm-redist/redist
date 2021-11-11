##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/11/07
## Purpose: redist constraints for a tidy workflow
##############################################


#######################
# constructors and reconstructors

# internal constructor
new_redist_constr = function(constr=list(), data=tibble()) {
    constr = reconstruct.redist_constr(constr, NULL)

    if (length(constr) > 0) {
        if (is.null(names(constr))) cli_abort("Null names.")
        if (any(names(constr) == "")) cli_abort("Empty names.")
        for (el in constr) {
            if (class(el) != "list") cli_abort("Not a nested list")
            classes = vapply(el, class, character(1))
            if (length(classes) == 0 || any(classes != "list"))
                cli_abort("Not a nested list")
        }
    }

    stopifnot(is.data.frame(data))
    attr(constr, "data") = data

    constr
}

validate_redist_constr = function(constr) {
    if (!is.list(constr)) stop("Not a list")
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
#' `r paste0("* [", ls("package:redist")[grep("add_constr_", ls("package:redist"))], "()]", collapse="\n")`
#'
#' More information about each constraint can be found on the relevant constraint page.
#'
#' @param map a [redist_map()] object; the map that will be used in sampling
#'
#' @returns a `redist_constr` object, which is just a list with a certain nested structure.
#'
#' @examples
#' constr = redist_constr()
#' constr = add_constr_splits(constr, strength=1.5)
#' print(constr)
#'
#' @md
#' @concept simulate
#' @export
redist_constr = function(map=tibble()) {
    new_redist_constr(data=map)
}

# properly nest
reconstruct.redist_constr = function(constr, old) {
    # handle deprecated constraint names
    if ("vra" %in% names(constr)) {
        .Deprecated("grp_pow", old="vra")
        constr = add_to_constr(constr, "grp_pow", constr$vra)
        constr$vra = NULL
    }
    if ("hinge" %in% names(constr)) {
        .Deprecated("grp_hinge", old="hinge")
        constr = add_to_constr(constr, "grp_hinge", constr$hinge)
        constr$hinge = NULL
    }

    for (i in seq_along(constr)) {
        classes = vapply(constr[[i]], class, character(1))
        if (any(classes != "list")) {
            constr[[i]] = list(constr[[i]])
        }
    }

    class(constr) = c("redist_constr", "list")
    constr
}

#######################
# constraint helpers

# helper
add_to_constr = function(constr, name, new_constr) {
    cur_constr = constr[[name]]

    if (!is.null(cur_constr)) {
        classes = vapply(cur_constr, class, character(1))
        if (any(classes != "list")) {
            constr[[name]] = list(cur_constr, new_constr)
        } else {
            constr[[name]] = c(cur_constr, list(new_constr))
        }
    } else {
        constr[[name]] = list(new_constr)
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
#'
#' The `grp_pow` constraint (for expert use) adds a term of the form
#' \eqn{(|tgtgroup-grouppct||tgtother-grouppct|)^{pow})}, which
#' encourages districts to have group shares near either `tgt_group`
#' or `tgt_other`.
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
#' The `custom` constraint allows the user to specify their own constraint using
#' a function which evaluates districts one at a time. The provided function
#' `fn` should take two arguments: a vector describing the current plan
#' assignment (which may be incomplete, in the case of SMC), and an integer
#' describing the district which to evaluate. An example is provided below.
#' The flexibility of this constraint comes with an additional computational
#' cost, since the other constraints are written in C++ and so are more
#' performant.
#'
#' @param constr A [redist_constr()] object
#' @param strength The strength of the constraint. Higher values mean a more restrictive constraint.
#'
#' @examples
#' data(iowa)
#' iowa_map = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.05)
#' constr = redist_constr(iowa_map)
#' constr = add_constr_splits(constr, strength=1.5)
#' constr = add_constr_grp_hinge(constr, strength=100,
#'                               dem_08, tot_08, tgts_group=c(0.5, 0.6))
#' # encourage districts to have the same number of counties
#' constr = add_constr_custom(constr, strength=1000, fn=function(plan, distr) {
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
add_constr_status_quo = function(constr, strength, current) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data = attr(constr, "data")
    if (missing(current)) current = get_existing(data)

    new_constr = list(strength=strength,
                      current=eval_tidy(enquo(current), data))
    if (is.null(current) || length(new_constr$current) != nrow(data))
        cli_abort("{.arg current} must be provided, and must have as many
                  precincts as the {.cls redist_map}")
    new_constr$n_current = max(new_constr$current)

    add_to_constr(constr, "status_quo", new_constr)
}

#' @param group_pop A vector of group population
#' @param total_pop A vector of total population. Defaults to the population vector used for sampling.
#' @param tgt_group,tgt_other Target group shares for the power-type constraint.
#' @param pow The exponent for the power-type constraint.
#'
#' @rdname constraints
#' @export
add_constr_grp_pow = function(constr, strength, group_pop, total_pop=NULL,
                              tgt_group=0.5, tgt_other=0.5, pow=1.0) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data = attr(constr, "data")

    new_constr = list(strength=strength,
                      group_pop=eval_tidy(enquo(group_pop), data),
                      total_pop=eval_tidy(enquo(total_pop), data),
                      tgt_group=tgt_group, tgt_other=tgt_other, pow=pow)
    if (is.null(new_constr$total_pop)) {
        if (!is.null(attr(data, "pop_col"))) {
            new_constr$total_pop = data[[attr(data, "pop_col")]]
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
add_constr_grp_hinge = function(constr, strength, group_pop, total_pop=NULL,
                                tgts_group=c(0.55)) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data = attr(constr, "data")

    new_constr = list(strength=strength,
                      group_pop=eval_tidy(enquo(group_pop), data),
                      total_pop=eval_tidy(enquo(total_pop), data),
                      tgts_group=tgts_group)
    if (is.null(new_constr$total_pop)) {
        if (!is.null(attr(data, "pop_col"))) {
            new_constr$total_pop = data[[attr(data, "pop_col")]]
        } else {
            cli_abort("{.arg total_pop} missing.")
        }
    }

    stopifnot(length(new_constr$group_pop) == nrow(data))
    stopifnot(length(new_constr$total_pop) == nrow(data))

    add_to_constr(constr, "grp_hinge", new_constr)
}

#' @param dvote,rvote A vector of Democratic or Republican vote counts
#' @rdname constraints
#' @export
add_constr_compet = function(constr, strength, dvote, rvote, pow=0.5) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data = attr(constr, "data")

    new_constr = list(strength=strength,
                      dvote=eval_tidy(enquo(dvote), data),
                      rvote=eval_tidy(enquo(rvote), data),
                      pow=pow)
    stopifnot(length(new_constr$dvote) == nrow(data))
    stopifnot(length(new_constr$rvote) == nrow(data))

    add_to_constr(constr, "compet", new_constr)
}

#' @param incumbents A vector of precinct indicess for incumbents
#' @rdname constraints
#' @export
add_constr_incumbency = function(constr, strength, incumbents) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")
    data = attr(constr, "data")

    new_constr = list(strength=strength,
                      incumbents=eval_tidy(enquo(incumbents), data))

    add_to_constr(constr, "incumbency", new_constr)
}

#' @rdname constraints
#' @export
add_constr_splits = function(constr, strength) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")

    new_constr = list(strength=strength)
    add_to_constr(constr, "splits", new_constr)
}

#' @rdname constraints
#' @export
add_constr_multisplits = function(constr, strength) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")

    new_constr = list(strength=strength)
    add_to_constr(constr, "multisplits", new_constr)
}

#' @param fn A function
#' @rdname constraints
#' @export
add_constr_custom = function(constr, strength, fn) {
    if (!inherits(constr, "redist_constr")) cli_abort("Not a {.cls redist_constr} object")
    if (strength <= 0) cli_warn("Nonpositive strength may lead to unexpected results")

    args = rlang::fn_fmls(fn)
    if (length(args) != 2) cli_abort("Function must take two arguments.")

    new_constr = list(strength=strength, fn=fn)
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
print.redist_constr = function(x, header=TRUE, details=TRUE, ...) {
    if (header)
        cli_text("A {.cls redist_constr} object with {length(x)} constraint{?s}")

    print_constr = function(x) {
        if (details) {
            idx_strength = which(names(x) == "strength")
            str(x[-idx_strength], no.list=T, comp.str="   ", give.attr=FALSE)
        }
    }

    x = unlist(x, recursive=FALSE)
    for (nm in names(x)) {
        if (startsWith(nm, "status_quo")) {
            cli::cli_bullets(c("*"="A status quo constraint of strength {x[[nm]]$strength}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "grp_pow")) {
            cli::cli_bullets(c("*"="A (power-type) group share constraint of strength {x[[nm]]$strength}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "grp_hinge")) {
            cli::cli_bullets(c("*"="A (hinge-type) group share constraint of strength {x[[nm]]$strength}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "compet")) {
            cli::cli_bullets(c("*"="A competitiveness constraint of strength {x[[nm]]$strength}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "incumbency")) {
            cli::cli_bullets(c("*"="An incumbency constraint of strength {x[[nm]]$strength}"))
            print_constr(x[[nm]])
        } else if (startsWith(nm, "splits")) {
            cli::cli_bullets(c("*"="A splits constraint of strength {x[[nm]]$strength}"))
        } else if (startsWith(nm, "multisplits")) {
            cli::cli_bullets(c("*"="A multisplits constraint of strength {x[[nm]]$strength}"))
        } else if (startsWith(nm, "custom")) {
            cli::cli_bullets(c("*"="A custom constraint of strength {x[[nm]]$strength}"))
            print_constr(x[[nm]])
        } else {
            cli::cli_bullets(c("*"="An unknown constraint {.var {nm}}"))
            print_constr(x[[nm]])
        }
    }
}

#' @method str redist_constr
#' @export
str.redist_constr = function(object, give.attr=FALSE, ...) {
    NextMethod("str", object, give.attr=give.attr, ...)
}

#' @method as.list redist_constr
#' @export
as.list.redist_constr = function(x, ...) {
    class(x) = "list"
    attr(x, "data") = NULL
    x
}
