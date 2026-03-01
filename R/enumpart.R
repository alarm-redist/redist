#' Enumerate Partitions
#'
#' Enumerates all valid partitions (or a uniform random subset) for a given
#' redistricting problem.
#'
#' Path arguments default to `tempdir()`, so no file management is required
#' for typical use. Supply explicit paths to retain intermediate files across
#' sessions. This is recommended
#'
#' If a file already exists at `out_path`, `enumpart` is not run and
#' the existing file is read directly. This makes it easy to resume or
#' re-read a previous enumeration.
#'
#' @param map A [redist_map] object
#' @param ndists The number of districts. Defaults to the number specified
#'   in `map`
#' @param n Number of partitions to draw uniformly at random. If `NULL`
#'   (default), all valid partitions are enumerated
#' @param total_pop A vector of precinct populations. Defaults to the
#'   population column in `map`
#' @param lower Lower bound on each district's total population. Defaults to
#'   the lower bound in `map` if one is set
#' @param upper Upper bound on each district's total population. Defaults to
#'   the upper bound in `map` if one is set
#' @param ordered_path Path (without `".dat"`) for the ordered edge file.
#'   Defaults to a path in `tempdir()`
#' @param out_path Path (without `".dat"`) for the enumerated output file.
#'   Defaults to a path in `tempdir()`. If the file already exists, it is read
#'    directly
#' @param weight_path Path (without `".dat"`) for the vertex weights file.
#'   Defaults to a path in `tempdir()` when `lower` or `upper` are non-`NULL`.
#' @param read If `TRUE` (default), reads and returns the enumerated plans as a
#'   [redist_plans] object. If `FALSE`, runs enumpart and returns `out_path`
#'   invisibly.
#' @param init If `TRUE`, builds enumpart before running. Must be
#'   run once after installing the package. Defaults to `FALSE`.
#'
#' @return A [redist_plans] object, or `out_path` invisibly if `read = FALSE`.
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @concept enumerate
#' @export
#'
#' @examples \dontrun{
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#' plans <- redist_enumpart(fl_map)
#'
#' # Write to disk without reading (useful for large enumerations)
#' path <- redist_enumpart(fl_map, read = FALSE)
#' plans <- redist_enumpart(fl_map, out_path = path)  # reads existing file
#' }
redist_enumpart <- function(map,
                            ndists = attr(map, "ndists"),
                            n = NULL,
                            total_pop = map[[attr(map, "pop_col")]],
                            lower = attr(map, "pop_bounds")[1L],
                            upper = attr(map, "pop_bounds")[3L],
                            ordered_path = tempfile("enum_ord"),
                            out_path = tempfile("enum_out"),
                            weight_path = if (!is.null(lower) || !is.null(upper)) tempfile("enum_wgt") else NULL,
                            read = TRUE,
                            init = FALSE) {
    map <- validate_redist_map(map)
    adj <- get_adj(map)
    n_prec <- nrow(map)

    if (!file.exists(paste0(out_path, ".dat"))) {
        if (init) {
            if (Sys.info()[["sysname"]] == "Windows") {
                makecontent <- readLines(system.file("enumpart/Makefile", package = "redist"))
                makecontent[7] <- "\tg++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG -lpsapi"
                writeLines(text = makecontent, con = system.file("enumpart/Makefile", package = "redist"))
            }
            sys::exec_wait("make", args = c("-C", system.file("enumpart", package = "redist")))
            if (Sys.info()[["sysname"]] == "Windows") {
                makecontent <- readLines(system.file("enumpart/Makefile", package = "redist"))
                makecontent[7] <- "\tg++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG"
                writeLines(text = makecontent, con = system.file("enumpart/Makefile", package = "redist"))
            }
        }

        rlang::check_installed("igraph")
        ordered_mat <- do.call(rbind, ndscut(lapply(adj, unique)))
        utils::write.table(
            data.frame(ordered_mat),
            file = paste0(ordered_path, ".dat"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
        )

        if (!is.null(weight_path)) {
            utils::write.table(
                t(total_pop),
                file = paste0(weight_path, ".dat"),
                quote = FALSE, row.names = FALSE, col.names = FALSE
            )
        }

        if (is.null(n)) {
            ep_args <- c("-k", as.integer(ndists), "-comp", "-allsols")
        } else {
            ep_args <- c("-k", as.integer(ndists), "-comp", "-sample", as.integer(n))
        }
        if (!is.null(lower)) {
            ep_args <- c(ep_args, "-lower", as.character(lower))
        }
        if (!is.null(upper)) {
            ep_args <- c(ep_args, "-upper", as.character(upper))
        }
        ep_args <- c(
            paste0(ordered_path, ".dat"),
            if (!is.null(weight_path)) paste0(weight_path, ".dat"),
            ep_args
        )

        sys::exec_wait(
            paste0(system.file("enumpart", package = "redist"), "/enumpart"),
            args = ep_args,
            std_out = paste0(out_path, ".dat"),
            std_err = TRUE
        )
    }

    if (!read) {
        return(invisible(out_path))
    }

    vals <- scan(paste0(out_path, ".dat"), integer(), quiet = TRUE)
    plans <- matrix(vals, nrow = n_prec) + 1L
    n_plans <- ncol(plans)

    new_redist_plans(
        plans, map,
        "enumpart",
        wgt = rep(1, n_plans), resampled = FALSE, ndists = ndists
    )
}


#' Calculate Graph Frontier Size
#'
#' Computes the frontier sizes along an edge ordering, which determines the
#' computational complexity of the enumpart algorithm. Smaller frontier sizes
#' lead to faster enumeration. This is primarily a helper for `enumpart`.
#'
#' @param adj An adjacency list (as returned by [redist.adjacency()]). Ignored
#'   if `ordered_path` is provided.
#' @param ordered_path Path (without `".dat"` extension) to an ordered edge
#'   file written by [redist_enumpart()]. If provided, `adj` is ignored.
#'
#' @return A list with four elements:
#' - `max`: maximum frontier size
#' - `average`: average frontier size
#' - `average_sq`: average of squared frontier sizes
#' - `sequence`: numeric vector of frontier sizes at each step
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @concept enumerate
#' @export
#'
#' @examples
#' data(fl25)
#' adj <- redist.adjacency(fl25)
#' redist_enumpart_frontier(adj)
redist_enumpart_frontier <- function(adj = NULL, ordered_path = NULL) {
    if (!is.null(ordered_path)) {
        vals <- scan(paste0(ordered_path, ".dat"), integer(), quiet = TRUE)
        edges_unsort <- matrix(vals, ncol = 2, byrow = TRUE)
    } else if (!is.null(adj)) {
        edges_unsort <- do.call(rbind, lapply(seq_along(adj), function(i) {
            nbrs <- adj[[i]] + 1L  # 0-indexed -> 1-indexed
            nbrs <- nbrs[nbrs > i]
            if (length(nbrs) == 0L) return(NULL)
            cbind(i, nbrs)
        }))
    } else {
        cli::cli_abort("Must provide either {.arg adj} or {.arg ordered_path}.")
    }

    edges <- cbind(pmin(edges_unsort[, 1], edges_unsort[, 2]),
                   pmax(edges_unsort[, 1], edges_unsort[, 2]))

    n <- nrow(edges)
    frontier <- rep(FALSE, max(edges))
    frontier_sizes <- numeric(n + 1L)

    for (i in seq_len(n)) {
        e1 <- edges[i, 1]
        e2 <- edges[i, 2]
        frontier[e1] <- TRUE
        frontier[e2] <- TRUE

        if (i == n || !any(edges[(i + 1):n, ] == e1)) frontier[e1] <- FALSE
        if (i == n || !any(edges[(i + 1):n, ] == e2)) frontier[e2] <- FALSE

        frontier_sizes[i + 1L] <- sum(frontier)
    }

    list(
        max = max(frontier_sizes),
        average = mean(frontier_sizes),
        average_sq = mean(frontier_sizes^2),
        sequence = frontier_sizes
    )
}


# deprecations ----

#' Initialize enumpart
#'
#' @description
#' **Deprecated**: Use [redist_enumpart()] with `init = TRUE` instead.
#'
#' @returns 0 on success.
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @export
#' @concept enumerate
#' @keywords internal
redist.init.enumpart <- function() {
    .Deprecated("redist_enumpart")
    if (Sys.info()[["sysname"]] == "Windows") {
        makecontent <- readLines(system.file("enumpart/Makefile", package = "redist"))
        makecontent[7] <- "\tg++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG -lpsapi"
        writeLines(text = makecontent, con = system.file("enumpart/Makefile", package = "redist"))
    }
    sys::exec_wait("make", args = c("-C", system.file("enumpart", package = "redist")))
    if (Sys.info()[["sysname"]] == "Windows") {
        makecontent <- readLines(system.file("enumpart/Makefile", package = "redist"))
        makecontent[7] <- "\tg++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG"
        writeLines(text = makecontent, con = system.file("enumpart/Makefile", package = "redist"))
    }
    0
}


#' Prepare enumpart edge ordering
#'
#' @description
#' **Deprecated**: Use [redist_enumpart()] instead, which handles preparation
#' automatically.
#'
#' @param adj Zero-indexed adjacency list.
#' @param ordered_path Valid path to output the ordered adjacency map to.
#' @param weight_path A path (not including `".dat"`) to store vertex weights.
#'   Only supply with `total_pop`.
#' @param total_pop The vector of precinct populations. Only supply with
#'   `weight_path`.
#' @param unordered_path Deprecated, ignored.
#'
#' @returns 0 on success.
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @export
#' @concept enumerate
#' @keywords internal
redist.prep.enumpart <- function(adj, ordered_path, weight_path = NULL,
                                 total_pop = NULL, unordered_path) {
    .Deprecated("redist_enumpart")
    if (!missing(unordered_path)) {
        cli_warn("{.arg unordered_path} is deprecated and will be ignored.")
    }
    if (is.null(weight_path) + is.null(total_pop) == 1L) {
        cli::cli_abort("You must provide both of {.arg weight_path} and {.arg total_pop} or neither.")
    }
    rlang::check_installed("igraph")
    ordered_mat <- do.call(rbind, ndscut(lapply(adj, unique)))
    utils::write.table(data.frame(ordered_mat),
                       file = paste0(ordered_path, ".dat"),
                       quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    if (!is.null(weight_path)) {
        utils::write.table(t(total_pop),
                           file = paste0(weight_path, ".dat"),
                           quote = FALSE, row.names = FALSE, col.names = FALSE
        )
    }
    0L
}


#' Run the enumpart algorithm
#'
#' @description
#' **Deprecated**: Use [redist_enumpart()] instead.
#'
#' @param ordered_path Path used in `redist.prep.enumpart()` (not including
#'   `".dat"`).
#' @param out_path Valid path to output the enumerated districts.
#' @param ndists Number of districts to enumerate.
#' @param all If `TRUE`, outputs all districts. If `FALSE`, draws `n`
#'   districts uniformly at random.
#' @param n Number of districts to output if `all` is `FALSE`.
#' @param weight_path A path (not including `".dat"`) to a vertex weights file.
#' @param lower A lower bound on each partition's total weight.
#' @param upper An upper bound on each partition's total weight.
#' @param options Additional enumpart arguments. Not recommended.
#'
#' @returns 0 on success.
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @export
#' @concept enumerate
#' @keywords internal
redist.run.enumpart <- function(ordered_path, out_path, ndists = 2,
                                all = TRUE, n = NULL, weight_path = NULL,
                                lower = NULL, upper = NULL, options = NULL) {
    .Deprecated("redist_enumpart")
    if (is.null(options)) {
        if (all) {
            options <- c("-k", as.integer(ndists), "-comp", "-allsols")
        } else {
            if (is.null(n)) {
                cli::cli_abort("{.arg n} must be specified when {.code all = FALSE}.")
            }
            options <- c("-k", as.integer(ndists), "-comp", "-sample", as.integer(n))
        }
    }
    if (!is.null(lower)) options <- c(options, "-lower", as.character(lower))
    if (!is.null(upper)) options <- c(options, "-upper", as.character(upper))
    options <- c(paste0(ordered_path, ".dat"),
                 if (!is.null(weight_path)) paste0(weight_path, ".dat"),
                 options)
    sys::exec_wait(paste0(system.file("enumpart", package = "redist"), "/enumpart"),
                   args = options,
                   std_out = paste0(out_path, ".dat"), std_err = TRUE
    )
}


#' Read Results from enumpart
#'
#' @description
#' **Deprecated**: Use [redist_enumpart()] instead, which reads results
#' automatically.
#'
#' @param out_path `out_path` specified in `redist.run.enumpart()`.
#' @param skip Number of lines to skip.
#' @param n_max Maximum number of lines to read.
#'
#' @returns A district membership matrix (1-indexed).
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @export
#' @concept enumerate
#' @keywords internal
redist.read.enumpart <- function(out_path, skip = 0, n_max = -1L) {
    .Deprecated("redist_enumpart")
    sols <- readLines(paste0(out_path, ".dat"), n = n_max)
    if (skip > 0) sols <- sols[-seq_len(skip)]
    apply(do.call("cbind", strsplit(sols, " ")), 2, as.integer) + 1L
}


#' Calculate Frontier Size
#'
#' @description
#' **Deprecated**: Use [redist_enumpart_frontier()] instead.
#'
#' @param ordered_path Path to ordered edge file created by
#'   `redist.prep.enumpart()`.
#'
#' @returns A list with `max`, `average`, `average_sq`, and `sequence`.
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @export
#' @concept enumerate
#' @keywords internal
redist.calc.frontier.size <- function(ordered_path) {
    .Deprecated("redist_enumpart_frontier")
    redist_enumpart_frontier(ordered_path = ordered_path)
}


#' Enumerate All Partitions
#'
#' @description
#' **Deprecated**: Use [redist_enumpart()] instead.
#'
#' Single function for standard enumeration analysis, using ZDD methodology
#' (Fifield, Imai, Kawahara, and Kenny 2020).
#'
#' @param adj Zero-indexed adjacency list.
#' @param ordered_path Valid path to output the ordered adjacency map to.
#' @param out_path Valid path to output the enumerated districts.
#' @param ndists Number of districts to enumerate.
#' @param all If `TRUE`, outputs all districts. If `FALSE`, draws `n`
#'   districts uniformly at random.
#' @param n Number of districts to output if `all` is `FALSE`.
#' @param weight_path A path (not including `".dat"`) to a vertex weights file.
#' @param lower A lower bound on each partition's total weight.
#' @param upper An upper bound on each partition's total weight.
#' @param init If `TRUE`, runs initialization first.
#' @param read If `TRUE`, reads and returns the results.
#' @param total_pop The vector of precinct populations.
#' @param unordered_path Deprecated, ignored.
#'
#' @returns A list with `plans` (district membership matrix) and `parity`.
#'
#' @references
#' Fifield, B., Imai, K., Kawahara, J., & Kenny, C. T. (2020).
#' The essential role of empirical validation in legislative redistricting
#' simulation. *Statistics and Public Policy*, 7(1), 52-68.
#'
#' @concept enumerate
#' @export
#' @keywords internal
redist.enumpart <- function(adj, ordered_path, out_path, ndists = 2,
                            all = TRUE, n = NULL, weight_path = NULL,
                            lower = NULL, upper = NULL, init = FALSE,
                            read = TRUE, total_pop = NULL, unordered_path) {
    .Deprecated("redist_enumpart")
    if (!missing(unordered_path)) {
        cli_warn("{.arg unordered_path} is deprecated and will be ignored.")
    }
    if (init) suppressWarnings(redist.init.enumpart())
    suppressWarnings(redist.prep.enumpart(adj, ordered_path, weight_path, total_pop))
    suppressWarnings(redist.run.enumpart(ordered_path, out_path, ndists, all, n,
                                         weight_path, lower, upper))
    if (!read) return(0)
    cds <- suppressWarnings(redist.read.enumpart(out_path))
    par <- if (!is.null(total_pop)) redist.parity(cds, total_pop) else rep(NA_real_, ncol(cds))
    list(plans = cds, parity = par)
}
