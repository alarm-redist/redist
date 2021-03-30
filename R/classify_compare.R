#' Make a comparison between two sets of plans
#'
#' This function provides one way to identify the structural differences between
#' two sets of redistricting plans. It operates by computing the precinct
#' co-occurrence matrix (a symmetric matrix where the i,j-th entry is the
#' fraction of plans where precinct i and j are in the same district) for each
#' set, and then computing the first eigenvalue of the difference in these two
#' matrices (in each direction). These eigenvalues identify the important parts
#' of the map.
#'
#' The co-occurence matrices are regularized with a \eqn{Beta(1/ndists, 1-1/ndists)}
#' prior, which is useful for when either `set1` or `set2` is small.
#'
#' @param plans a [redist_plans] object
#' @param set1,set2 [`<data-masking>`][dplyr::dplyr_data_masking] indexing vectors
#'   for the plan draws to compare. Must be mutually exclusive.
#' @param shp a shapefile for plotting.
#' @param plot If `plot="line"`, display a plot for each set showing the set of
#'   boundaries which most distinguish it from the other set (the squared
#'   differences in the eigenvector values across the boundary).  If
#'   `plot="fill"`, plot the eigenvector for each set as a choropleth. See below
#'   for more information.  Set to `FALSE` to disable plotting (or leave out
#'   `shp`).
#' @param thresh the value to threshold the eigenvector at in determining the
#'   relevant set of precincts for comparison.
#' @param labs the names of the panels in the plot.
#'
#' @return If possible, makes a comparison plot according to `plot`. Otherwise
#'   returns the following list:
#' \item{eigen1}{A numeric vector containing the first eigenvector of
#'   \code{p1 - p2}, where \code{p1} and \code{p2} are the co-occurrence matrices
#'   for \code{set1} and \code{set2}, respectively.}
#' \item{eigen2}{A numeric vector containing the first eigenvector of
#'   \code{p2 - p1}, where \code{p1} and \code{p2} are the co-occurrence matrices
#'   for \code{set1} and \code{set2}, respectively.}
#' \item{group_1a, group_1b}{Lists of precincts. Compared to \code{set2}, in the
#'   \code{set1} plans these precincts were much more likely to be in separate
#'   districts. Computed by thresholding \code{eigen1} at \code{thresh}.}
#' \item{group_2a, group_2b}{Lists of precincts. Compared to \code{set1}, in the
#'   \code{set2} plans these precincts were much more likely to be in separate
#'   districts. Computed by thresholding \code{eigen2} at \code{thresh}.}
#' \item{cooccur_sep_1}{The difference in the average co-occurrence of precincts
#'   in \code{group_1a} and \code{group_1b} between \code{set2} and \code{set1}.
#'   Higher indicates better separation.}
#' \item{cooccur_sep_2}{The difference in the average co-occurrence of precincts
#'   in \code{group_2a} and \code{group_2b} between \code{set1} and \code{set2}.
#'   Higher indicates better separation.}
#'
#'
#' @md
#' @concept analyze
#' @export
compare_plans = function(plans, set1, set2, shp=NULL, plot="line", thresh=0.1,
                         labs=c("Set 1", "Set 2")) {
    stopifnot(inherits(plans, "redist_plans"))

    set1 = eval_tidy(enquo(set1), plans)
    set2 = eval_tidy(enquo(set2), plans)
    if (is.logical(set1)) set1 = which(set1)
    if (is.logical(set2)) set2 = which(set2)
    if (length(intersect(set1, set2)) > 0)
        stop("`set1` and `set2` must be mutually exclusive.")
    n1 = length(set1)
    n2 = length(set2)

    pm = get_plans_matrix(plans)
    base_co = 1 / max(pm[, 1]) # baseline coccurence

    p1 = (n1*prec_cooccur(pm, set1) + base_co) / (n1 + 1)
    p2 = (n2*prec_cooccur(pm, set2) + base_co) / (n2 + 1)

    evecs = eigen(p1 - p2, symmetric=TRUE)$vectors
    evec1 = evecs[, 1]
    evec2 = evecs[, nrow(pm)]

    group_1a = which(evec1 >= thresh)
    group_1b = which(evec1 <= -thresh)
    group_2a = which(evec2 >= thresh)
    group_2b = which(evec2 <= -thresh)

    out = list(eigen1 = evec1,
               eigen2 = evec2,
               group_1a = group_1a,
               group_1b = group_1b,
               group_2a = group_2a,
               group_2b = group_2b,
               cooccur_sep_1 = mean(p2[group_1a, group_1b]) -
               mean(p1[group_1a, group_1b]),
               cooccur_sep_2 = mean(p1[group_2a, group_2b]) -
               mean(p2[group_2a, group_2b]))

    if (inherits(shp, "sf")) {
        if (plot == "line") {
            edges = dplyr::select(shp, attr(shp, "sf_column")) %>%
                sf::st_intersection() %>%
                dplyr::as_tibble() %>%
                dplyr::filter(.data$n.overlaps == 2) %>%
                dplyr::mutate(from = sapply(.data$origins, function(x) x[1]),
                              to = sapply(.data$origins, function(x) x[2]),
                              wgt1 = (evec1[.data$from] - evec1[.data$to])^2,
                              wgt2 = (evec2[.data$from] - evec2[.data$to])^2) %>%
                dplyr::filter(sf::st_dimension(.data$geometry) == 1) %>%
                sf::st_as_sf()

            make_plot = function(x, lab) {
                ggplot(edges, aes(size=x)) +
                    geom_sf() +
                    ggplot2::guides(size=F) +
                    ggplot2::scale_size_continuous(range=c(0, 3)) +
                    labs(title=lab) +
                    theme_void()
            }
            p1 = make_plot(edges$wgt1, labs[1])
            p2 = make_plot(edges$wgt2, labs[2])

            p1 + p2 + patchwork::plot_annotation(title="Distinguishing edges")
        } else if (plot == "fill") {
            make_plot = function(x, lab) {
                ggplot(shp, aes(fill=x)) +
                    geom_sf(size=0) +
                    ggplot2::guides(fill=F) +
                    labs(title=lab) +
                    theme_void()
            }
            p1 = make_plot(evec1, labs[1])
            p2 = make_plot(evec2, labs[2])

            p1 + p2 + patchwork::plot_annotation(title="Eigenvectors")
        }
    } else {
        out
    }
}

make_classif_lbl = function(idxs) {
    n = length(idxs)
    out = character(n)
    opts = list(c("I", "II", "III", "IV", "V", "VI", "VII", "VIII"),
                c("A", "B", "C", "D", "E", "F", "G", "H"),
                c("1", "2", "3", "4", "5", "6", "7", "8"),
                c("a", "b", "c", "d", "e", "f", "g", "h"),
                c("i", "ii", "iii", "iv", "v", "vi", "vii", "viii"))
    n_opts = length(opts)
    for (i in seq_len(n)) out[i] = opts[[i]][idxs[i]]
    paste0(out, collapse=".")
}

#' Hierarchically classify a set of redistricting plans
#'
#' Applies hierarchical clustering to a distance matrix computed from a set of
#' plans and takes the first `k` splits.
#'
#' @param dist_mat a distance matrix, the output of [plan_distances()]
#' @param k the number of groupings to create
#' @param method the clustering method to use. See [hclust()] for options.
#'
#' @return An object of class `redist_classified`, which is a list with two
#'   elements:
#' \item{groups}{A character vector of group labels of the form `"I.A.1.a.i"`,
#' one for each plan.}
#' \item{splits}{A list of splits in the hierarchical clustering. Each list
#' element is a list of two mutually exclusive vectors of plan indices, labeled
#' by their group classification, indicating the plans on each side of the split.}
#' Use [plot.redist_classified()] for a visual summary.
#'
#' @concept analyze
#' @md
#' @export
classify_plans = function(dist_mat, k=8, method="complete") {
    stopifnot(is.matrix(dist_mat))
    stopifnot(isSymmetric(dist_mat))
    stopifnot(k >= 2)

    cl = stats::hclust(stats::as.dist(dist_mat), method)
    tr = stats::as.dendrogram(cl)
    cut_height = utils::tail(cl$height, k)[1]

    queue = list(1L, 2L)
    groups = character(nrow(dist_mat))
    splits = list(list(I=labels(tr[[1]]), II=labels(tr[[2]])))
    while (length(queue) > 0) {
        node_idx = queue[[1]]
        node = tr[[node_idx]]
        queue = queue[-1]

        if (attr(node, "height") > cut_height) {
            left_idx = c(node_idx, 1)
            right_idx = c(node_idx, 2)
            split_obj = list()
            split_obj[[make_classif_lbl(left_idx)]] = labels(tr[[left_idx]])
            split_obj[[make_classif_lbl(right_idx)]] = labels(tr[[right_idx]])
            splits = c(splits, list(split_obj))

            queue = c(queue, list(left_idx, right_idx))
        } else {
            groups[labels(node)] = make_classif_lbl(node_idx)
        }
    }

    out = list(groups = groups,
               splits = splits)
    class(out) = "redist_classified"
    out
}


#' @export
print.redist_classified = function(x, ...) {
    n_split = length(x$splits)
    cat(length(x$groups), "plans classified into", n_split+1L, "groups.\n")
    cat("Group assignment:", utils::capture.output(str(x$group, vec.len=3)),
        "\n", sep="")
    for (i in seq_len(n_split)) {
        split = x$splits[[i]]
        cat("Split ", i, ":\n", sep="")
        cat("    ", names(split)[2], ": ",
            utils::capture.output(str(split[[1]], vec.len=3)), "\n", sep="")
        cat("    ", names(split)[2], ": ",
            utils::capture.output(str(split[[2]], vec.len=3)), "\n", sep="")
    }
}


#' Plot a plan classification
#'
#' @param x a `redist_classified` object, the output of [classify_plans()].
#' @param plans a [redist_plans] object.
#' @param shp a shapefile or [redist_map] object.
#' @param type either `"line"` or `"fill"`. Passed on to [compare_plans()] as
#'   `plot`.
#' @param ... passed on to [compare_plans()]
#'
#' @concept analyze
#' @md
#' @export
plot.redist_classified = function(x, plans, shp, type="line", ...) {
    stopifnot(inherits(plans, "redist_plans"))
    stopifnot(inherits(shp, "sf"))

    plots = lapply(x$splits, function(split) {
        compare_plans(plans, split[[1]], split[[2]], shp, plot=type,
                      ..., labs=names(split))
    })
    patchwork::wrap_plots(plots, ncol=1)
}