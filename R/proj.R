#' Calculate Projective Distributions, Averages, and Contrasts for a Summary Statistic
#'
#' @description
#' The *projective distribution* of a district-level summary statistic (McCartan
#' 2024) is the distribution of values of that statistic across a set of plans
#' for the district each precinct belongs to. The *projective average* of a
#' statistic is the average value of the projective distribution in each
#' precinct. A *projective contrast* is the difference between the projective
#' average for a single plan and the projective average for an ensemble of
#' sampled plans.
#'
#' It is very important to properly account for variation in the projective
#' distribution when looking at projective contrasts. The `pfdr` argument to
#' `proj_contr()` will calculate q-values for each precinct that can be used to
#' control the positive false discovery rate (pFDR) to avoid being misled by
#' this variation. See [redist.plot.contr_pfdr()] for a way to automatically
#' plot projective contrasts with this false discovery rate control.
#'
#' @param plans A [redist_plans] object.
#' @param x A district-level summary statistic calculated from the `plans`
#'   object. Tidy-evaluated in `plans`.
#' @param draws which draws/samples to include in the projective distribution.
#'   `NULL` will include all draws, including reference plans. The special value
#'   `NA` will include all sampled (non-reference) draws. An integer, logical,
#'   or character vector indicating specific draws may also be provided.
#'
#' @references
#' McCartan, C. (2024). Projective Averages for Summarizing Redistricting
#' Ensembles. *arXiv preprint*. Available at <https://arxiv.org/abs/>.
#'
#' @examples
#' data(iowa)
#' map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
#' plans <- redist_smc(map, 50, silent = TRUE)
#' plans$dem <- group_frac(map, dem_08, tot_08, plans)
#'
#' proj_distr(plans, dem) # a 99-by-50 matrix
#' plot(map, proj_avg(plans, dem))
#' plot(map, proj_contr(plans, dem))
#' plot(map, proj_contr(plans, dem, comp="cd_2010"))
#'
#' @concept analyze
#' @name proj
NULL


#' @returns `proj_distr`: A matrix with a row for each precinct (row in the map
#'   object) and a column for every draw described by `draws`.
#'
#' @rdname proj
#' @export
proj_distr <- function(plans, x, draws=NA) {
    plans_m <- get_plans_matrix(plans)

    n_ref <- 0
    # copied from get_n_ref()
    if (!is.null(colnames(plans_m))) {
        n_ref <- sum(nchar(colnames(plans_m)) > 0)
    }

    if (is.null(draws)) {
        draw_idx <- seq_len(ncol(plans_m))
    } else if (length(draws) == 1 && is.na(draws)) {
        if (n_ref > 0) {
            draw_idx <- seq_len(ncol(plans_m))[-seq_len(n_ref)]
        } else {
            draw_idx <- seq_len(ncol(plans_m))
        }
    } else if (is.logical(draws)) {
        draw_idx <- which(draws)
    } else {
        draw_idx <- match(as.character(draws), levels(plans$draw))
    }

    plans <- vctrs::vec_slice(plans, vctrs::vec_order(cbind(plans$draw, plans$district)))
    n_distr <- attr(plans, "ndists")
    x <- rlang::eval_tidy(rlang::enquo(x), plans)

    proj_distr_m(plans_m, x, draw_idx, n_distr)
}

#' @returns `proj_avg`: A numeric vector of length matching the number of
#'   precincts.
#'
#' @rdname proj
#' @export
proj_avg <- function(plans, x, draws=NA) {
    rowMeans(proj_distr(plans, {{ x }}, draws))
}

#' @param compare The plan to compare to the rest of the ensemble (which is
#'   controlled by `draws`). Defaults to the first reference plan, if any exists
#' @param norm If `TRUE`, normalize the contrast by the standard deviation of
#'   the projective distribution, precinct-wise. This will make the projective
#'   contrast in terms of z-scores.
#' @param pfdr If `TRUE`, calculate q-values for each precinct that can be
#'   used to control the positive false discovery rate (pFDR) at a given level
#'   by thresholding the q-values at that level. Q-values are stored as the
#'   `"q"` attribute on the returned vector. Requires the `matrixStats` package
#'   be installed.
#'
#' @returns `proj_contr`: A numeric vector of length matching the number of
#'   precincts, optionally with a `"q"` attribute containing q-values.
#'
#' @rdname proj
#' @export
proj_contr <- function(plans, x, compare=NA, draws=NA, norm=FALSE, pfdr=FALSE) {
    plans_m <- get_plans_matrix(plans)

    n_ref <- 0
    # copied from get_n_ref()
    if (!is.null(colnames(plans_m))) {
        n_ref <- sum(nchar(colnames(plans_m)) > 0)
    }

    if (length(compare) != 1) {
        cli_abort("{.arg compare} must refer to a single plan")
    }
    comp_idx <- if (is.na(compare)) {
        if (n_ref > 0) {
            1
        } else {
            cli_abort("If there are no reference plans, {.arg comp} must be specified and not NA.")
        }
    } else if (is.logical(compare)) {
        which(compare)
    } else {
        match(as.character(draws), levels(plans$draw))
    }

    draw_idx <- if (is.null(draws)) {
         seq_len(ncol(plans_m))
    } else if (length(draws) == 1 && is.na(draws)) {
        # VS ABOVE: also remove comp_idx
        if (n_ref > 0) {
            seq_len(ncol(plans_m))[-c(comp_idx, seq_len(n_ref))]
        } else {
            seq_len(ncol(plans_m))[-comp_idx]
        }
    } else if (is.logical(draws)) {
        which(draws)
    } else {
        match(as.character(draws), levels(plans$draw))
    }

    plans <- vctrs::vec_slice(plans, vctrs::vec_order(cbind(plans$draw, plans$district)))
    nd <- attr(plans, "ndists")
    x <- rlang::eval_tidy(rlang::enquo(x), plans)
    pd = proj_distr_m(plans_m, x, draw_idx, nd)
    pa = rowMeans(pd)

    contr = x[(comp_idx - 1)*nd + 1:nd][plans_m[, comp_idx]] - pa

    if (isTRUE(norm)) {
        sds = sqrt(rowMeans((pd - rowMeans(pd))^2))
        contr = contr / sds
    }

    if (isTRUE(pfdr)) {
        rlang::check_installed("matrixStats", "for pFDR control.")
        pd = pd - pa
        m_pv = matrixStats::rowRanks(abs(pd), ties.method="max") / ncol(pd)
        pv_draw = matrixStats::rowMeans2(abs(pd) >= abs(contr)) + 1 / ncol(pd)

        attr(contr, "q") = est_pfdr(pv_draw, m_pv) |>
            qvalues(pv_draw)
    }

    contr
}

# Algorithm 1 of Storey & Tibshirani (2001)
est_pfdr <- function(p, ref, lambda=0.5, n_alpha=15) {
    alphas = quantile(p, probs=seq(0, 1, length.out=n_alpha), type=1, names=FALSE)
    m = length(p)
    pi0 = (m - sum(p <= lambda)) / (m - mean(colSums(ref <= lambda)))
    pi0 = min(pi0, 1)
    fdrs = vapply(alphas, function(a) {
        pi0 *  mean(colSums(ref <= a)) / sum(p <= a)
    }, 0.0)
    ok = !is.nan(fdrs)

    list(
        alpha = alphas,
        pfdr = dplyr::coalesce(fdrs, 0.0)
    )
}
# Produce q values with spline
qvalues <- function(ests, p) {
    idx = vctrs::vec_unique_loc(ests$alpha)
    splinefun(ests$alpha[idx], rev(cummin(rev(ests$pfdr[idx]))), method="monoH.FC")(p)
}




#' Plot a Projective Contrast with positive False Discovery Rate (pFDR) Control
#'
#' Plot a projective contrast on a map with areas selected by the pFDR control
#' procedure hatched.
#'
#' @param map A [redist_map] object
#' @param contr The output of [proj_contr()] with `pfdr=TRUE`: A vector containing
#'   the contrast and an attribute `"q"` containing the q-values.
#' @param level The positive false discovery rate level to control.
#' @param density The density of the hatching (roughly what portion is shaded).
#' @param spacing The spacing of the hatches.
#'
#' @returns A ggplot.
#'
#' @examples
#' # example code
#' set.seed(1812)
#' data(iowa)
#' map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
#' plans <- redist_smc(map, 50, silent = TRUE)
#' plans$dem <- group_frac(map, dem_08, tot_08, plans)
#'
#' pc = proj_contr(plans, dem, pfdr=TRUE)
#' redist.plot.contr_pfdr(map, pc, level=0.4) # high `level` just to demonstrate
#'
#'
#' @concept plot
#' @export
redist.plot.contr_pfdr <- function(map, contr, level=0.05, density=0.2, spacing=0.015) {
    if (is.null(attr(contr, "q"))) {
        cli_abort("Must provide {.arg pfdr=TRUE} to {.fn proj_contr} to use {.fn redist.plot.contr_pfdr}.")
    }

    p = plot(map, contr)
    q = attr(contr, "q")
    if (!is.numeric(level) || level < 0 || level > 1) {
        cli_abort("{.arg level} must be a number between 0 and 1.")
    }
    if (!any(q <= level)) return(p)
    rlang::check_installed("ggpattern", "for pFDR-controlling projective contrast plots.")

    geom_sel = map |>
        dplyr::filter(q <= level) |>
        dplyr::summarize(is_coverage=TRUE) |>
        suppressWarnings()

    p +
        ggplot2::geom_sf(data=geom_sel, inherit.aes=FALSE,
                         fill=NA, linewidth=0.15, color="black") +
        ggpattern::geom_sf_pattern(data=geom_sel, inherit.aes=FALSE,
                                   fill=NA, linewidth=0.0, color=NA,
                                   pattern_fill="#00000070", pattern_color=NA,
                                   pattern_density=density, pattern_spacing=spacing)
}
