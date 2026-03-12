#' Diagnose Single-Tree MMSS Partitions
#'
#' Samples plans from `map`, chooses one merged region `R` per sampled plan by
#' selecting `l` adjacent districts, then draws spanning trees on each induced
#' subgraph `G[R]` and counts the valid `l`-partitions of each tree.
#'
#' For each sampled plan, the diagnostic:
#' - draws one plan with [redist_smc()]
#' - selects `l` connected districts using the same random expansion rule used by
#'   MMSS
#' - forms the merged region `R` from those districts
#' - draws `n_trees_per_region` uniform spanning trees on `G[R]`
#' - counts the number `m(T)` of valid `l`-partitions of each tree
#' - records the reverse-proposal weight contribution
#'   `Z = 1[xi_ref in P(T)] / m(T)`
#'
#' A valid `l`-partition is a set of `l - 1` tree edges whose removal gives `l`
#' connected components, each with population inside
#' `[(1 - pop_tol) * pop(R) / l, (1 + pop_tol) * pop(R) / l]`.
#'
#' @param map A [redist_map] object.
#' @param l Number of districts to merge into each sampled region `R`.
#' @param n_plans Number of plans to draw with [redist_smc()]. Defaults to `50`.
#' @param n_trees_per_region Number of spanning trees to draw for each sampled
#'   merged region. Defaults to `200`.
#'
#' @return A list with elements:
#' - `hit_rate`: `Pr[m(T) >= 1]` across all sampled `(R, T)` pairs.
#' - `mean_m`: `E[m(T)]`.
#' - `max_m`: the maximum observed value of `m(T)`.
#' - `dist_m`: the empirical distribution of `m(T)` as a named integer vector.
#' - `cv2_Z`: `var(Z) / mean(Z)^2` for the sampled reference-plan reverse
#'   proposal weight.
#'
#' @export
diag_single_tree_partitions <- function(map, l, n_plans = 50, n_trees_per_region = 200) {
    map <- validate_redist_map(map)
    l <- as.integer(l)
    n_plans <- as.integer(n_plans)
    n_trees_per_region <- as.integer(n_trees_per_region)

    if (is.na(l) || l < 1L || l > attr(map, "ndists")) {
        cli::cli_abort("{.arg l} must be between 1 and the number of districts in {.arg map}.")
    }
    if (is.na(n_plans) || n_plans < 1L) {
        cli::cli_abort("{.arg n_plans} must be a positive integer.")
    }
    if (is.na(n_trees_per_region) || n_trees_per_region < 1L) {
        cli::cli_abort("{.arg n_trees_per_region} must be a positive integer.")
    }

    sampled <- redist_smc(
        map,
        nsims = n_plans,
        silent = TRUE,
        verbose = FALSE,
        ncores = 1,
        ref_name = FALSE
    )

    diag_single_tree_partitions_impl(
        get_adj(map),
        map[[attr(map, "pop_col")]],
        get_plans_matrix(sampled),
        attr(map, "ndists"),
        get_pop_tol(map),
        l,
        n_trees_per_region
    )
}
