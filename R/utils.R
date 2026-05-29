#' Load `redist` in parallel workers, respecting `pkgload::load_all()`.
#'
#' When `redist` is being used via `devtools::load_all()` (i.e., not installed
#' in the workers' library path, or installed at a different version), calling
#' `library(redist)` in a worker resolves to the wrong copy. This helper
#' detects a dev-loaded package and calls `pkgload::load_all()` on the workers
#' instead.
#'
#' @param cl a cluster object from `parallel::makeCluster()`
#' @noRd
init_workers <- function(cl) {
    dev <- requireNamespace("pkgload", quietly = TRUE) &&
        isTRUE(pkgload::is_dev_package("redist"))
    if (dev) {
        pkg_path <- getNamespaceInfo("redist", "path")
        parallel::clusterCall(
            cl,
            function(p) {
                pkgload::load_all(
                    p,
                    quiet = TRUE,
                    export_all = FALSE,
                    helpers = FALSE,
                    attach_testthat = FALSE
                )
            },
            pkg_path
        )
    } else {
        parallel::clusterEvalQ(cl, library(redist))
    }
    invisible(NULL)
}
